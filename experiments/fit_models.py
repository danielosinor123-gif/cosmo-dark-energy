"""fit_models.py — independent re-analysis of the DESI DR2 evolving-dark-energy hint.

Pipeline (numpy/scipy only, our own code in src/cosmo.py):
  data   : DESI DR2 BAO (all tracers) + DES-Dovekie SN Ia + Planck 2018 compressed CMB prior
  models : flat LCDM (Om, Ob, h, rd) | flat wCDM (+w) | flat w0waCDM (+w0, wa)
  sampler: adaptive Metropolis (src/cosmo.py)

Usage:
  python experiments/fit_models.py main      # 3 models x main combo + leave-one-out combos
  python experiments/fit_models.py tracers   # per-tracer BAO consistency (w0wa vs LCDM)
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import Background, bao_chi2, sn_hd_chi2, sn_mB_chi2, cmb_compressed_chi2, metropolis

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")
os.makedirs(RESULTS, exist_ok=True)

# ----------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------
print("loading data ...", flush=True)
bao_all = cosmo.load_desi_bao(which="ALL")[0]
bao_tracers = cosmo.load_desi_bao(which="tracers")
z_sn, mu_sn, invC_sn = cosmo.load_dovekie()
print(f"  BAO ALL: {len(bao_all['z'])} rows; tracers: {[t['name'] for t in bao_tracers]}")
print(f"  DES-Dovekie SN: {len(z_sn)} SNe", flush=True)

# ----------------------------------------------------------------------------
# Parameterisations:  theta = [Om, Ob, h, rd] (+ w) (+ w0, wa)
# ----------------------------------------------------------------------------
PRIORS = {
    "Om": (0.05, 0.9), "Ob": (0.005, 0.10), "h": (0.55, 0.85), "rd": (120.0, 165.0),
    "w": (-1.6, -0.4), "w0": (-1.8, -0.2), "wa": (-3.0, 3.0),
}
MODELS = {
    "lcdm": ["Om", "Ob", "h", "rd"],
    "wcdm": ["Om", "Ob", "h", "rd", "w"],
    "w0wa": ["Om", "Ob", "h", "rd", "w0", "wa"],
}
X0 = {"Om": 0.31, "Ob": 0.022 / 0.67 ** 2, "h": 0.67, "rd": 147.0, "w": -1.0, "w0": -0.9, "wa": 0.3}


def make_logpost(model, combo, bao_rows=None, nsteps_hint=1):
    names = MODELS[model]
    bao_rows = bao_rows if bao_rows is not None else bao_all
    use_bao = "bao" in combo
    use_sn = "sn" in combo
    use_cmb = "cmb" in combo

    def logpost(x):
        d = dict(zip(names, x))
        for n in names:
            lo, hi = PRIORS[n]
            if not (lo <= d[n] <= hi):
                return -np.inf
        Om, Ob, h, rd = d["Om"], d["Ob"], d["h"], d["rd"]
        w0 = d.get("w0", -1.0)
        wa = d.get("wa", 0.0)
        if "w" in d:
            w0, wa = d["w"], 0.0
        H0 = 100.0 * h
        chi2 = 0.0
        if use_bao:
            chi2 += bao_chi2(Background(Om, w0, wa), H0, rd, bao_rows)
        if use_sn:
            chi2 += sn_hd_chi2(Background(Om, w0, wa), H0, z_sn, mu_sn, invC_sn)
        if use_cmb:
            chi2 += cmb_compressed_chi2(Om, Ob, h, Background(Om, w0, wa))
        return -0.5 * chi2

    return logpost, names


def run(model, combo, nsteps=24000, burn=6000, seed=0, bao_rows=None, tag=""):
    logpost, names = make_logpost(model, combo, bao_rows)
    x0 = np.array([X0[n] for n in names])
    step = np.array([0.01 if n != "rd" else 0.5 for n in names])
    t0 = time.time()
    samples, acc = metropolis(logpost, x0, nsteps=nsteps, burn=burn, step=step, seed=seed)
    # best-fit = max logpost over a subsample (kept consistent for saving)
    stride = max(1, len(samples) // 2000)
    sub = samples[::stride]
    lps = np.array([logpost(s) for s in sub])
    i_best = int(np.argmax(lps))
    best = sub[i_best]
    chi2_map = -2.0 * lps[i_best]
    summ = {
        "model": model, "combo": combo, "tag": tag, "acceptance": float(acc),
        "chi2_map": float(chi2_map), "best": dict(zip(names, map(float, best))),
        "mean": dict(zip(names, map(float, samples.mean(axis=0)))),
        "std": dict(zip(names, map(float, samples.std(axis=0)))),
        "nsamples": int(len(samples)), "seconds": round(time.time() - t0, 1),
    }
    np.save(os.path.join(RESULTS, f"chain_{model}_{combo}{('_' + tag) if tag else ''}_s{seed}.npy"),
            np.column_stack([sub, lps]))
    print(f"  [{combo}{'/' + tag if tag else ''}] {model}: chi2_map={chi2_map:.2f} acc={acc:.2f} "
          f"({summ['seconds']}s) " + " ".join(f"{n}={m:.3f}+-{s:.3f}" for n, m, s in
          zip(names, samples.mean(axis=0), samples.std(axis=0))), flush=True)
    return summ


def sigma_from_dchi2(dchi2, k):
    """Equivalent Gaussian sigma for a chi2 difference with k extra params.
    dchi2 < 0 means the more complex model is preferred; report |dchi2| significance."""
    from scipy.stats import chi2 as chi2dist, norm
    p = chi2dist.sf(abs(dchi2), k)
    return float(norm.isf(p)) * np.sign(dchi2), float(p)


def main_mode():
    combos = ["bao_sn_cmb", "bao_sn", "bao_cmb", "bao", "sn_cmb"]
    out = []
    for combo in combos:
        res = {}
        for model in ["lcdm", "wcdm", "w0wa"]:
            res[model] = run(model, combo, nsteps=24000, seed=1)
        for a, b, k in [("lcdm", "wcdm", 1), ("lcdm", "w0wa", 2), ("wcdm", "w0wa", 1)]:
            d = res[b]["chi2_map"] - res[a]["chi2_map"]
            sig, p = sigma_from_dchi2(d, k)
            res[f"{b}_vs_{a}"] = dict(dchi2=float(d), sigma=sig, p=float(p),
                                      note=("prefers " + b) if d < 0 else ("prefers " + a))
            print(f"  >> {combo}: {b} vs {a}: dchi2={d:+.2f} -> {abs(sig):.2f} sigma "
                  f"({'prefers ' + b if d < 0 else 'prefers ' + a}, p={p:.2e})", flush=True)
        out.append({"combo": combo, "runs": res})
        with open(os.path.join(RESULTS, "main_results.json"), "w") as f:
            json.dump(out, f, indent=2)
    print("MAIN DONE", flush=True)


def tracers_mode():
    out = []
    for tr in bao_tracers:
        res = {}
        for model in ["lcdm", "w0wa"]:
            res[model] = run(model, "bao_sn_cmb", nsteps=18000, seed=2, bao_rows=tr, tag=tr["name"])
        d = res["w0wa"]["chi2_map"] - res["lcdm"]["chi2_map"]
        sig, p = sigma_from_dchi2(d, 2)
        res["w0wa_vs_lcdm"] = dict(dchi2=float(d), sigma=sig, p=float(p))
        print(f"  >> tracer {tr['name']}: dchi2={d:+.2f} -> {abs(sig):.2f} sigma "
              f"({'prefers w0wa' if d < 0 else 'prefers lcdm'})", flush=True)
        out.append({"tracer": tr["name"], "runs": res})
        with open(os.path.join(RESULTS, "tracer_results.json"), "w") as f:
            json.dump(out, f, indent=2)
    print("TRACERS DONE", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "main"
    if mode == "main":
        main_mode()
    elif mode == "tracers":
        tracers_mode()
