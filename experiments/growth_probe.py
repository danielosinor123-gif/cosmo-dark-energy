"""growth_probe.py — growth-vs-geometry test: add public f sigma8 (RSD/PV) data to the
BAO+SN+CMB pipeline and see whether the evolving-DE (w0wa) hint survives an independent
growth probe. sigma8_0 is a free nuisance parameter.

Data: data/fs8 (public f sigma8 benchmark catalogue, 24 points, published WiggleZ/eBOSS
covariance blocks) — see src/cosmo.py load_fs8.

Models: lcdm/wcdm/w0wa + s8.  Combos: bao_sn_cmb_fs8, bao_sn_fs8, bao_cmb_fs8, bao_fs8, fs8.
Usage: python experiments/growth_probe.py
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import (Background, bao_chi2, sn_hd_chi2, cmb_compressed_chi2, fs8_chi2,
                   metropolis)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")
os.makedirs(RESULTS, exist_ok=True)

print("loading data ...", flush=True)
bao_all = cosmo.load_desi_bao(which="ALL")[0]
z_sn, mu_sn, invC_sn = cosmo.load_dovekie()
z_f, fs8v, invC_f = cosmo.load_fs8()
print(f"  BAO {len(bao_all['z'])} rows, SN {len(z_sn)}, fs8 {len(z_f)} pts", flush=True)

PRIORS = {
    "Om": (0.05, 0.9), "Ob": (0.005, 0.10), "h": (0.55, 0.85), "rd": (120.0, 165.0),
    "s8": (0.2, 1.4),
    "w": (-1.6, -0.4), "w0": (-1.8, -0.2), "wa": (-3.0, 3.0),
}
MODELS = {
    "lcdm": ["Om", "Ob", "h", "rd", "s8"],
    "wcdm": ["Om", "Ob", "h", "rd", "s8", "w"],
    "w0wa": ["Om", "Ob", "h", "rd", "s8", "w0", "wa"],
}
X0 = {"Om": 0.31, "Ob": 0.049, "h": 0.67, "rd": 147.0, "s8": 0.81, "w": -1.0, "w0": -0.9, "wa": -0.3}


def make_logpost(model, combo):
    names = MODELS[model]
    use_bao = "bao" in combo
    use_sn = "sn" in combo
    use_cmb = "cmb" in combo
    use_fs8 = "fs8" in combo

    def logpost(x):
        d = dict(zip(names, x))
        for n in names:
            lo, hi = PRIORS[n]
            if not (lo <= d[n] <= hi):
                return -np.inf
        Om, Ob, h, rd, s8 = d["Om"], d["Ob"], d["h"], d["rd"], d["s8"]
        w0 = d.get("w0", -1.0)
        wa = d.get("wa", 0.0)
        if "w" in d:
            w0, wa = d["w"], 0.0
        H0 = 100.0 * h
        bg = Background(Om, w0, wa)
        chi2 = 0.0
        if use_bao:
            chi2 += bao_chi2(bg, H0, rd, bao_all)
        if use_sn:
            chi2 += sn_hd_chi2(bg, H0, z_sn, mu_sn, invC_sn)
        if use_cmb:
            chi2 += cmb_compressed_chi2(Om, Ob, h, bg)
        if use_fs8:
            chi2 += fs8_chi2(Om, s8, z_f, fs8v, invC_f, w0, wa)
        return -0.5 * chi2

    return logpost, names


def sigma_from_dchi2(dchi2, k):
    from scipy.stats import chi2 as chi2dist, norm
    p = chi2dist.sf(abs(dchi2), k)
    z = float(norm.isf(p)) if p < 0.5 else 0.0
    return z * float(np.sign(dchi2)) if dchi2 != 0 else 0.0, float(p)


def run(model, combo, nsteps=18000, burn=5000, seed=9):
    logpost, names = make_logpost(model, combo)
    x0 = np.array([X0[n] for n in names])
    step = np.array([0.01 if n not in ("rd", "s8") else 0.5 if n == "rd" else 0.02 for n in names])
    t0 = time.time()
    samples, acc = metropolis(logpost, x0, nsteps=nsteps, burn=burn, step=step, seed=seed)
    stride = max(1, len(samples) // 2500)
    sub = samples[::stride]
    lps = np.array([logpost(s) for s in sub])
    chi2_map = float(-2.0 * lps.max())
    summ = {"model": model, "combo": combo, "acceptance": float(acc), "chi2_map": chi2_map,
            "mean": dict(zip(names, map(float, samples.mean(axis=0)))),
            "std": dict(zip(names, map(float, samples.std(axis=0)))),
            "seconds": round(time.time() - t0, 1)}
    np.save(os.path.join(RESULTS, f"chain_fs8_{model}_{combo}_s{seed}.npy"),
            np.column_stack([sub, lps]))
    print(f"  [{combo}] {model}: chi2_map={chi2_map:.2f} acc={acc:.2f} ({summ['seconds']}s) "
          + " ".join(f"{n}={m:.3f}+-{s:.3f}" for n, m, s in
                     zip(names, samples.mean(axis=0), samples.std(axis=0))), flush=True)
    return summ


def main():
    combos = sys.argv[1:] if len(sys.argv) > 1 else \
        ["bao_sn_cmb_fs8", "bao_sn_fs8", "bao_cmb_fs8", "bao_fs8", "fs8"]
    try:
        with open(os.path.join(RESULTS, "growth_results.json")) as f:
            out = [e for e in json.load(f) if e["combo"] not in combos]
    except Exception:
        out = []
    for combo in combos:
        res = {}
        for model in ["lcdm", "wcdm", "w0wa"]:
            res[model] = run(model, combo)
        for a, b, k in [("lcdm", "wcdm", 1), ("lcdm", "w0wa", 2), ("wcdm", "w0wa", 1)]:
            d = res[b]["chi2_map"] - res[a]["chi2_map"]
            sig, p = sigma_from_dchi2(d, k)
            res[f"{b}_vs_{a}"] = dict(dchi2=float(d), sigma=sig, p=float(p))
            print(f"  >> {combo}: {b} vs {a}: dchi2={d:+.2f} -> {abs(sig):.2f} sigma "
                  f"({'prefers ' + b if d < 0 else 'prefers ' + a})", flush=True)
        out.append({"combo": combo, "runs": res})
        with open(os.path.join(RESULTS, "growth_results.json"), "w") as f:
            json.dump(out, f, indent=2)
    print("GROWTH DONE", flush=True)


if __name__ == "__main__":
    main()
