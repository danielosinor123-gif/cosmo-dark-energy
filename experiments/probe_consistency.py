"""probe_consistency.py — probe-by-probe consistency certificates for w(z).

Idea (the novel test): reconstruct w(z) in 3 redshift bins independently from different
probe subsets and quantify where they disagree. Disagreement localizes systematics —
this is what decides whether the 3-sigma evolving-DE hint is real.

Model: flat, w piecewise-constant in bins [0,0.75), [0.75,1.5), [1.5,inf); params
[Om, Ob, h, rd, w0, w1, w2]. Data: DESI DR2 BAO (all or per-tracer) + DES-Dovekie SN +
Planck-2018 compressed (R, wb).

Outputs: results/probe_results.json — per combo: w_i posterior mean/std, chi2_map, and
pairwise tension matrix T_ij = (mu_i^A - mu_i^B)/sqrt(sA^2 + sB^2) per bin, plus w(z) at
pivot redshifts with 1-sigma bands.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import BackgroundBinned, bao_chi2, sn_hd_chi2, cmb_compressed_chi2, metropolis

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")
os.makedirs(RESULTS, exist_ok=True)

EDGES = [0.0, 0.75, 1.5, np.inf]
NAMES = ["Om", "Ob", "h", "rd", "w0", "w1", "w2"]
PRIORS = {"Om": (0.05, 0.9), "Ob": (0.005, 0.10), "h": (0.55, 0.85), "rd": (120.0, 165.0),
          "w0": (-3.0, 1.0), "w1": (-3.0, 1.0), "w2": (-3.0, 1.0)}
X0 = {"Om": 0.31, "Ob": 0.022 / 0.67 ** 2, "h": 0.67, "rd": 147.0,
      "w0": -0.9, "w1": -0.8, "w2": -0.6}

print("loading data ...", flush=True)
bao_all = cosmo.load_desi_bao(which="ALL")[0]
bao_tracers = cosmo.load_desi_bao(which="tracers")
z_sn, mu_sn, invC_sn = cosmo.load_dovekie()


def make_logpost(combo):
    bao_rows = combo.get("bao")          # None, "ALL", or a tracer dict
    use_sn = combo.get("sn", False)
    use_cmb = combo.get("cmb", False)
    if bao_rows == "ALL":
        bao_rows = bao_all

    def logpost(x):
        d = dict(zip(NAMES, x))
        for n in NAMES:
            lo, hi = PRIORS[n]
            if not (lo <= d[n] <= hi):
                return -np.inf
        Om, Ob, h, rd = d["Om"], d["Ob"], d["h"], d["rd"]
        wb = [d["w0"], d["w1"], d["w2"]]
        H0 = 100.0 * h
        chi2 = 0.0
        if bao_rows is not None:
            chi2 += bao_chi2(BackgroundBinned(Om, wb, EDGES), H0, rd, bao_rows)
        if use_sn:
            chi2 += sn_hd_chi2(BackgroundBinned(Om, wb, EDGES), H0, z_sn, mu_sn, invC_sn)
        if use_cmb:
            chi2 += cmb_compressed_chi2(Om, Ob, h, BackgroundBinned(Om, wb, EDGES))
        return -0.5 * chi2

    return logpost


def w_of_z(z, w_bins):
    """piecewise-constant w at redshift z."""
    for i, wi in enumerate(w_bins):
        a = EDGES[i]
        b = EDGES[i + 1] if i + 1 < len(EDGES) else np.inf
        if a <= z < b:
            return wi
    return w_bins[-1]


def run(tag, combo, nsteps=24000, burn=6000, seed=3):
    logpost = make_logpost(combo)
    x0 = np.array([X0[n] for n in NAMES])
    step = np.array([0.01, 0.001, 0.005, 0.5, 0.05, 0.05, 0.05])
    t0 = time.time()
    samples, acc = metropolis(logpost, x0, nsteps=nsteps, burn=burn, step=step, seed=seed)
    stride = max(1, len(samples) // 3000)
    sub = samples[::stride]
    lps = np.array([logpost(s) for s in sub])
    chi2_map = float(-2.0 * lps.max())
    wi = samples[:, 4:7]
    piv = {}
    for zp in [0.3, 0.6, 1.0, 1.4, 2.0]:
        vals = np.array([w_of_z(zp, s) for s in wi])
        piv[f"z={zp}"] = [float(vals.mean()), float(vals.std())]
    summ = {
        "tag": tag, "combo": {k: (v if isinstance(v, str) else bool(v)) for k, v in combo.items() if k != "bao"} if isinstance(combo, dict) else combo,
        "acceptance": float(acc), "chi2_map": chi2_map,
        "w_bins_mean": wi.mean(axis=0).tolist(), "w_bins_std": wi.std(axis=0).tolist(),
        "w_pivots": piv, "nsamples": int(len(samples)), "seconds": round(time.time() - t0, 1),
        "params_mean": dict(zip(NAMES, map(float, samples.mean(axis=0)))),
    }
    np.save(os.path.join(RESULTS, f"chain_binned_{tag}.npy"), np.column_stack([sub, lps]))
    print(f"  [{tag}] chi2_map={chi2_map:.2f} acc={acc:.2f} w_bins="
          + " ".join(f"{m:.2f}+-{s:.2f}" for m, s in zip(wi.mean(axis=0), wi.std(axis=0)))
          + f" ({summ['seconds']}s)", flush=True)
    return summ


def main():
    combos = [
        ("joint_all", dict(bao="ALL", sn=True, cmb=True)),
        ("bao_only_cmb", dict(bao="ALL", sn=False, cmb=True)),
        ("sn_only_cmb", dict(bao=None, sn=True, cmb=True)),
        ("bao_sn_nocmb", dict(bao="ALL", sn=True, cmb=False)),
    ]
    for tr in bao_tracers:
        short = tr["name"].replace("_GCcomb", "").replace("BGS_BRIGHT-21.35", "BGS") \
                          .replace("LOPnotqso", "").replace("LRG+ELG", "LRGELG")
        combos.append((f"tracer_{short}", dict(bao=tr, sn=True, cmb=True)))

    out = []
    for tag, combo in combos:
        out.append(run(tag, combo))
        with open(os.path.join(RESULTS, "probe_results.json"), "w") as f:
            json.dump(out, f, indent=2)

    # pairwise tension matrix per w-bin between probe subsets (vs joint)
    joint = out[0]
    print("\n=== pairwise tension vs joint (per w-bin) ===", flush=True)
    for r in out[1:]:
        T = [(a - b) / np.hypot(sa, sb) for a, b, sa, sb in zip(
            r["w_bins_mean"], joint["w_bins_mean"], r["w_bins_std"], joint["w_bins_std"])]
        r["tension_vs_joint"] = T
        print(f"  {r['tag']}: " + " ".join(f"{t:+.2f}" for t in T), flush=True)
    with open(os.path.join(RESULTS, "probe_results.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("PROBE DONE", flush=True)


if __name__ == "__main__":
    main()
