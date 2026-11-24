"""session9.py — CAPSTONE: the full systematics-inclusive joint fit.

w0wa vs LCDM, both with ALL nuisance/systematic params free at once:
  - dM: free magnitude step for SNe with z<0.1   (low-z calibration, §16)
  - Cepheid/SH0ES H0 prior 73.04 +- 1.04          (anchor, §12)
  - s8: free, against the 24-pt f sigma8 growth   (growth, §7)
on BAO + Dovekie + CMB-compressed. Question: with every known systematic allowed,
is wa still < 0?  -> results/capstone_results.json

Usage: python experiments/session9.py capstone
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import Background, bao_chi2, sn_hd_chi2, cmb_compressed_chi2, fs8_chi2, metropolis

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

print("loading data ...", flush=True)
bao_all = cosmo.load_desi_bao(which="ALL")[0]
z_sn, mu_sn, invC_sn = cosmo.load_dovekie()
z_f, fs8v, invC_f = cosmo.load_fs8()
mask = (z_sn < 0.1).astype(float)

PRIORS = {"Om": (0.05, 0.9), "Ob": (0.005, 0.10), "h": (0.55, 0.85), "rd": (120.0, 165.0),
          "s8": (0.2, 1.4), "w0": (-1.8, -0.2), "wa": (-3.0, 3.0), "dm": (-0.3, 0.3)}
X0 = {"Om": 0.31, "Ob": 0.049, "h": 0.70, "rd": 144.0, "s8": 0.80, "w0": -0.9, "wa": -0.3,
      "dm": 0.02}
H0_CEPHEID = (73.04, 1.04)


def sigma_from_dchi2(dchi2, k):
    from scipy.stats import chi2 as chi2dist, norm
    p = chi2dist.sf(abs(dchi2), k)
    z = float(norm.isf(p)) if p < 0.5 else 0.0
    return z * float(np.sign(dchi2)) if dchi2 != 0 else 0.0, float(p)


def make_logpost(model):
    names = ["Om", "Ob", "h", "rd", "dm", "s8"] + (["w0", "wa"] if model == "w0wa" else [])

    def logpost(x):
        d = dict(zip(names, x))
        for n in names:
            lo, hi = PRIORS[n]
            if not (lo <= d[n] <= hi):
                return -np.inf
        Om, Ob, h, rd, dm, s8 = d["Om"], d["Ob"], d["h"], d["rd"], d["dm"], d["s8"]
        bg = Background(Om, d.get("w0", -1.0), d.get("wa", 0.0))
        H0 = 100.0 * h
        r = mu_sn - bg.mu(z_sn, H0) - dm * mask
        chi2 = bao_chi2(bg, H0, rd, bao_all) + float(r @ invC_sn @ r) \
            + cmb_compressed_chi2(Om, Ob, h, bg) \
            + fs8_chi2(Om, s8, z_f, fs8v, invC_f, d.get("w0", -1.0), d.get("wa", 0.0)) \
            + ((H0 - H0_CEPHEID[0]) / H0_CEPHEID[1]) ** 2
        return -0.5 * chi2

    return logpost, names


def capstone_mode():
    out = []
    res = {}
    for model in ["lcdm", "w0wa"]:
        logpost, names = make_logpost(model)
        x0 = np.array([X0[n] for n in names])
        step = np.array([0.01 if n == "Om" else 0.005 if n == "Ob" else 0.01 if n == "h"
                         else 0.5 if n == "rd" else 0.02 if n in ("s8", "dm") else 0.05
                         for n in names])
        t0 = time.time()
        samples, acc = metropolis(logpost, x0, nsteps=18000, burn=5000, step=step, seed=91)
        stride = max(1, len(samples) // 3000)
        sub = samples[::stride]
        lps = np.array([logpost(s) for s in sub])
        chi2_map = float(-2.0 * lps.max())
        r = {"chi2_map": chi2_map, "acceptance": float(acc),
             "mean": dict(zip(names, map(float, samples.mean(axis=0)))),
             "std": dict(zip(names, map(float, samples.std(axis=0)))),
             "seconds": round(time.time() - t0, 1)}
        if model == "w0wa":
            wa = samples[:, 6]
            r["P_wa_lt_0"] = float((wa < 0).mean())
            r["P_wa_lt_-0.3"] = float((wa < -0.3).mean())
        np.save(os.path.join(RESULTS, f"chain_s9_capstone_{model}.npy"),
                np.column_stack([sub, lps]))
        res[model] = r
        print(f"  [{model}] chi2_map={chi2_map:.2f} acc={acc:.2f} "
              + " ".join(f"{n}={m:.3f}+-{s:.3f}" for n, m, s in
                         zip(names, samples.mean(axis=0), samples.std(axis=0)))
              + (f" P(wa<0)={r['P_wa_lt_0']:.3f}" if model == "w0wa" else "")
              + f" ({r['seconds']}s)", flush=True)
    d = res["w0wa"]["chi2_map"] - res["lcdm"]["chi2_map"]
    sig, p = sigma_from_dchi2(d, 2)
    entry = {"nuisances": ["dm(z<0.1) step", "Cepheid H0 prior 73.04+-1.04", "fs8 growth"],
             "combo": "BAO+Dovekie+CMB+fs8", "models": res,
             "w0wa_vs_lcdm": {"dchi2": float(d), "sigma": float(sig), "p": float(p)}}
    print(f"  >> capstone w0wa vs lcdm: dchi2={d:+.2f} -> {abs(sig):.2f} sigma", flush=True)
    with open(os.path.join(RESULTS, "capstone_results.json"), "w") as f:
        json.dump(entry, f, indent=2)
    print("CAPSTONE DONE", flush=True)


if __name__ == "__main__":
    capstone_mode()
