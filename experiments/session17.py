"""session17.py — cert5 without Cepheid calibrators: does the z~1 dip survive?

5-bin w(z) + growth on BAO + Pantheon+ (IS_CALIBRATOR==0 only) + CMB + fs8.
Compare w bins vs the full-sample cert5 (§24). If the z 0.75-1.5 dip survives with
non-calibrator SNe only, it is not an anchor artifact.

Usage: python experiments/session17.py cert5noc
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import BackgroundBinned, lnF_binned, bao_chi2, sn_hd_chi2, \
    cmb_compressed_chi2, fs8_chi2, metropolis

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

print("loading data ...", flush=True)
bao_all = cosmo.load_desi_bao(which="ALL")[0]
z_p, mB_p, invC_p = cosmo.load_pantheonplus()
raw = np.genfromtxt(os.path.join(cosmo.DATA, "sn_data", "PantheonPlus", "Pantheon+SH0ES.dat"),
                    names=True, dtype=None, encoding="utf-8")
IS_CAL = np.asarray(raw["IS_CALIBRATOR"], float)
keep = IS_CAL <= 0.5
z_sn = z_p[keep]
mB_sn = mB_p[keep]
# correct marginalization: submatrix of the COVARIANCE, then invert
C_full = np.linalg.inv(invC_p)
idx = np.where(keep)[0]
invC_sn = np.linalg.inv(C_full[np.ix_(idx, idx)])
z_f, fs8v, invC_f = cosmo.load_fs8()
print(f"  Pantheon+ non-calibrator SNe: {keep.sum()} of {len(z_p)}", flush=True)

EDGES5 = [0.0, 0.3, 0.75, 1.5, 2.5, np.inf]
PRIORS = {"Om": (0.05, 0.9), "Ob": (0.005, 0.10), "h": (0.55, 0.85), "rd": (120.0, 165.0),
          "s8": (0.2, 1.4)}
for i in range(5):
    PRIORS[f"w{i}"] = (-3.0, 1.0)
names = ["Om", "Ob", "h", "rd", "s8"] + [f"w{i}" for i in range(5)]


def logpost(x):
    d = dict(zip(names, x))
    for n in names:
        lo, hi = PRIORS[n]
        if not (lo <= d[n] <= hi):
            return -np.inf
    Om, Ob, h, rd, s8 = d["Om"], d["Ob"], d["h"], d["rd"], d["s8"]
    wb = [d[f"w{i}"] for i in range(5)]
    bg = BackgroundBinned(Om, wb, EDGES5)
    lnF = lambda z: lnF_binned(z, EDGES5, wb)
    H0 = 100.0 * h
    chi2 = bao_chi2(bg, H0, rd, bao_all) + sn_hd_chi2(bg, H0, z_sn, mB_sn, invC_sn) \
        + cmb_compressed_chi2(Om, Ob, h, bg) \
        + fs8_chi2(Om, s8, z_f, fs8v, invC_f, lnF=lnF)
    return -0.5 * chi2


x0 = np.array([0.31, 0.049, 0.68, 147.0, 0.81] + [-0.93, -0.97, -1.2, -1.0, -1.1])
step = np.array([0.01, 0.005, 0.01, 0.5, 0.02] + [0.06] * 5)
t0 = time.time()
samples, acc = metropolis(logpost, x0, nsteps=24000, burn=6000, step=step, seed=141)
stride = max(1, len(samples) // 3000)
sub = samples[::stride]
lps = np.array([logpost(s) for s in sub])
wi = samples[:, 5:10]
entry = {"sample": "Pantheon+ non-calibrator", "edges": EDGES5,
         "chi2_map": float(-2.0 * lps.max()), "acceptance": float(acc),
         "w_bins_mean": wi.mean(axis=0).tolist(), "w_bins_std": wi.std(axis=0).tolist(),
         "seconds": round(time.time() - t0, 1)}
np.save(os.path.join(RESULTS, "chain_s17_cert5noc.npy"), np.column_stack([sub, lps]))
lbl = ["0-0.3", "0.3-0.75", "0.75-1.5", "1.5-2.5", ">2.5"]
print("  cert5-nocal w bins: " + " | ".join(
    f"{l}: {m:+.3f}+-{s:.3f}" for l, m, s in zip(lbl, wi.mean(axis=0), wi.std(axis=0)))
    + f"  chi2_map={entry['chi2_map']:.2f} acc={acc:.2f} ({entry['seconds']}s)", flush=True)
with open(os.path.join(RESULTS, "cert5noc_results.json"), "w") as f:
    json.dump(entry, f, indent=2)
print("CERT5NOC DONE", flush=True)
