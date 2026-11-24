"""session15.py — calibrator-host test: is the low-z SN offset concentrated in the
Cepheid-calibrator hosts?

Pantheon+ flags: IS_CALIBRATOR (host has Cepheid distance), USED_IN_SH0ES_HF.
Geometry-only LCDM (BAO+CMB), monopole marginalized per subsample; compare low-z
(z<0.1, z<0.3) residual means for calibrator vs non-calibrator hosts. If the +0.03 mag
offset is a Cepheid-anchor systematic, it should appear (or be larger) in the
calibrator subset; if cosmic (bubble), it must be the same in both.

Usage: python experiments/session15.py calib
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import Background, bao_chi2, cmb_compressed_chi2

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

print("loading data ...", flush=True)
bao_all = cosmo.load_desi_bao(which="ALL")[0]
z_p, mB_p, invC_p = cosmo.load_pantheonplus()
raw = np.genfromtxt(os.path.join(cosmo.DATA, "sn_data", "PantheonPlus", "Pantheon+SH0ES.dat"),
                    names=True, dtype=None, encoding="utf-8")
IS_CAL = np.asarray(raw["IS_CALIBRATOR"], float)
USED_HF = np.asarray(raw["USED_IN_SH0ES_HF"], float)
print(f"  Pantheon+ {len(z_p)}; calibrators {int(IS_CAL.sum())}, used-in-HF {int(USED_HF.sum())}",
      flush=True)

from scipy.optimize import minimize


def lp(x):
    Om, Ob, h, rd = x
    if not (0.05 < Om < 0.9 and 0.005 < Ob < 0.10 and 0.55 < h < 0.85 and 120 < rd < 165):
        return -np.inf
    bg = Background(Om)
    return -0.5 * (bao_chi2(bg, 100 * h, rd, bao_all) + cmb_compressed_chi2(Om, Ob, h, bg))


r0 = minimize(lambda x: -lp(x), [0.31, 0.049, 0.68, 147.0], method="Nelder-Mead",
              options={"maxiter": 3000})
bg = Background(r0.x[0])
H0 = 100 * r0.x[2]
print(f"  geometry LCDM H0={H0:.2f}", flush=True)

out = {"geometry_H0": float(H0)}
groups = {"calibrator": IS_CAL > 0.5, "non_calibrator": IS_CAL <= 0.5,
          "used_in_SH0ES_HF": USED_HF > 0.5, "not_used_HF": USED_HF <= 0.5}
for zcut in [0.1, 0.3]:
    out[f"z<{zcut}"] = {}
    for gname, gmask in groups.items():
        m = (z_p < zcut) & gmask
        idx = np.where(m)[0]
        if len(idx) < 5:
            continue
        C = invC_p[np.ix_(idx, idx)]
        r = mB_p[idx] - bg.mu(z_p[idx], H0)
        one = np.ones(len(idx))
        w = C @ one
        mean = float(r @ w / (one @ w))
        sig = float(np.sqrt(1.0 / (one @ w)))
        out[f"z<{zcut}"][gname] = {"n": int(len(idx)), "mean_resid": mean, "sigma": sig}
        print(f"  z<{zcut} {gname:18s} n={len(idx):4d} resid={mean:+.4f}+-{sig:.4f} mag",
              flush=True)
    a = out[f"z<{zcut}"].get("calibrator")
    b = out[f"z<{zcut}"].get("non_calibrator")
    if a and b:
        d = a["mean_resid"] - b["mean_resid"]
        s = np.hypot(a["sigma"], b["sigma"])
        out[f"z<{zcut}"]["cal_minus_noncal"] = {"delta": float(d), "sigma": float(s),
                                                "tension": float(abs(d) / s)}
        print(f"  >> z<{zcut} calibrator-minus-noncalibrator = {d:+.4f}+-{s:.4f} "
              f"({abs(d)/s:.2f}s)", flush=True)

with open(os.path.join(RESULTS, "calib_results.json"), "w") as f:
    json.dump(out, f, indent=2)
print("CALIB DONE", flush=True)
