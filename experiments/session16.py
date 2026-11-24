"""session16.py — peculiar-velocity sensitivity of the calibrator offset (caveat of §27).

Recompute the calibrator-vs-noncalibrator low-z offset using different redshift columns:
zHD (default, v_pec-corrected), zCMB (no v_pec correction), zHEL (heliocentric).
If the +0.20 mag calibrator offset is a v_pec artifact, it should change a lot when
switching columns; if it is the Cepheid anchor, it should be stable.

Usage: python experiments/session16.py vpec
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

bao_all = cosmo.load_desi_bao(which="ALL")[0]
z_p, mB_p, invC_p = cosmo.load_pantheonplus()
raw = np.genfromtxt(os.path.join(cosmo.DATA, "sn_data", "PantheonPlus", "Pantheon+SH0ES.dat"),
                    names=True, dtype=None, encoding="utf-8")
IS_CAL = np.asarray(raw["IS_CALIBRATOR"], float)
Z = {"zHD": np.asarray(raw["zHD"], float), "zCMB": np.asarray(raw["zCMB"], float),
     "zHEL": np.asarray(raw["zHEL"], float)}

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
print(f"geometry LCDM H0={H0:.2f}", flush=True)

out = {"geometry_H0": float(H0)}
for zname, zz in Z.items():
    for zcut in [0.1, 0.3]:
        m_all = zz < zcut
        idx_all = np.where(m_all)[0]
        C = invC_p[np.ix_(idx_all, idx_all)]
        r = mB_p[idx_all] - bg.mu(zz[idx_all], H0)
        one = np.ones(len(idx_all))
        w = C @ one
        r = r - (r @ w) / (one @ w)
        res = {}
        for gname, gmask in [("calibrator", IS_CAL > 0.5), ("non_calibrator", IS_CAL <= 0.5)]:
            ii = np.where(gmask[idx_all])[0]
            Cii = C[np.ix_(ii, ii)]
            obi = np.ones(len(ii))
            wi = Cii @ obi
            res[gname] = {"n": int(len(ii)),
                          "mean": float(r[ii] @ wi / (obi @ wi)),
                          "sigma": float(np.sqrt(1.0 / (obi @ wi)))}
        d = res["calibrator"]["mean"] - res["non_calibrator"]["mean"]
        s = float(np.hypot(res["calibrator"]["sigma"], res["non_calibrator"]["sigma"]))
        res["cal_minus_noncal"] = {"delta": float(d), "sigma": s, "tension": float(abs(d) / s)}
        out[f"{zname}_z<{zcut}"] = res
        print(f"  {zname} z<{zcut}: cal-noncal = {d:+.4f}+-{s:.4f} ({abs(d)/s:.2f}s) "
              f"(cal n={res['calibrator']['n']})", flush=True)

with open(os.path.join(RESULTS, "vpec_results.json"), "w") as f:
    json.dump(out, f, indent=2)
print("VPEC DONE", flush=True)
