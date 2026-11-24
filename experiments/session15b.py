"""session15b.py — calibrator-host test, done carefully: ONE global monopole M
marginalized over the whole z<0.3 sample, then compare calibrator vs non-calibrator
group means (and vs used-in-SH0ES-HF). Cross-check with mB column too.
Also print z distributions to check for selection effects.
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
USED_HF = np.asarray(raw["USED_IN_SH0ES_HF"], float)
mB_raw = np.asarray(raw["mB"], float)
zHD = np.asarray(raw["zHD"], float)

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
print(f"  calibrator z range: {zHD[IS_CAL>0.5].min():.4f}..{zHD[IS_CAL>0.5].max():.4f} "
      f"median {np.median(zHD[IS_CAL>0.5]):.4f}", flush=True)

out = {"geometry_H0": float(H0)}
for zcut in [0.1, 0.3]:
    m_all = z_p < zcut
    idx_all = np.where(m_all)[0]
    C = invC_p[np.ix_(idx_all, idx_all)]
    for col_name, col in [("m_b_corr", mB_p), ("mB_raw", mB_raw)]:
        r = col[idx_all] - bg.mu(z_p[idx_all], H0)
        one = np.ones(len(idx_all))
        w = C @ one
        r = r - (r @ w) / (one @ w)  # ONE global monopole
        res = {}
        for gname, gmask in [("calibrator", IS_CAL > 0.5), ("non_calibrator", IS_CAL <= 0.5),
                             ("used_HF", USED_HF > 0.5), ("not_used_HF", USED_HF <= 0.5)]:
            ii = np.where(gmask[idx_all])[0]
            if len(ii) < 5:
                continue
            Cii = C[np.ix_(ii, ii)]
            obi = np.ones(len(ii))
            wi = Cii @ obi
            res[gname] = {"n": int(len(ii)),
                          "mean_resid": float(r[ii] @ wi / (obi @ wi)),
                          "sigma": float(np.sqrt(1.0 / (obi @ wi)))}
        a, b = res.get("calibrator"), res.get("non_calibrator")
        if a and b:
            d = a["mean_resid"] - b["mean_resid"]
            s = float(np.hypot(a["sigma"], b["sigma"]))
            res["cal_minus_noncal"] = {"delta": float(d), "sigma": s, "tension": float(abs(d) / s)}
        out[f"z<{zcut}_{col_name}"] = res
        print(f"  z<{zcut} {col_name}: " + " | ".join(
            f"{k}: {v['mean_resid']:+.4f}+-{v['sigma']:.4f} (n={v['n']})"
            for k, v in res.items() if k != "cal_minus_noncal")
            + (f" || cal-noncal={res['cal_minus_noncal']['delta']:+.4f}"
               f"+-{res['cal_minus_noncal']['sigma']:.4f} "
               f"({res['cal_minus_noncal']['tension']:.2f}s)" if "cal_minus_noncal" in res else ""),
            flush=True)

with open(os.path.join(RESULTS, "calib_results.json"), "w") as f:
    json.dump(out, f, indent=2)
print("CALIB2 DONE", flush=True)
