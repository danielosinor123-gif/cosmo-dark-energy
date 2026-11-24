"""session18.py — DR3 forecast: how much BAO precision breaks the w0wa <-> dM-step
degeneracy?

Scale the ALL BAO inverse covariance by f^2 (f = effective BAO precision gain of
DR3-like data), keep Dovekie SN + CMB-compressed fixed, and profile-fit
  lcdm+dM  vs  w0wa+dM  vs  w0wa (no dM)
Report delta-chi2 as a function of f. When delta-chi2(w0wa+dM vs lcdm+dM) < -11.8
(3 sigma, k=2), evolving DE becomes distinguishable from a calibration step.

Usage: python experiments/session18.py forecast
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import Background, bao_chi2, sn_hd_chi2, cmb_compressed_chi2

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

print("loading data ...", flush=True)
bao_all = cosmo.load_desi_bao(which="ALL")[0]
z_sn, mu_sn, invC_sn = cosmo.load_dovekie()
mask = (z_sn < 0.1).astype(float)

from scipy.optimize import minimize


def chi2(x, model, f):
    if model == "lcdm":
        Om, Ob, h, rd, dm = x
        w0 = wa = 0.0
        if not (0.05 < Om < 0.9 and 0.005 < Ob < 0.10 and 0.55 < h < 0.85 and 120 < rd < 165
                and -0.3 < dm < 0.3):
            return 1e12
    elif model == "w0wa":
        Om, Ob, h, rd, w0, wa, dm = x
        if not (0.05 < Om < 0.9 and 0.005 < Ob < 0.10 and 0.55 < h < 0.85 and 120 < rd < 165
                and -1.8 < w0 < -0.2 and -3 < wa < 3 and -0.3 < dm < 0.3):
            return 1e12
    else:  # w0wa_nodm
        Om, Ob, h, rd, w0, wa = x
        dm = 0.0
        if not (0.05 < Om < 0.9 and 0.005 < Ob < 0.10 and 0.55 < h < 0.85 and 120 < rd < 165
                and -1.8 < w0 < -0.2 and -3 < wa < 3):
            return 1e12
    bg = Background(Om, w0 if model != "lcdm" else -1.0, wa if model != "lcdm" else 0.0)
    H0 = 100.0 * h
    r = mu_sn - bg.mu(z_sn, H0) - dm * mask
    # scale BAO precision by f: chi2_bao * f^2
    return float(f * f * bao_chi2(bg, H0, rd, bao_all) + r @ invC_sn @ r
                 + cmb_compressed_chi2(Om, Ob, h, bg))


def fit(model, f):
    if model == "lcdm":
        x0s = [[0.31, 0.049, 0.68, 147.0, 0.0], [0.30, 0.045, 0.70, 144.0, 0.03]]
    elif model == "w0wa":
        x0s = [[0.31, 0.049, 0.68, 147.0, -0.9, -0.4, 0.02], [0.30, 0.046, 0.70, 144.0, -0.85, -0.5, 0.0]]
    else:
        x0s = [[0.31, 0.049, 0.68, 147.0, -0.9, -0.4], [0.30, 0.046, 0.70, 144.0, -0.85, -0.5]]
    best = None
    for x0 in x0s:
        res = minimize(lambda x: chi2(x, model, f), np.array(x0, float), method="Nelder-Mead",
                       options={"maxiter": 4000, "xatol": 1e-5, "fatol": 1e-7})
        if best is None or res.fun < best.fun:
            best = res
    return float(best.fun), best.x


def forecast_mode():
    rows = []
    for f in [1.0, 1.5, 2.0, 3.0, 4.0, 6.0]:
        c_lcdm, _ = fit("lcdm", f)
        c_w0wa, xw = fit("w0wa", f)
        c_w0wa_nodm, _ = fit("w0wa_nodm", f)
        d1 = c_w0wa - c_lcdm          # w0wa+dM vs lcdm+dM
        d2 = c_w0wa_nodm - c_lcdm     # w0wa (no step) vs lcdm+dM
        rows.append({"f_bao": f, "chi2_lcdm_dM": c_lcdm, "chi2_w0wa_dM": c_w0wa,
                     "chi2_w0wa_nodM": c_w0wa_nodm,
                     "dchi2_w0wa_dM_vs_lcdm_dM": d1, "dchi2_w0wa_nodM_vs_lcdm_dM": d2,
                     "wa_map": float(xw[5]) if len(xw) > 5 else None,
                     "dm_map": float(xw[6]) if len(xw) > 6 else 0.0})
        print(f"  f={f:4.1f}: dchi2(w0wa+dM vs lcdm+dM)={d1:+.2f}  "
              f"dchi2(w0wa no-step vs lcdm+dM)={d2:+.2f}  wa_map={rows[-1]['wa_map']:+.3f} "
              f"dM_map={rows[-1]['dm_map']:+.4f}", flush=True)
    entry = {"setup": "BAO ALL cov scaled by f^-2; Dovekie+CMB fixed; dM step z<0.1 free in lcdm/w0wa",
             "discovery_threshold_dchi2_k2_3sigma": -11.83, "rows": rows}
    with open(os.path.join(RESULTS, "forecast_dr3.json"), "w") as f:
        json.dump(entry, f, indent=2)
    print("FORECAST DONE", flush=True)


if __name__ == "__main__":
    forecast_mode()
