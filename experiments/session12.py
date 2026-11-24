"""session12.py — edge-sharpness test: is the low-z SN step a sharp wall or a smooth drift?

Model: mu_model(z) = mu_LCDM/geometry(z) + dM * sigmoid((z_step - z)/sz)
with z_step fixed at 0.1, dM and logistic width sz free (sz -> 0 = sharp step).
Fit on Dovekie (full) + BAO + CMB-compressed, profile chi2 over sz grid; compare
sharp-step vs smooth (sz large) vs no step (dM=0). If data prefer sz ~ 0, that's the
signature of a phase-transition bubble wall (N4/N5) rather than smooth calibration drift.

Usage: python experiments/session12.py edge
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

from scipy.optimize import minimize

Z_STEP = 0.1


def chi2_at(sz, dM, x):
    Om, Ob, h, rd = x
    if not (0.05 < Om < 0.9 and 0.005 < Ob < 0.10 and 0.55 < h < 0.85 and 120 < rd < 165):
        return 1e12
    bg = Background(Om)
    H0 = 100.0 * h
    r = mu_sn - bg.mu(z_sn, H0) - dM / (1.0 + np.exp((z_sn - Z_STEP) / max(sz, 1e-6)))
    return float(bao_chi2(bg, H0, rd, bao_all) + r @ invC_sn @ r
                 + cmb_compressed_chi2(Om, Ob, h, bg))


def edge_mode():
    rows = []
    for sz in [1e-4, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 2.0]:
        best = None
        for x0 in ([0.31, 0.049, 0.68, 147.0, 0.03], [0.31, 0.049, 0.70, 147.0, 0.0]):
            res = minimize(lambda v: chi2_at(sz, v[4], v[:4]), np.array(x0, float),
                           method="Nelder-Mead", options={"maxiter": 3000, "xatol": 1e-5})
            if res.fun < (best[0] if best else np.inf):
                best = (res.fun, res.x)
        rows.append({"sz": sz, "chi2": float(best[0]), "dM": float(best[1][4]),
                     "Om": float(best[1][0]), "h": float(best[1][2])})
        print(f"  sz={sz:6.3f} chi2={best[0]:.2f} dM={best[1][4]:+.4f} "
              f"Om={best[1][0]:.3f} h={best[1][2]:.3f}", flush=True)
    # no-step reference
    res0 = minimize(lambda v: chi2_at(1e-4, 0.0, v), [0.31, 0.049, 0.68, 147.0],
                    method="Nelder-Mead", options={"maxiter": 3000})
    chi2_nostep = float(res0.fun)
    best = min(rows, key=lambda r: r["chi2"])
    entry = {"z_step": Z_STEP, "chi2_nostep": chi2_nostep, "rows": rows,
             "best": best, "dchi2_best_vs_nostep": best["chi2"] - chi2_nostep}
    print(f"  >> no-step chi2={chi2_nostep:.2f}; best sz={best['sz']} dM={best['dM']:+.4f} "
          f"dchi2={entry['dchi2_best_vs_nostep']:+.2f}", flush=True)
    with open(os.path.join(RESULTS, "edge_results.json"), "w") as f:
        json.dump(entry, f, indent=2)
    print("EDGE DONE", flush=True)


if __name__ == "__main__":
    edge_mode()
