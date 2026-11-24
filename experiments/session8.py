"""session8.py — N4 test: is vacuum energy locally coupled to matter (environment)?

Hypothesis N4 ("screenable / environment-dependent vacuum energy"): if rho_vac is an
order parameter that matter can locally tilt (screen/unscreen), low-z SN residuals
should correlate with host-galaxy stellar mass. Test on Pantheon+ (has HOST_LOGMASS):

  geometry-only LCDM from BAO+CMB -> mu_model(z); per-sample offset M marginalized;
  residuals r = m_b_corr - mu_model - M; then in z<0.3 (and z>0.3 control):
    - mean residual in host logM<10 vs >10 bins, difference + tension
    - weighted slope dr/dlogM and its significance

Usage: python experiments/session8.py env
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
d = os.path.join(cosmo.DATA, "sn_data", "PantheonPlus")
raw = np.genfromtxt(os.path.join(d, "Pantheon+SH0ES.dat"), names=True, dtype=None,
                    encoding="utf-8")
logM = np.asarray(raw["HOST_LOGMASS"], float)
print(f"  Pantheon+ {len(z_p)} rows, host logM range {np.nanmin(logM):.1f}..{np.nanmax(logM):.1f}",
      flush=True)

from scipy.optimize import minimize


def lp_geo(x):
    Om, Ob, h, rd = x
    if not (0.05 < Om < 0.9 and 0.005 < Ob < 0.10 and 0.55 < h < 0.85 and 120 < rd < 165):
        return -np.inf
    bg = Background(Om)
    return -0.5 * (bao_chi2(bg, 100 * h, rd, bao_all) + cmb_compressed_chi2(Om, Ob, h, bg))


r0 = minimize(lambda x: -lp_geo(x), [0.31, 0.049, 0.68, 147.0], method="Nelder-Mead",
              options={"maxiter": 3000})
Om, Ob, h, rd = r0.x
bg = Background(Om)
H0 = 100 * h
print(f"  geometry LCDM: Om={Om:.4f} H0={H0:.2f} chi2={-2*lp_geo(r0.x):.2f}", flush=True)

r = mB_p - bg.mu(z_p, H0)
one = np.ones(len(z_p))
C1 = invC_p @ one
r_off = r - (r @ C1) / (one @ C1) * one   # marginalize global M

out = {"geometry": {"Om": float(Om), "H0": float(H0), "chi2": float(-2 * lp_geo(r0.x))},
       "bins": {}}
for zlo, zhi, tag in [(0.0, 0.3, "z<0.3"), (0.3, 99.0, "z>0.3 (control)")]:
    m = (z_p >= zlo) & (z_p < zhi)
    idx = np.where(m)[0]
    Cbb = invC_p[np.ix_(idx, idx)]
    rbb = r_off[idx]
    lm = logM[idx]
    ob = np.ones(len(idx))
    w = Cbb @ ob
    mean = float(rbb @ w / (ob @ w))
    sig = float(np.sqrt(1.0 / (ob @ w)))
    # host-mass split at 10.0
    lo = lm < 10.0
    hi = ~lo
    res = {}
    for name, mm in [("logM<10", lo), ("logM>=10", hi)]:
        ii = np.where(mm)[0]
        Cii = Cbb[np.ix_(ii, ii)]
        obi = np.ones(len(ii))
        wi = Cii @ obi
        res[name] = {"mean": float(rbb[ii] @ wi / (obi @ wi)),
                     "sigma": float(np.sqrt(1.0 / (obi @ wi))), "n": int(len(ii))}
    dd = res["logM>=10"]["mean"] - res["logM<10"]["mean"]
    ss = np.hypot(res["logM>=10"]["sigma"], res["logM<10"]["sigma"])
    # weighted slope dr/dlogM
    X = np.column_stack([ob, lm - 10.0])
    XtWX = X.T @ Cbb @ X
    beta = np.linalg.solve(XtWX, X.T @ Cbb @ rbb)
    covb = np.linalg.inv(XtWX)
    slope, slope_sig = float(beta[1]), float(np.sqrt(covb[1, 1]))
    out["bins"][tag] = {"mean_resid": mean, "sigma": sig, "n": int(m.sum()),
                        "host_split": res, "delta_high_minus_low": float(dd),
                        "delta_sigma": float(ss), "delta_tension": float(abs(dd) / ss),
                        "slope_dr_dlogM": slope, "slope_sigma": slope_sig,
                        "slope_tension": float(abs(slope) / slope_sig)}
    print(f"  [{tag}] n={m.sum()} mean={mean:+.4f}+-{sig:.4f} | host split: "
          f"lowM {res['logM<10']['mean']:+.4f}+-{res['logM<10']['sigma']:.4f} (n={res['logM<10']['n']}), "
          f"highM {res['logM>=10']['mean']:+.4f}+-{res['logM>=10']['sigma']:.4f} (n={res['logM>=10']['n']}) "
          f"-> d={dd:+.4f}+-{ss:.4f} ({abs(dd)/ss:.2f}s) | slope dr/dlogM={slope:+.4f}+-{slope_sig:.4f} "
          f"({abs(slope)/slope_sig:.2f}s)", flush=True)

with open(os.path.join(RESULTS, "env_results.json"), "w") as f:
    json.dump(out, f, indent=2)
print("ENV DONE", flush=True)
