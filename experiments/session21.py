"""session21.py — fold DESI DR2 Lyman-alpha full-shape AP (arXiv:2607.27410) into the
pipeline. AP measurement at z_eff=2.33 (1% precision, 2x tighter than Ly-a BAO):
    DH/rd = 8.600 +- 0.066
    DM/rd = 39.32 +- 0.33
    rho(DM, DH) = 0.225
Paper headline: Ly-a+BBN gives H0 = 66.5 +- 1.3 (LambdaCDM-consistent), Omega_m =
0.325 +- 0.018 (1.4 sigma above DESI BAO) — the high-z anchor moved TOWARD LambdaCDM.

Here: w0wa vs LCDM with Dovekie SN + CMB fixed, BAO varied:
  all      : ALL 13 rows (baseline, dchi2=-6.49 from S14)
  all+ap   : ALL 13 + Lya DR2 AP 2
  noLya+ap : ALL minus Lya rows (11) + Lya DR2 AP 2  (AP replaces Y3 Ly-a BAO)
  ap_only  : Lya DR2 AP alone + SN + CMB

Usage: python experiments/session21.py lyaap
"""
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import Background, bao_chi2, sn_hd_chi2, cmb_compressed_chi2

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

# --- DESI DR2 Lya full-shape AP (arXiv:2607.27410, Table in Sec. V) -------------
Z_AP = 2.33
DH_AP, S_DH = 8.600, 0.066
DM_AP, S_DM = 39.32, 0.33
RHO = 0.225
cov_ap = np.array([[S_DH * S_DH, RHO * S_DH * S_DM],
                   [RHO * S_DH * S_DM, S_DM * S_DM]])
LYA_AP = dict(name="Lya_DR2_AP", z=np.array([Z_AP, Z_AP]),
              vals=np.array([DH_AP, DM_AP]), qty=["DH_over_rs", "DM_over_rs"],
              inv_cov=np.linalg.inv(cov_ap))

# persist the public numbers for reproducibility
d = os.path.join(cosmo.DATA, "bao_data", "desi_lya_dr2_ap")
os.makedirs(d, exist_ok=True)
with open(os.path.join(d, "lya_dr2_ap_mean.txt"), "w") as f:
    f.write("# DESI DR2 Lya full-shape AP, arXiv:2607.27410 (z_eff=2.33)\n")
    f.write("# z value quantity\n")
    f.write(f"{Z_AP} {DH_AP} DH_over_rs\n{Z_AP} {DM_AP} DM_over_rs\n")
np.savetxt(os.path.join(d, "lya_dr2_ap_cov.txt"), cov_ap)


def concat_rows(a, b):
    n = len(a["z"])
    inv = np.zeros((n + len(b["z"]),) * 2)
    inv[:n, :n] = a["inv_cov"]
    inv[n:, n:] = b["inv_cov"]
    return dict(name=a["name"] + "+" + b["name"], z=np.concatenate([a["z"], b["z"]]),
                vals=np.concatenate([a["vals"], b["vals"]]),
                qty=list(a["qty"]) + list(b["qty"]), inv_cov=inv)


def drop_lya(rows):
    keep = [i for i, q in enumerate(rows["qty"]) if True]
    # Lya rows in ALL are the z=2.33 entries
    idx = [i for i, z in enumerate(rows["z"]) if abs(z - 2.33) > 1e-6]
    n = len(idx)
    inv = rows["inv_cov"][np.ix_(idx, idx)]
    return dict(name="ALL-minus-Lya", z=rows["z"][idx], vals=rows["vals"][idx],
                qty=[rows["qty"][i] for i in idx], inv_cov=inv)


def fit(bg_rows, model, z_sn, mu_sn, invC_sn):
    if model == "lcdm":
        def neglp(x):
            Om, Ob, h, rd = x
            if not (0.05 < Om < 0.9 and 0.005 < Ob < 0.10 and 0.55 < h < 0.85 and 120 < rd < 165):
                return 1e12
            bg = Background(Om)
            return bao_chi2(bg, 100 * h, rd, bg_rows) + sn_hd_chi2(bg, 100 * h, z_sn, mu_sn, invC_sn) \
                + cmb_compressed_chi2(Om, Ob, h, bg)
        x0s = [[0.31, 0.049, 0.68, 147.0], [0.30, 0.046, 0.70, 144.0]]
    else:
        def neglp(x):
            Om, Ob, h, rd, w0, wa = x
            if not (0.05 < Om < 0.9 and 0.005 < Ob < 0.10 and 0.55 < h < 0.85 and 120 < rd < 165
                    and -1.8 < w0 < -0.2 and -3 < wa < 3):
                return 1e12
            bg = Background(Om, w0, wa)
            return bao_chi2(bg, 100 * h, rd, bg_rows) + sn_hd_chi2(bg, 100 * h, z_sn, mu_sn, invC_sn) \
                + cmb_compressed_chi2(Om, Ob, h, bg)
        x0s = [[0.31, 0.049, 0.68, 147.0, -0.9, -0.4], [0.30, 0.046, 0.70, 144.0, -0.85, -0.5],
               [0.31, 0.049, 0.70, 147.0, -1.0, 0.0]]
    best = None
    for x0 in x0s:
        r = minimize(neglp, np.array(x0, float), method="Nelder-Mead",
                     options={"maxiter": 6000, "xatol": 1e-6, "fatol": 1e-8})
        if best is None or r.fun < best.fun:
            best = r
    return float(best.fun), best.x


def lyaap_mode():
    print("loading data ...", flush=True)
    bao_all = cosmo.load_desi_bao(which="ALL")[0]
    z_sn, mu_sn, invC_sn = cosmo.load_dovekie()
    nomlya = drop_lya(bao_all)
    print(f"  ALL rows={len(bao_all['z'])}, minus-Lya={len(nomlya['z'])}", flush=True)

    combos = {
        "all": bao_all,
        "all+lyaAP": concat_rows(bao_all, LYA_AP),
        "noLya+lyaAP": concat_rows(nomlya, LYA_AP),
        "lyaAP_only": LYA_AP,
    }
    out = {"lya_ap": {"z": Z_AP, "DH_over_rd": [DH_AP, S_DH], "DM_over_rd": [DM_AP, S_DM],
                      "rho": RHO, "source": "arXiv:2607.27410"}, "combos": {}}
    for cname, rows in combos.items():
        c_l, xl = fit(rows, "lcdm", z_sn, mu_sn, invC_sn)
        c_w, xw = fit(rows, "w0wa", z_sn, mu_sn, invC_sn)
        d = c_w - c_l
        entry = {"chi2_lcdm": c_l, "chi2_w0wa": c_w, "dchi2": d,
                 "h_lcdm": float(100 * xl[2]), "h_w0wa": float(100 * xw[2]),
                 "w0": float(xw[4]), "wa": float(xw[5])}
        out["combos"][cname] = entry
        print(f"  {cname:14s}: dchi2={d:+.2f}  lcdm chi2={c_l:.2f} w0wa chi2={c_w:.2f} "
              f"wa={xw[5]:+.3f} w0={xw[4]:+.3f} H0_w0wa={100*xw[2]:.2f}", flush=True)

    # also: pure Lya-AP geometry test vs LCDM fiducial (no SN/CMB) — H0 content
    bg0 = Background(0.315)
    r = cosmo.bao_model_vec(bg0, 67.7, 147.1, LYA_AP) - LYA_AP["vals"]
    out["lyaAP_only_vs_lcdm_H067.7_chi2"] = float(r @ LYA_AP["inv_cov"] @ r)
    print(f"  LyaAP-only chi2 vs LCDM(Om=0.315,H0=67.7,rd=147.1): {out['lyaAP_only_vs_lcdm_H067.7_chi2']:.2f}",
          flush=True)

    with open(os.path.join(RESULTS, "lyaap_results.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("LYAAP DONE", flush=True)


if __name__ == "__main__":
    lyaap_mode()
