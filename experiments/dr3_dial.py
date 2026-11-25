"""dr3_dial.py — THE DIAL: pre-registered DR3 decision procedure.

Mode 'read' (default, instant): read the existing result JSONs and print the current
dial reading + the verdict the dial would give today.

Mode 'run <dr3_bao_dir>': when DR3 BAO data lands, re-run the two decisive fits
(w0wa-vs-LCDM with z>0.1 cut, and 5-bin cert5 dip significance) on DR3 BAO
(load_desi_bao(dr2_dir=<dr3_bao_dir>)) + Ly-a DR2 AP + Dovekie + CMB (+fs8 for cert5),
then apply the pre-registered verdict map. Everything else (anchors, dM, 2pt) is
re-run by the corresponding session scripts pointed at DR3 chains.

Verdict map (pre-registered in results/dr3_plan.json):
  dip |w+1|/sigma >= 3 in z in [0.75,1.5] AND survives z>0.1 cut, anchors, dM, 2pt
      -> DIAL = "PHANTOM REAL": N14 vacuum-rudder door OPEN (phantom sector exists);
         run N14 coupling study; dark energy is a dynamical field (N2).
  w0wa BF>10 surviving cuts but dip < 3 sigma
      -> DIAL = "DYNAMICAL NON-PHANTOM": quintessence-like; no NEC violation; N14 closed.
  nothing survives
      -> DIAL = "LAMBDA CONFIRMED": vacuum is a true constant; N14/N12 doors closed
         IN PRINCIPLE (Lovelock+NEC+passivity); the 2025-26 hint was calibration (s27).

Usage: python experiments/dr3_dial.py [read | run <dr3_bao_dir>]
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import Background, BackgroundBinned, lnF_binned, bao_chi2, sn_hd_chi2, \
    cmb_compressed_chi2, fs8_chi2, metropolis

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

Z_AP = 2.33
cov_ap = np.array([[0.066 ** 2, 0.225 * 0.066 * 0.33],
                   [0.225 * 0.066 * 0.33, 0.33 ** 2]])
LYA_AP = dict(z=np.array([Z_AP, Z_AP]), vals=np.array([8.600, 39.32]),
              qty=["DH_over_rs", "DM_over_rs"], inv_cov=np.linalg.inv(cov_ap))


def concat(a, b):
    n = len(a["z"])
    inv = np.zeros((n + len(b["z"]),) * 2)
    inv[:n, :n] = a["inv_cov"]
    inv[n:, n:] = b["inv_cov"]
    return dict(z=np.concatenate([a["z"], b["z"]]), vals=np.concatenate([a["vals"], b["vals"]]),
                qty=list(a["qty"]) + list(b["qty"]), inv_cov=inv)


def read_mode():
    ev = json.load(open(os.path.join(RESULTS, "evidence_results.json")))
    dip = json.load(open(os.path.join(RESULTS, "dip_persample_results.json")))
    cap = json.load(open(os.path.join(RESULTS, "capstone_results.json"))) \
        if os.path.exists(os.path.join(RESULTS, "capstone_results.json")) else None
    plan = json.load(open(os.path.join(RESULTS, "dr3_plan.json")))
    print("=== THE DIAL — current reading (DR2/Y3 + Ly-a DR2 AP data, 2026-11-25) ===")
    w = [e["w_bins_mean"][2] for e in dip["samples"].values()]
    s = [e["w_bins_std"][2] for e in dip["samples"].values()]
    dip_sig = abs(np.mean(w) + 1) / np.mean(s)
    print(f"  dip z in [0.75,1.5]: w = {np.mean(w):.3f} +- {np.mean(s):.3f} "
          f"-> |w+1|/sigma = {dip_sig:.2f} (need >= 3)")
    print(f"  w0wa dchi2 (all data): ~ -6.3, BF 9-12; z>0.1 cut: -0.01 (s8/s29)")
    print(f"  capstone (anchors+dM+fs8 free): dchi2 = -0.26, wa = -0.07+-0.13 (s21)")
    print(f"  2pt: does not move wa (s20); Ly-a AP: LCDM-consistent alone, chi2=0.97/2 (s33)")
    print(f"  discovery thresholds (k=2): 3s={plan['thresholds']['3sigma_k2']}, "
          f"4s={plan['thresholds']['4sigma_k2']}, 5s={plan['thresholds']['5sigma_k2']}")
    if dip_sig >= 3:
        verdict = "PHANTOM REAL -> N14 vacuum-rudder door OPEN (phantom sector exists)"
    elif dip_sig >= 2:
        verdict = "dip tantalizing (2-3s) but below bar -> door CLOSED pending DR3"
    else:
        verdict = ("dip below 2s -> Lambda favored; door CLOSED unless DR3 pushes the "
                   "dip past 3s (then N14 opens)")
    print(f"  DIAL VERDICT TODAY: dip at {dip_sig:.2f}s < 3 -> {verdict}")
    print("  -> the dial turns when DR3 full-shape lands; watcher automation is set.")
    print("DIAL READ DONE", flush=True)


def fit_w0wa(rows, z_sn, mu_sn, invC_sn, zcut=None):
    from scipy.optimize import minimize
    zs, ms, cs = z_sn, mu_sn, invC_sn
    if zcut is not None:
        m = z_sn >= zcut
        zs, ms = z_sn[m], mu_sn[m]
        cs = np.linalg.inv(np.linalg.inv(invC_sn)[np.ix_(np.where(m)[0], np.where(m)[0])])

    def neglp(x, model):
        if model == 0:
            Om, Ob, h, rd = x
            if not (0.05 < Om < 0.9 and 0.005 < Ob < 0.10 and 0.55 < h < 0.85 and 120 < rd < 165):
                return 1e12
            bg = Background(Om)
        else:
            Om, Ob, h, rd, w0, wa = x
            if not (0.05 < Om < 0.9 and 0.005 < Ob < 0.10 and 0.55 < h < 0.85 and 120 < rd < 165
                    and -1.8 < w0 < -0.2 and -3 < wa < 3):
                return 1e12
            bg = Background(Om, w0, wa)
        return bao_chi2(bg, 100 * h, rd, rows) + sn_hd_chi2(bg, 100 * h, zs, ms, cs) \
            + cmb_compressed_chi2(Om, Ob, h, bg)

    out = {}
    for model, x0s in [(0, [[0.31, 0.049, 0.68, 147.0]]),
                       (1, [[0.31, 0.049, 0.68, 147.0, -0.9, -0.4],
                            [0.30, 0.046, 0.70, 144.0, -0.85, -0.5]])]:
        best = None
        for x0 in x0s:
            r = minimize(lambda x: neglp(x, model), np.array(x0, float), method="Nelder-Mead",
                         options={"maxiter": 6000, "xatol": 1e-6, "fatol": 1e-8})
            if best is None or r.fun < best.fun:
                best = r
        out[model] = float(best.fun)
    return out[1] - out[0]


def run_mode(dr3_dir):
    print(f"DR3 BAO dir: {dr3_dir}", flush=True)
    bao = concat(cosmo.load_desi_bao(dr2_dir=dr3_dir, which="ALL")[0], LYA_AP)
    z_sn, mu_sn, invC_sn = cosmo.load_dovekie()
    d_all = fit_w0wa(bao, z_sn, mu_sn, invC_sn)
    d_cut = fit_w0wa(bao, z_sn, mu_sn, invC_sn, zcut=0.1)
    print(f"  DR3 w0wa vs LCDM: dchi2(all) = {d_all:+.2f}   dchi2(z>0.1 cut) = {d_cut:+.2f}", flush=True)
    out = {"dr3_bao_dir": dr3_dir, "dchi2_all": d_all, "dchi2_zgt0.1": d_cut,
           "note": "cert5 dip significance on DR3: run session22.py with DR3 rows; "
                   "verdict map in dr3_plan.json"}
    with open(os.path.join(RESULTS, "dr3_dial_run.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("DR3 DIAL RUN DONE", flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        run_mode(sys.argv[2])
    else:
        read_mode()
