"""session23.py — two falsifiable pushes on the surviving z~1 dip + a "crazier" one.

Mode 'edges': does the dip stay at z~[0.75,1.5] under different bin edges?
  cert5-style 5-bin w(z)+growth on Dovekie + BAO(ALL+Lya DR2 AP) + CMB + fs8 with
  edge sets A=[0,0.3,0.75,1.5,2.5,inf], B=[0,0.5,1.0,2.0,3.0,inf], C=[0,0.3,1.0,1.75,2.5,inf].
  If the dip is a bin-edge artifact it moves/disappears; if physical it persists.

Mode 'dipole': Hooper-style (arXiv:1806.05203) harvesting-civilization imprint search —
  fit a dipole (and quadrupole) in the Pantheon+ SN Hubble residuals vs BAO+CMB LCDM
  geometry, full sample and z<0.1, per hemisphere. A civilization harvesting vacuum/stellar
  energy would imprint a sky asymmetry in the local expansion. Expect null.

Usage: python experiments/session23.py edges|dipole
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import Background, BackgroundBinned, lnF_binned, bao_chi2, sn_hd_chi2, sn_mB_chi2, \
    cmb_compressed_chi2, fs8_chi2, metropolis

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

Z_AP = 2.33
cov_ap = np.array([[0.066 ** 2, 0.225 * 0.066 * 0.33],
                   [0.225 * 0.066 * 0.33, 0.33 ** 2]])
LYA_AP = dict(name="Lya_DR2_AP", z=np.array([Z_AP, Z_AP]),
              vals=np.array([8.600, 39.32]), qty=["DH_over_rs", "DM_over_rs"],
              inv_cov=np.linalg.inv(cov_ap))


def concat_rows(a, b):
    n = len(a["z"])
    inv = np.zeros((n + len(b["z"]),) * 2)
    inv[:n, :n] = a["inv_cov"]
    inv[n:, n:] = b["inv_cov"]
    return dict(z=np.concatenate([a["z"], b["z"]]), vals=np.concatenate([a["vals"], b["vals"]]),
                qty=list(a["qty"]) + list(b["qty"]), inv_cov=inv)


def edges_mode():
    print("loading data ...", flush=True)
    bao = concat_rows(cosmo.load_desi_bao(which="ALL")[0], LYA_AP)
    z_d, mu_d, invC_d = cosmo.load_dovekie()
    z_f, fs8v, invC_f = cosmo.load_fs8()

    EDGE_SETS = {
        "A_0.3_0.75_1.5_2.5": [0.0, 0.3, 0.75, 1.5, 2.5, np.inf],
        "B_0.5_1.0_2.0_3.0": [0.0, 0.5, 1.0, 2.0, 3.0, np.inf],
        "C_0.3_1.0_1.75_2.5": [0.0, 0.3, 1.0, 1.75, 2.5, np.inf],
    }
    PRIORS = {"Om": (0.05, 0.9), "Ob": (0.005, 0.10), "h": (0.55, 0.85), "rd": (120.0, 165.0),
              "s8": (0.2, 1.4)}
    for i in range(5):
        PRIORS[f"w{i}"] = (-3.0, 1.0)
    names = ["Om", "Ob", "h", "rd", "s8"] + [f"w{i}" for i in range(5)]
    x0 = np.array([0.31, 0.049, 0.68, 147.0, 0.81] + [-0.93, -0.97, -1.2, -1.0, -1.1])
    step = np.array([0.01, 0.005, 0.01, 0.5, 0.02] + [0.06] * 5)
    out = {"sample": "Dovekie + BAO(ALL+lyaAP) + CMB + fs8", "edge_sets": {}}
    for ename, edges in EDGE_SETS.items():
        def logpost(x, edges=edges):
            d = dict(zip(names, x))
            for n in names:
                lo, hi = PRIORS[n]
                if not (lo <= d[n] <= hi):
                    return -np.inf
            Om, Ob, h, rd, s8 = d["Om"], d["Ob"], d["h"], d["rd"], d["s8"]
            wb = [d[f"w{i}"] for i in range(5)]
            bg = BackgroundBinned(Om, wb, edges)
            lnF = lambda z: lnF_binned(z, edges, wb)
            H0 = 100.0 * h
            return -0.5 * (bao_chi2(bg, H0, rd, bao) + sn_hd_chi2(bg, H0, z_d, mu_d, invC_d)
                           + cmb_compressed_chi2(Om, Ob, h, bg)
                           + fs8_chi2(Om, s8, z_f, fs8v, invC_f, lnF=lnF))
        t0 = time.time()
        samples, acc = metropolis(logpost, x0, nsteps=22000, burn=6000, step=step, seed=233)
        stride = max(1, len(samples) // 2500)
        sub = samples[::stride]
        lps = np.array([logpost(s) for s in sub])
        wi = sub[:, 5:10]
        lbl = [f"[{edges[i]:g},{edges[i+1]:g})" for i in range(5)]
        entry = {"edges": edges, "chi2_map": float(-2 * lps.max()), "acceptance": float(acc),
                 "w_bins_mean": wi.mean(axis=0).tolist(), "w_bins_std": wi.std(axis=0).tolist(),
                 "seconds": round(time.time() - t0, 1)}
        out["edge_sets"][ename] = entry
        print(f"  [{ename}] chi2_map={entry['chi2_map']:.2f}: "
              + " | ".join(f"{l}: {m:+.3f}+-{s:.3f}" for l, m, s in
                          zip(lbl, wi.mean(axis=0), wi.std(axis=0))), flush=True)
        np.save(os.path.join(RESULTS, f"chain_s23_edges_{ename}.npy"),
                np.column_stack([sub, lps]))
    with open(os.path.join(RESULTS, "edges_results.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("EDGES DONE", flush=True)


def dipole_mode():
    print("loading data ...", flush=True)
    bao = concat_rows(cosmo.load_desi_bao(which="ALL")[0], LYA_AP)
    z_p, mB_p, invC_p = cosmo.load_pantheonplus()
    raw = np.genfromtxt(os.path.join(cosmo.DATA, "sn_data", "PantheonPlus", "Pantheon+SH0ES.dat"),
                        names=True, dtype=None, encoding="utf-8")
    RA = np.deg2rad(np.asarray(raw["RA"], float))
    DEC = np.deg2rad(np.asarray(raw["DEC"], float))
    z = z_p

    # LCDM geometry from BAO+LyaAP+CMB (as in session16)
    from scipy.optimize import minimize

    def lp(x):
        Om, Ob, h, rd = x
        if not (0.05 < Om < 0.9 and 0.005 < Ob < 0.10 and 0.55 < h < 0.85 and 120 < rd < 165):
            return -np.inf
        bg = Background(Om)
        return -0.5 * (bao_chi2(bg, 100 * h, rd, bao) + cmb_compressed_chi2(Om, Ob, h, bg))

    r0 = minimize(lambda x: -lp(x), [0.31, 0.049, 0.68, 147.0], method="Nelder-Mead",
                  options={"maxiter": 3000})
    bg = Background(r0.x[0])
    H0 = 100 * r0.x[2]
    print(f"  geometry LCDM H0={H0:.2f}", flush=True)

    # M-marginalized residuals per SN (analytic, same as sn_mB_chi2)
    r_all = mB_p - bg.mu(z, H0)
    one = np.ones(len(z))
    C1 = invC_p @ one
    r_all = r_all - (r_all @ C1) / (one @ C1)   # M marginalized

    # design: dipole (3 params) + quadrupole (5 params, l=2 spherical harmonics real)
    def design(ra, dec):
        x = np.sin(dec) * np.cos(ra)
        y = np.sin(dec) * np.sin(ra)
        zz = np.cos(dec)
        l2 = [np.sqrt(3) * x * y, np.sqrt(3) * y * zz, np.sqrt(3) * (zz ** 2 - 1.0 / 3.0),
              np.sqrt(3) * x * zz, np.sqrt(5.0 / 4.0) * (x ** 2 - y ** 2)]
        return np.column_stack([np.ones_like(ra), x, y, zz] + l2)

    out = {"geometry_H0": float(H0), "tests": {}}
    for label, mask in [("full", np.ones(len(z), bool)), ("z<0.1", z < 0.1), ("z>0.1", z >= 0.1)]:
        ii = np.where(mask)[0]
        C = invC_p[np.ix_(ii, ii)]
        X = design(RA[ii], DEC[ii])
        # generalized least squares with covariance
        CiX = np.linalg.solve(C, X)
        A = X.T @ CiX
        b = X.T @ np.linalg.solve(C, r_all[ii])
        beta = np.linalg.solve(A, b)
        covb = np.linalg.inv(A)
        names_d = ["const", "dx", "dy", "dz", "q_xy", "q_yz", "q_zz", "q_xz", "q_xx-yy"]
        # dipole amplitude and significance
        dip = np.linalg.norm(beta[1:4])
        sig_dip = np.sqrt(np.trace(covb[1:4, 1:4]))   # sigma on |d| (approx, for small d)
        # chi2 improvement of dipole+quadrupole vs constant
        r0v = r_all[ii]
        chi2_0 = r0v @ np.linalg.solve(C, r0v)
        chi2_m = chi2_0 - b @ beta
        entry = {"n": int(len(ii)), "dipole_vec": beta[1:4].tolist(),
                 "dipole_amp_mag": float(dip), "dipole_amp_sigma_approx": float(sig_dip),
                 "dipole_tension": float(dip / sig_dip),
                 "chi2_const": float(chi2_0), "chi2_dipquad": float(chi2_m),
                 "dchi2_8dof": float(chi2_m - chi2_0),
                 "params": dict(zip(names_d, map(float, beta)))}
        out["tests"][label] = entry
        print(f"  [{label}] n={len(ii)}: dipole amp = {dip:.4f}+-{sig_dip:.4f} mag "
              f"({dip/sig_dip:.2f}s), dchi2(8 dof) = {chi2_m-chi2_0:+.2f}", flush=True)
    with open(os.path.join(RESULTS, "dipole_results.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("DIPOLE DONE", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "dipole"
    if mode == "edges":
        edges_mode()
    else:
        dipole_mode()
