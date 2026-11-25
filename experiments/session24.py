"""session24.py — the harnessing ceilings, computed (N12/N13 support) + DESI-chain dip
cross-check.

Mode 'caps': every number behind the "harness or redistribute dark energy" concepts:
  rho_Lam c^2, continuous cap H0*rho*c^2, phantom-turbine power 3H(z)|1+w|rho c^2
  (posterior from the s22 cert5 Dovekie chain at z=1, w in bin [0.75,1.5)),
  Casimir gap d where |rho_Casimir| = rho_Lam, plate area to shield 1 m^3,
  Unruh bath (a for T=1 K, flux), de Sitter horizon (T_dS, total power),
  Alcubierre rho at R=100 m vs best Casimir density, expansion-piston power.

Mode 'desichain': w(z) posterior from DESI's OWN published w0wa chains
(w(z) = w0 + wa z/(1+z)) — does the collaboration's chain see the z~1 dip?

Usage: python experiments/session24.py caps|desichain
"""
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

HBAR = 1.054571817e-34
KB = 1.380649e-23
C = 299792458.0
G = 6.67430e-11


def caps_mode():
    H0 = 67.4  # km/s/Mpc, our geometry fit
    h = H0 / 100.0
    Om = 0.315
    OL = 1 - Om
    H0_si = H0 * 1e3 / 3.08567758e22          # 1/s
    rhoL = 3 * H0_si ** 2 * OL / (8 * np.pi * G)   # kg/m^3
    rhoL_c2 = rhoL * C * C                          # J/m^3
    out = {"H0_kmsMpc": H0, "rho_Lambda_c2_Jm3": float(rhoL_c2),
           "continuous_cap_Wm3_H0rho": float(H0_si * rhoL_c2),
           "one_shot_full_Jm3": float(rhoL_c2)}

    # phantom-turbine: posterior of 3 H(z=1) |1+w| rho_DE(z=1) c^2 from s22 cert5 chain
    ch = np.load(os.path.join(RESULTS, "chain_s22_cert5_dovekie.npy"))
    Om_s, Ob_s, h_s, rd_s, s8_s = ch[:, 0], ch[:, 2], ch[:, 2], ch[:, 3], ch[:, 4]
    h_s = ch[:, 2]
    w2 = ch[:, 7]  # bin [0.75,1.5)
    z1 = 1.0
    Ez1 = np.sqrt(Om_s * (1 + z1) ** 3 + (1 - Om_s) * np.exp(3 * (1 + w2) * np.log(1 + z1)))
    Hz = 100 * h_s * Ez1 * 1e3 / 3.08567758e22
    rhoDE_z1_c2 = rhoL_c2 * np.exp(3 * (1 + w2) * np.log(1 + z1))
    P_turb = 3 * Hz * np.abs(1 + w2) * rhoDE_z1_c2   # W/m^3 (rho grows if w<-1)
    out["phantom_turbine"] = {
        "w_bin_0.75_1.5_mean": float(w2.mean()), "w_bin_std": float(w2.std()),
        "P_Wm3_mean": float(P_turb.mean()), "P_Wm3_p16": float(np.percentile(P_turb, 16)),
        "P_Wm3_p84": float(np.percentile(P_turb, 84)),
        "note": "power density if the z~1 phantom dip is real: d(rho c^2)/dt = -3H(1+w) rho c^2 > 0"}

    # Casimir: rho(d) = -pi^2 hbar c / (240 d^4)
    pref = np.pi ** 2 * HBAR * C / 240.0
    d_eq = (pref / rhoL_c2) ** 0.25
    out["casimir"] = {
        "rho_at_10nm_Jm3": float(-pref / (1e-8) ** 4),
        "gap_d_for_full_cancel_m": float(d_eq),
        "pressure_at_d_eq_Pa": float(rhoL_c2),   # equals rho_Lam c^2
        "plate_area_for_1m3_shield_m2": float(1.0 / d_eq),
        "one_shot_per_collapse_J_per_m3": float(rhoL_c2),
        "note": "resetting the cavity costs the same Casimir binding energy -> net zero; "
                "only expansion (3H rho c^2) is net-positive"}

    # Unruh bath: T = hbar a / (2 pi c kB)
    a_1K = 2 * np.pi * C * KB * 1.0 / HBAR
    sigSB = 5.670374419e-8
    out["unruh_bath"] = {"a_for_T1K_ms2": float(a_1K),
                         "flux_1K_1m2_W": float(sigSB * 1.0 ** 4),
                         "note": "radiation is paid for by the accelerating agent (2nd law); no gain"}

    # de Sitter horizon
    T_dS = HBAR * H0_si / (2 * np.pi * KB)
    A_dS = 4 * np.pi * (C / H0_si) ** 2
    out["desitter_horizon"] = {"T_dS_K": float(T_dS), "flux_Wm2": float(sigSB * T_dS ** 4),
                               "total_power_W": float(sigSB * T_dS ** 4 * A_dS)}

    # Alcubierre: rho ~ -c^4/(8 pi G R^2)
    R = 100.0
    rho_alc = C ** 4 / (8 * np.pi * G * R * R)
    out["alcubierre"] = {"R_m": R, "rho_needed_Jm3": float(rho_alc),
                         "vs_casimir_10nm_ratio": float(rho_alc / (pref / 1e-16)),
                         "note": "warp needs ~1e33x more negative density than the best Casimir gap"}

    # expansion piston: dW = rho c^2 dV, comoving dV/V = 3H dt
    out["expansion_piston_Wm3"] = float(3 * H0_si * rhoL_c2)

    with open(os.path.join(RESULTS, "harness_caps.json"), "w") as f:
        json.dump(out, f, indent=2)
    for k, v in out.items():
        print(f"  {k}: {v}", flush=True)
    print("CAPS DONE", flush=True)


def desichain_mode():
    d = os.path.join(cosmo.DATA, "desi_y3_chains")
    out = {}
    zg = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.33]
    for f in sorted(glob.glob(os.path.join(d, "w0wa__*.txt"))):
        name = os.path.basename(f).replace("w0wa__", "").replace(".txt", "")
        try:
            lines = open(f).readlines()
            hdr = [l for l in lines if l.startswith("#")][-1]
            toks = hdr.lstrip("#").split()
            iw = toks.index("w")
            iwa = toks.index("wa")
            iwt = toks.index("weight")
            data = np.loadtxt(f)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            wgt = data[:, iwt]
            if wgt.sum() <= 0 or not np.isfinite(wgt).all():
                raise ValueError("bad weights")
            w0 = data[:, iw]
            wa = data[:, iwa]
        except Exception as e:
            print(f"  [{name}] skipped ({e})", flush=True)
            continue
        entry = {"n": int(len(w0)), "w0": [float(np.average(w0, weights=wgt)), float(np.std(w0))],
                 "wa": [float(np.average(wa, weights=wgt)), float(np.std(wa))], "wz": {}}
        for z in zg:
            wz = w0 + wa * z / (1 + z)
            m = np.average(wz, weights=wgt)
            s = np.sqrt(np.average((wz - m) ** 2, weights=wgt))
            entry["wz"][str(z)] = [float(m), float(s), float((m + 1) / s)]  # tension vs -1
        out[name] = entry
        print(f"  [{name}] w0={entry['w0'][0]:+.3f}+-{entry['w0'][1]:.3f} "
              f"wa={entry['wa'][0]:+.3f}+-{entry['wa'][1]:.3f} | "
              + " ".join(f"w({z})={entry['wz'][str(z)][0]:+.3f}+-{entry['wz'][str(z)][1]:.3f}"
                         f"({entry['wz'][str(z)][2]:+.1f}s vs -1)" for z in zg), flush=True)
    with open(os.path.join(RESULTS, "desichain_wz.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("DESICHAIN DONE", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "caps"
    if mode == "desichain":
        desichain_mode()
    else:
        caps_mode()
