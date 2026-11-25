""""session31.py — concept designs for the three mechanisms + levitation + phantom dial.

Computes the engineering envelopes for:
  G1  Lambda-pressure Casimir gauge at d* (sensitivity, patch-potential budget)
  S2  vacuum spring / Casimir levitator load capacity vs gap
  E3  expansion tether force at fixed proper separation + 1 W volumes
  P5  phantom-dip rudder 1 W volume (if DR3 dip >=5s)
Saved to results/concept_designs.json. Fast, no MCMC.

Usage: python experiments/session31.py concepts
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")
HBAR = 1.054571817e-34
C = 299792458.0
G = 6.67430e-11
MPC = 3.08567758e22
EPS0 = 8.8541878128e-12
KB = 1.380649e-23
AU = 1.495978707e11


def concepts_mode():
    H0 = 67.4
    Om = 0.315
    H0_si = H0 * 1e3 / MPC
    rhoL = 3 * H0_si ** 2 * (1 - Om) / (8 * np.pi * G) * C * C
    K = np.pi ** 2 * HBAR * C / 240.0
    d_star = (K / rhoL) ** 0.25
    P_star = rhoL
    out = {"rho_Lambda_c2": float(rhoL), "d_star_m": float(d_star), "P_star_Pa": float(P_star)}

    # G1 gauge sensitivity: sigma_P = F_noise/A ; relative = sigma_P/P_star
    gauge = []
    for A, Fn in [(1e-6, 1e-15), (1e-4, 1e-15), (1.0, 1e-15), (1.0, 1e-12), (1e6, 1e-9)]:
        sP = Fn / A
        rel = sP / P_star
        dd = d_star * rel / 4.0            # since P ~ d^-4, rel_P = 4 rel_d
        # patch potential whose electrostatic pressure equals sigma_P: sP = eps0 V^2/(2 d^2)
        V = np.sqrt(2 * d_star ** 2 * sP / EPS0)
        gauge.append({"area_m2": A, "force_noise_N": Fn, "sigma_P_Pa": float(sP),
                      "sigma_rho_over_rho": float(rel), "gap_resolution_m": float(dd),
                      "max_patch_potential_V_for_this_noise": float(V)})
    out["G1_gauge"] = gauge
    # drift over 1 yr if w != -1 today: drho/rho = 3(1+w) H dt
    yr = 3.15576e7
    out["G1_drift"] = {f"w={w}": float(3 * (1 + w) * H0_si * yr) for w in [-1.0, -0.9, -1.27]}

    # S2 levitator: load per area = P/g (attractive Casimir; repulsive Lifshitz same order)
    g0 = 9.81
    lev = []
    for d in [1e-8, 1e-7, 1e-6, d_star, 1e-4]:
        P = K / d ** 4
        load = P / g0                     # kg/m^2
        lev.append({"gap_m": float(d), "P_Pa": float(P), "load_kg_per_m2": float(load),
                    "area_for_1kg_m2": float(1.0 / load)})
    out["S2_levitator"] = lev
    out["S2_note"] = ("attractive Casimir is unstable (force falls as gap opens); need "
                      "repulsive Lifshitz fluid/metamaterial or active feedback. d* is the "
                      "Lambda-matching sensor point (load ~5e-11 kg/m2); real load needs "
                      "nm gaps where Casimir >> Lambda, so dark energy is irrelevant to the force.")

    # E3 expansion tether holding fixed proper separation r in de Sitter: F = m * 8 pi G rhoL r/(3 c^2)
    def tether_F(m, r):
        return m * 8 * np.pi * G * rhoL * r / (3 * C * C)
    out["E3_tether"] = {f"r={r:g}": float(tether_F(1.0, r))
                        for r in [1.0, 1e3, AU, MPC]}
    pwr_exp = 3 * H0_si * rhoL
    P_phantom = 8.082320983263925e-28
    out["E3_P5_1W_volumes"] = {
        "expansion_W_per_m3": float(pwr_exp), "V_for_1W_m3": float(1 / pwr_exp),
        "side_for_1W_m": float((1 / pwr_exp) ** (1 / 3)),
        "phantom_W_per_m3": P_phantom, "phantom_V_for_1W_m3": float(1 / P_phantom),
        "phantom_side_for_1W_m": float((1 / P_phantom) ** (1 / 3))}

    with open(os.path.join(RESULTS, "concept_designs.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("d_star =", d_star, "P_star =", P_star)
    print("G1 gauge rows (A,Fn,sigmaP,rel,dd,Vpatch):")
    for r in gauge:
        print("  %.0e %.0e %.3e %.3e %.3e %.3e" % (r["area_m2"], r["force_noise_N"], r["sigma_P_Pa"],
                                                     r["sigma_rho_over_rho"], r["gap_resolution_m"],
                                                     r["max_patch_potential_V_for_this_noise"]))
    print("G1 drift/yr:", out["G1_drift"])
    print("S2 levitator (gap,P,load kg/m2,area 1kg):")
    for r in lev:
        print("  %.3e %.3e %.3e %.3e" % (r["gap_m"], r["P_Pa"], r["load_kg_per_m2"], r["area_for_1kg_m2"]))
    print("E3 tether F (1 kg):", out["E3_tether"])
    print("1W volumes:", out["E3_P5_1W_volumes"])
    print("CONCEPTS DONE", flush=True)


if __name__ == "__main__":
    concepts_mode()
