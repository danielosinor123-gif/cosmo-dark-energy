""""session33.py — universal-coupling levitation feasibility frontier (assistive 'invisible hand').

Question: what universal property of ALL solids can a non-contact field grab?
Answer computed here: charge/polarizability (electric), mass/density (acoustic/gravity),
magnetic susceptibility (diamagnetic — universal but weak; ferro needs iron).
Dark matter/energy are NOT usable handles: DM ~ gravitational-only/unknown, DE uniform
repulsive and cannot be localized (see REDISTRIBUTION_THEOREM.md).

Computes force density / max payload vs safe limits for gravity, electrostatic (air/oil),
DEP (air/oil), acoustic (SPL-limited), magnetic diamagnetic (B), ion wind power.
Saved to results/universal_levitation.json.

Usage: python experiments/session33.py universal
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
MU0 = 4 * np.pi * 1e-7
RHO_AIR = 1.2
C_AIR = 343.0


def universal_mode():
    H0 = 67.4
    Om = 0.315
    H0_si = H0 * 1e3 / MPC
    rhoL = 3 * H0_si ** 2 * (1 - Om) / (8 * np.pi * G) * C * C
    g = 9.81
    out = {"rho_Lambda_c2": float(rhoL), "note": "assistive invisible-hand feasibility; not a weapon; BCI is setpoint-only with hard safety limits"}

    # 1 gravity: mass overhead at distance r to give upward a=g
    out["gravity"] = {f"r={r:g} m": {"M_kg_for_g": float(g * r * r / G)} for r in [0.05, 0.1, 0.5, 1.0, 10.0]}

    # 2 electrostatic in air at breakdown, and in oil (higher eps, higher E_bd)
    E_bd_air = 3e6
    P_es_air = 0.5 * EPS0 * E_bd_air ** 2
    load_es_air = P_es_air / g
    masses = {"10 g": 0.01, "100 g": 0.1, "1 kg": 1.0, "70 kg person": 70.0}
    out["electrostatic_air"] = {
        "E_bd_V_m": E_bd_air, "P_Pa": float(P_es_air), "load_kg_m2": float(load_es_air),
        "V_for_gap": {f"d={d:g} m": float(E_bd_air * d) for d in [1e-3, 1e-2, 0.1]},
        "area_for_mass_m2": {k: float(m * g / P_es_air) for k, m in masses.items()},
        "mass_hand_pad_0.05m2": float(load_es_air * 0.05),
        "mass_panel_0.25m2": float(load_es_air * 0.25),
    }
    eps_oil, E_bd_oil = 3.0, 20e6
    P_es_oil = 0.5 * EPS0 * eps_oil * E_bd_oil ** 2
    out["electrostatic_oil"] = {"P_Pa": float(P_es_oil), "load_kg_m2": float(P_es_oil / g),
                                 "mass_hand_pad_0.05m2": float(P_es_oil / g * 0.05)}

    # 3 DEP: sphere CM=(ep-em)/(ep+2em); force density approx 0.5 eps0 em CM E^2
    def dep(eps_p, eps_m, E):
        cm = (eps_p - eps_m) / (eps_p + 2 * eps_m)
        P = 0.5 * EPS0 * eps_m * cm * E ** 2
        return float(cm), float(P)
    out["DEP"] = {
        "air_plastic_eps3": dep(3.0, 1.0, E_bd_air),
        "air_conductor_CM1": dep(1e30, 1.0, E_bd_air),
        "oil_eps2_in_eps3_REPULSIVE": dep(2.0, 3.0, E_bd_oil),
        "note": "in air all solids have eps>1 -> positive DEP = attraction to electrode; repulsion needs a higher-eps medium (negative CM) or a push mechanism (acoustic/ion)"}

    # 4 acoustic: radiation pressure on reflector P ~ p^2/(rho c^2); p from SPL
    def aco(spl):
        p = 20e-6 * 10 ** (spl / 20.0)
        I = p * p / (2 * RHO_AIR * C_AIR)
        P = 2 * I / C_AIR
        return float(p), float(I), float(P), float(P / g)
    out["acoustic"] = {f"{spl} dB": dict(zip(["p_Pa", "I_W_m2", "P_Pa", "load_kg_m2"], aco(spl)))
                       for spl in [120, 140, 160, 183]}
    # SPL needed for payload m over area A (perfect reflector, order of magnitude)
    need = {}
    for m, A in [(0.001, 1e-4), (0.01, 1e-3), (0.1, 0.01), (1.0, 0.05), (70.0, 0.1)]:
        P = m * g / A
        p = np.sqrt(P * RHO_AIR * C_AIR ** 2)
        spl = 20 * np.log10(p / 20e-6)
        need[f"m={m:g}kg A={A:g}m2"] = {"P_Pa": float(P), "SPL_dB": float(spl)}
    out["acoustic_SPL_needed"] = need

    # 5 magnetic diamagnetic (universal, weak): P=|chi| B^2/(2 mu0), chi~-1e-5
    chi = -1e-5
    out["magnetic_diamagnetic"] = {
        f"B={B:g} T": {"P_Pa": float(abs(chi) * B * B / (2 * MU0)),
                       "load_kg_m2": float(abs(chi) * B * B / (2 * MU0) / g),
                       "mass_pad_0.05m2": float(abs(chi) * B * B / (2 * MU0) / g * 0.05)}
        for B in [1.0, 10.0, 16.0]}

    # 6 ion wind: ideal power = 0.5 * thrust * exhaust velocity; real = /eta
    ve = 10.0
    out["ion_wind"] = {f"m={m:g} kg": {"thrust_N": float(m * g), "P_ideal_W": float(0.5 * m * g * ve),
                                       "P_real_eta0.2_W": float(0.5 * m * g * ve / 0.2)}
                       for m in [0.01, 0.1, 1.0, 70.0]}

    # 7 dark handles
    coeff = 8 * np.pi * G * rhoL / (3 * C * C)
    out["dark_handles"] = {
        "Lambda_accel_coeff": float(coeff), "r_for_g_m": float(g / coeff),
        "Lambda_pressure_Pa": float(rhoL), "Lambda_load_kg_m2": float(rhoL / g),
        "dark_matter": "no known non-gravitational coupling; cannot be localized/controlled; gravity-only => same weakness as gravity row",
        "verdict": "dark matter/energy are NOT usable handles; universal handle = charge/polarizability (electric) and density (acoustic/gravity)"}

    with open(os.path.join(RESULTS, "universal_levitation.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("electrostatic air:", out["electrostatic_air"])
    print("electrostatic oil:", out["electrostatic_oil"])
    print("DEP:", out["DEP"])
    print("acoustic SPL needed:", out["acoustic_SPL_needed"])
    print("magnetic:", out["magnetic_diamagnetic"])
    print("ion wind:", out["ion_wind"])
    print("gravity:", out["gravity"])
    print("dark:", out["dark_handles"])
    print("UNIVERSAL DONE", flush=True)


if __name__ == "__main__":
    universal_mode()
