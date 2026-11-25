""""session30.py — WHAT the 40-90 um redistribution channel can actually do.

Computes the applied envelope of a Lambda-matching Casimir/curvature cavity:
force/pressure at d_eq, expansion power vs volume, work between gaps, plate-cost
per joule, local gravitational signature of canceling Lambda, sensor area needed to
resolve the rho_Lambda pressure, and the (unachievable) homogeneous no-Lambda local
Hubble ratio. All numbers saved to results/redist_apps.json.

Usage: python experiments/session30.py apps
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


def apps_mode():
    H0 = 67.4
    Om = 0.315
    H0_si = H0 * 1e3 / MPC
    rhoL = 3 * H0_si ** 2 * (1 - Om) / (8 * np.pi * G) * C * C
    K = np.pi ** 2 * HBAR * C / 240.0
    d_eq = (K / rhoL) ** 0.25
    P = rhoL                      # Pa at the matching gap (by construction)
    out = {"rho_Lambda_c2_Jm3": float(rhoL), "d_eq_m": float(d_eq), "P_at_d_eq_Pa": float(P)}

    # force on areas
    out["force"] = {f"A={a:g} m2": {"F_N": float(P * a), "F_fN": float(P * a / 1e-15)}
                    for a in [1e-6, 1e-4, 1.0, 1e6]}
    # sensor area for force thresholds
    out["sensor_area_for_threshold"] = {
        "F=1 fN": float(1e-15 / P), "F=1 nN": float(1e-9 / P), "F=1 pN": float(1e-12 / P)}

    # expansion power = 3 H rho_Lambda c^2 per m^3 ; volume for 1 W
    pwr = 3 * H0_si * rhoL
    V1W = 1.0 / pwr
    out["expansion_power"] = {
        "W_per_m3": float(pwr),
        "V_for_1W_m3": float(V1W), "side_of_cube_for_1W_m": float(V1W ** (1 / 3)),
        "moon_volumes_for_1W": float(V1W / 2.2e19),
        "energy_per_m3_over_Hubble_time_J": float(H0_si * rhoL / H0_si),  # ~ rhoL
    }

    # work extractable by collapsing a gap from d_large to d_small: W/A = K/3 (1/d_small^3 - 1/d_large^3)
    gaps = [d_eq, 1e-5, 1e-6, 1e-7, 1e-8]
    wa = {}
    for d in gaps:
        wa[f"collapse_to_{d:g}"] = float(K / 3 * (1 / d ** 3 - 1 / d_eq ** 3))  # J/m^2 (d<d_eq)
    out["work_per_area_collapse_Jm2"] = wa
    d10 = 1e-8
    W_A = K / 3 * (1 / d10 ** 3 - 1 / d_eq ** 3)
    A_for_1J = 1.0 / W_A
    R_10nm = 2 * 2700 * C * C * 3e-10 * d10 ** 3 / K
    out["one_joule_from_10nm_gap_collapse"] = {
        "area_m2": float(A_for_1J), "plate_mass_energy_ratio_R": float(R_10nm),
        "material_cost_J": float(R_10nm * 1.0), "net_J": float(1.0 - R_10nm)}

    # local gravitational signature if Lambda were homogeneously canceled in a sphere r
    r = d_eq / 2
    g = 8 * np.pi * G * rhoL / (3 * C * C) * r
    out["local_gravity_of_cancelled_Lambda"] = {
        "r_m": float(r), "delta_g_ms2": float(g),
        "note": "acceleration from Lambda at the cavity radius; ~1e-40 m/s^2 — no propulsion"}
    out["homogeneous_noLambda_local_H_ratio"] = float(np.sqrt(Om))  # H_local/H if rho_Lambda=0 at z=0

    with open(os.path.join(RESULTS, "redist_apps.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("P_at_d_eq =", P, "Pa ; d_eq =", d_eq)
    print("force:", {k: v["F_N"] for k, v in out["force"].items()})
    print("sensor area F=1fN/1pN/1nN (m2):", out["sensor_area_for_threshold"])
    print("expansion W/m3 =", pwr, " V(1W)=", V1W, " side=", V1W ** (1 / 3), " Moons=", V1W / 2.2e19)
    print("W/A collapse:", wa)
    print("1 J from 10nm collapse: area=", A_for_1J, "R=", R_10nm, "net=", 1 - R_10nm)
    print("delta_g at r=d_eq/2 =", g, "m/s^2 ; homogeneous no-Lambda H ratio =", np.sqrt(Om))
    print("APPS DONE", flush=True)


if __name__ == "__main__":
    apps_mode()
