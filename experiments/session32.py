""""session32.py — dissect 'move bigger objects' scaling. Fast numbers, no MCMC.

Computes, vs gap d and target mass m:
  load = P/g = K/(g d^4); area needed A = m g d^4 / K
  gap tolerance for +/-10% force: dd = d/40
  plate mass for two plates thickness t over area A
  actuator count for cell size s: N = A/s^2
  collapse energy W = K A/(3 d^3)
Plus: vacuum-balloon volume to lift 1 kg (heuristic buoyancy in Lambda),
radius where Lambda acceleration equals g (can never lift), theoretical horizon-shell
power, and material-cost R at d*.
Saved to results/bigger_objects.json.

Usage: python experiments/session32.py scale
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


def scale_mode():
    H0 = 67.4
    Om = 0.315
    H0_si = H0 * 1e3 / MPC
    rhoL = 3 * H0_si ** 2 * (1 - Om) / (8 * np.pi * G) * C * C
    K = np.pi ** 2 * HBAR * C / 240.0
    g = 9.81
    rho_m = 2700.0
    t = 1e-6
    out = {"rho_Lambda_c2": float(rhoL), "K": float(K), "d_star_m": float((K / rhoL) ** 0.25)}

    gaps = [1e-8, 1e-7, 1e-6, 1e-5, (K / rhoL) ** 0.25, 1e-4]
    masses = {"1 mg": 1e-6, "1 g": 1e-3, "1 kg": 1.0, "1 ton": 1e3, "1000 tons": 1e6}
    rows = []
    for d in gaps:
        P = K / d ** 4
        load = P / g
        r = {"gap_m": float(d), "P_Pa": float(P), "load_kg_per_m2": float(load),
             "gap_tol_for_10pct_force_m": float(d / 40.0)}
        for name, m in masses.items():
            A = m * g / P
            r[f"area_for_{name}_m2"] = float(A)
            r[f"plate_mass_for_{name}_kg_t1um"] = float(rho_m * 2 * A * t)
            r[f"collapse_energy_for_{name}_J"] = float(K * A / (3 * d ** 3))
            r[f"actuators_1mm_cells_for_{name}"] = float(A / 1e-6)
        rows.append(r)
    out["gap_mass_scaling"] = rows

    # vacuum balloon heuristic: cancel Lambda inside volume V -> density contrast rhoL/c^2
    drho = rhoL / C ** 2
    V_1kg = 1.0 / (drho * g)
    out["vacuum_balloon"] = {
        "density_contrast_kg_m3": float(drho),
        "volume_to_lift_1kg_m3": float(V_1kg),
        "moon_volumes_for_1kg": float(V_1kg / 2.2e19),
        "note": "heuristic buoyancy in the Lambda background; shell material cost R>=1e9-1e14 kills it"}

    # radius where Lambda acceleration equals g (mass-independent)
    coeff = 8 * np.pi * G * rhoL / (3 * C * C)   # m/s^2 per meter of r
    r_lift = g / coeff
    out["lambda_acceleration"] = {
        "coeff_s^-2_per_m": float(coeff), "r_to_reach_g_m": float(r_lift),
        "r_lift_over_observable_radius": float(r_lift / (C / H0_si)),
        "note": "dark-energy gravity can never lift 1 kg at any accessible distance"}

    # theoretical horizon-shell expansion power (impossible boundary)
    R_H = C / H0_si
    V_obs = 4 / 3 * np.pi * R_H ** 3
    out["horizon_shell_theoretical"] = {
        "R_H_m": float(R_H), "V_m3": float(V_obs),
        "power_W_if_enclosed": float(3 * H0_si * rhoL * V_obs),
        "note": "requires a piston shell at the horizon; not physical/causal"}

    d_star = (K / rhoL) ** 0.25
    out["material_cost_R_at_dstar"] = float(2 * rho_m * C * C * 3e-10 * d_star ** 3 / K)

    with open(os.path.join(RESULTS, "bigger_objects.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("d_star um =", d_star * 1e6)
    for r in rows:
        print("gap=%.3e P=%.3e load=%.3e kg/m2 tol=%.2e m | 1kg area=%.3e platekg=%.3e W=%.3e N1mm=%.3e"
              % (r["gap_m"], r["P_Pa"], r["load_kg_per_m2"], r["gap_tol_for_10pct_force_m"],
                 r["area_for_1 kg_m2"], r["plate_mass_for_1 kg_kg_t1um"], r["collapse_energy_for_1 kg_J"],
                 r["actuators_1mm_cells_for_1 kg"]))
    print("vacuum balloon:", out["vacuum_balloon"])
    print("lambda accel:", out["lambda_acceleration"])
    print("horizon shell:", out["horizon_shell_theoretical"])
    print("R at dstar:", out["material_cost_R_at_dstar"])
    print("SCALE DONE", flush=True)


if __name__ == "__main__":
    scale_mode()
