""""session29.py — REDISTRIBUTION, derived from first principles (not a toy fit).

Computes the universal scale where boundary/curvature vacuum energy equals the
cosmological-constant density, and the two independent killers of macroscopic
redistribution:
  (A) passivity/zero-sum of the ground state (reset pays the binding energy);
  (B) MATERIAL COST — any plate that carries a Casimir energy equal to rho_Lambda*c^2
      over volume V must itself have mass-energy larger by the factor
          R = 2 rho_m c^2 t d^3 / K ,  K = pi^2 hbar c / 240 ,
      minimized over gap d and plate thickness t. We scan d from atomic to 100 um.
Also: black-hole curvature matching (horizon radius where curvature vacuum
polarization ~ rho_Lambda c^2), and the dynamical-Casimir statement that real
photons require drive work (energy conservation), so the vacuum cannot be a source.

Usage: python experiments/session29.py redist
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


def redist_mode():
    H0 = 67.4
    Om = 0.315
    H0_si = H0 * 1e3 / MPC
    rhoL = 3 * H0_si ** 2 * (1 - Om) / (8 * np.pi * G) * C * C   # rho_Lambda c^2 [J/m^3]
    K = np.pi ** 2 * HBAR * C / 240.0                              # Casimir coeff [J m]

    d_eq = (K / rhoL) ** 0.25                                      # parallel-plate full-cancel gap
    L_curv = (HBAR * C / rhoL) ** 0.25                             # curvature radius, coeff=1
    out = {
        "rho_Lambda_c2_Jm3": float(rhoL),
        "Casimir_coeff_K_Jm": float(K),
        "universal_scales": {
            "parallel_plate_gap_full_cancel_m": float(d_eq),
            "curvature_radius_coeff1_m": float(L_curv),
            "note": "order 40-90 um: the ONLY length where boundary/curvature vacuum energy "
                    "density equals Lambda; macroscopic volumes average it to ~0",
        },
        "black_hole_curvature_match": {},
        "gap_scan": [],
        "material_cost_min": {},
        "dynamical_casimir": {},
    }

    # black hole whose horizon curvature scale L_curv matches Lambda density
    Rs = L_curv
    M = Rs * C * C / (2 * G)
    P_haw = HBAR * C ** 6 / (15360 * np.pi * G ** 2 * M ** 2)
    tau_s = 5120 * np.pi * G ** 2 * M ** 3 / (HBAR * C ** 4)
    out["black_hole_curvature_match"] = {
        "Rs_m": float(Rs), "mass_kg": float(M), "mass_vs_moon": float(M / 7.35e22),
        "hawking_power_W": float(P_haw), "lifetime_yr": float(tau_s / 3.156e7),
        "note": "a horizon of ~90 um (mass ~0.8 Moon) is where curvature vacuum "
                "polarization ~ Lambda; Hawking power ~1e-13 W, lifetime ~1e45 yr — a "
                "redistributor in principle, useless in practice",
    }

    rho_m = 2700.0          # aluminum [kg/m^3]
    rho_m_c2 = rho_m * C * C
    t_atom = 3.0e-10        # ~one atomic layer [m]
    for d in [3e-10, 1e-9, 3e-9, 1e-8, 1e-7, 1e-6, 1e-5, d_eq, 1e-4]:
        P = K / d ** 4                       # Casimir pressure [Pa]
        rho_C = K / d ** 4                   # |energy density| [J/m^3]
        f = rhoL / rho_C                     # fill factor of gap volume to average-cancel Lambda
        A_per_V = f / d                      # plate area per unit shielded volume [1/m]
        R = 2 * rho_m_c2 * t_atom * d ** 3 / K   # plate mass-energy / canceled Lambda energy
        out["gap_scan"].append({
            "gap_m": float(d), "pressure_Pa": float(P), "rho_C_Jm3": float(rho_C),
            "fill_factor": float(f), "plate_area_per_m3": float(A_per_V),
            "plate_mass_energy_ratio_R": float(R),
        })
    # absolute minimum of R over d>=t>=t_atom (Casimir form extrapolated to atomic contact)
    R_min = 2 * rho_m_c2 * t_atom * t_atom ** 3 / K
    out["material_cost_min"] = {
        "R_at_atomic_contact": float(R_min),
        "meaning": "even with one-atom-thick plates at atomic gap, plate rest-energy is "
                   ">= ~1e9 x the Lambda energy being canceled; for any mechanically "
                   "stable gap (>=10 nm, pressure <= ~1 atm) R >= 1e14",
        "conclusion": "material cost kills boundary redistribution independently of the "
                      "zero-sum reset cost (passivity)",
    }

    # dynamical Casimir: a moving mirror converts virtual pairs to real photons, but the
    # drive supplies the energy. One 1 eV photon needs >=1 eV of work from the actuator;
    # vacuum contributes the virtual intermediate state only (energy conservation).
    eV = 1.602176634e-19
    out["dynamical_casimir"] = {
        "one_1eV_photon_min_drive_work_J": float(eV),
        "net_extractable_from_vacuum_J": 0.0,
        "statement": "photon creation by moving boundaries is real (dynamical Casimir) but "
                     "the vacuum is not the energy source; the actuator is. No redistribution "
                     "of Lambda follows.",
    }

    with open(os.path.join(RESULTS, "redistribution_theorem.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("rho_Lambda_c2 =", rhoL)
    print("K =", K)
    print("d_eq =", d_eq, "m ; L_curv =", L_curv, "m")
    print("BH match: Rs=", Rs, "M=", M, "kg =", M / 7.35e22, "Moon ; P_haw=", P_haw, "W")
    print("gap scan (gap_m, P_Pa, fill_factor, A_per_m3, R):")
    for r in out["gap_scan"]:
        print("  %.3e  %.3e  %.3e  %.3e  %.3e" % (r["gap_m"], r["pressure_Pa"],
                                                    r["fill_factor"], r["plate_area_per_m3"],
                                                    r["plate_mass_energy_ratio_R"]))
    print("R_min at atomic contact =", R_min)
    print("REDIST DONE", flush=True)


if __name__ == "__main__":
    redist_mode()
