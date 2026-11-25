"""session25.py — "can we manipulate the dark-energy reading?" quantified.

Mode 'manip': what it would take to FAKE each feature we read, as a local/global
change of Lambda, vs what our data/theorems allow:
  (1) fake the +0.03 mag low-z step: needed delta rho_Lam / rho_Lam inside z<0.1
      (via delta H0 / H0 ~ -delta mu / (5/ln10) at z~0.05, and delta rho/rho ~ 2 delta H/H)
      + the volume and total energy of the z<0.1 sphere.
  (2) fake the w0wa hint (wa=-0.5): needed fractional change of rho_DE(z) vs z=0
      at z~0.5-1 (w<-1 makes rho grow: delta rho/rho = 3|1+w| ln(1+z)).
  (3) what our N4 host-mass bound allows: <0.01 mag/dex -> delta rho/rho bound.
  (4) Casimir local cancel: gap 40 um, volume = gap*area; to tilt a z<0.1-scale
      reading you'd need the effect over ~ (400 Mpc)^3 — plate area implied.
  (5) global change of Lambda = boundary condition of the whole universe
      (unimodular): not a local operation; would shift ALL redshifts, not just z<0.1 —
      excluded by the z>0.1 data agreeing with Lambda (§8: z-cut kills the signal).

Usage: python experiments/session25.py manip
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

HBAR = 1.054571817e-34
C = 299792458.0
G = 6.67430e-11
MPC = 3.08567758e22


def manip_mode():
    H0 = 67.4
    h = H0 / 100.0
    Om = 0.315
    H0_si = H0 * 1e3 / MPC
    rhoL_c2 = 3 * H0_si ** 2 * (1 - Om) / (8 * np.pi * G) * C * C
    out = {"rho_Lambda_c2": float(rhoL_c2)}

    # (1) fake the +0.03 mag step at z~0.05
    dmu = 0.03
    dDL_DL = dmu * np.log(10) / 5.0            # 1.38%
    dH_H = -dDL_DL                              # distance ~ 1/H locally
    drho_rho = -2 * dH_H                        # rho_Lam ~ H0^2 (flat, Lambda-dominated low-z approx)
    z_sphere = 0.1
    R = C * z_sphere / H0_si                    # ~ c z / H0
    V = 4 / 3 * np.pi * R ** 3
    E = drho_rho * rhoL_c2 * V
    out["fake_lowz_step"] = {
        "delta_mu_mag": dmu, "delta_DL_over_DL": float(dDL_DL),
        "delta_H_over_H": float(dH_H), "delta_rhoLambda_over_rhoLambda": float(drho_rho),
        "z<0.1_sphere_radius_Mpc": float(R / MPC), "volume_m3": float(V),
        "energy_that_must_be_injected_J": float(E),
        "note": "a ~2.8% local change of Lambda inside a 430-Mpc sphere, with NO mechanism "
                "(equivalence principle: Lambda cannot be concentrated); N4 bound <0.01 mag/dex "
                "in host mass is far tighter"}

    # (2) fake wa=-0.5 at z~1: rho_DE(z)/rho_DE(0) = (1+z)^{3(1+w)} ; w=-1.5 -> growth
    w_fake = -1.5
    zf = 1.0
    ratio = (1 + zf) ** (3 * (1 + w_fake))        # rho(z=1)/rho0 ; w=-1 gives 1
    grow = 1.0 / ratio - 1.0                      # growth from z=1 to now vs const
    out["fake_wa"] = {"w_needed": w_fake, "rho_z1_over_rho0": float(ratio),
                      "growth_z1_to_now_vs_const": float(grow),
                      "note": "phantom w=-1.5 means the vacuum density grows by a factor "
                              "2.8 from z=1 to today; a GLOBAL field change visible at ALL z "
                              "— our data cap |w+1| in the dip bin at 0.27+-0.16 (1.7s)"}

    # (3) N4 bound: <0.01 mag per dex of host mass -> delta rho/rho
    dmu_dex = 0.01
    out["n4_bound"] = {"delta_mu_per_dex_mag": dmu_dex,
                       "delta_rho_over_rho_per_dex": float(2 * dmu_dex * np.log(10) / 5),
                       "note": "astrophysical environments cannot tilt Lambda beyond ~0.9%/dex (§19)"}

    # (4) Casimir: to cancel Lambda over the z<0.1 sphere you'd need plates covering
    #     V/d area at d=40 um
    pref = np.pi ** 2 * HBAR * C / 240.0
    d_eq = (pref / rhoL_c2) ** 0.25
    A_need = V / d_eq
    out["casimir_over_sphere"] = {"gap_m": float(d_eq), "plate_area_m2": float(A_need),
                                  "plate_area_vs_observable_universe_area": float(A_need / (4 * np.pi * (C / H0_si) ** 2)),
                                  "note": "shielding a cosmological volume is absurd; lab-scale "
                                          "cavities only affect the gap volume (zero-sum, N12)"}

    # (5) analytic steering: anchor swap numbers (from §15) — no physics change needed
    out["analytic_steering"] = {
        "Planck_anchor_sigma": 4.36, "Cepheid_anchor_sigma": 0.42, "TRGB": 1.65, "JAGB": 1.80,
        "z>0.1_cut": "-6.6 -> -0.0", "deltaM_step": "1.72s -> 0.00s",
        "note": "the '3-4 sigma evolving DE' reading is steered ENTIRELY by analysis choices "
                "(anchor, z-cut, step) — which is why dr3_plan.json pre-registers the cuts"}

    with open(os.path.join(RESULTS, "manip_results.json"), "w") as f:
        json.dump(out, f, indent=2)
    for k, v in out.items():
        print(f"  {k}: {v}", flush=True)
    print("MANIP DONE", flush=True)


if __name__ == "__main__":
    manip_mode()
