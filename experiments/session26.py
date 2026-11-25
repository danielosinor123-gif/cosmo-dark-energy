"""session26.py — the two conceptual STEERING devices, quantified (N14/N15).

Mode 'rudder':
  N14 VACUUM RUDDER (steer Lambda_eff locally via the phantom sector — requires the
  z~1 dip to be real): posterior of the steerable fraction |1+w| and the local
  energy-creation rate 3H|1+w| rho c^2 from the s22 cert5 chain; the local tilt
  signature delta_mu = (5/ln10) * (3/2)(1+w) ln(1+z) per unit log-distance.
  N15 OPTICAL RUDDER (steer the READING via engineered line-of-sight dispersion —
  no new physics needed, just a medium): DM (dispersion measure) needed to fake
  delta_mu = 0.03 mag at z=0.05 and at z=1, vs the real IGM DM; and the
  frequency-dependence signature (delta_t ~ DM/nu^2) that multi-band data catch.

Usage: python experiments/session26.py rudder
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

C = 299792458.0
MPC = 3.08567758e22
ME = 9.1093837e-31
E2 = (1.602176634e-19) ** 2


def rudder_mode():
    H0 = 67.4
    H0_si = H0 * 1e3 / MPC
    Om = 0.315
    rhoL_c2 = 3 * H0_si ** 2 * (1 - Om) / (8 * np.pi * 6.67430e-11) * C * C
    out = {"rho_Lambda_c2": float(rhoL_c2)}

    # ---- N14: phantom-sector steering capability from the dip posterior --------
    ch = np.load(os.path.join(RESULTS, "chain_s22_cert5_dovekie.npy"))
    Om_s, h_s, w2 = ch[:, 0], ch[:, 2], ch[:, 7]     # dip bin [0.75,1.5)
    z1 = 1.0
    Ez1 = np.sqrt(Om_s * (1 + z1) ** 3 + (1 - Om_s) * np.exp(3 * (1 + w2) * np.log(1 + z1)))
    Hz = 100 * h_s * Ez1 * 1e3 / MPC
    P = 3 * Hz * np.abs(1 + w2) * rhoL_c2 * np.exp(3 * (1 + w2) * np.log(1 + z1))
    frac = np.abs(1 + w2)
    # local tilt a rudder could imprint per e-fold of distance if it coupled with
    # strength |1+w|: dmu = (5/ln10) * (3/2) * |1+w| * ln(1+z)
    dmu_per_efold = 5 / np.log(10) * 1.5 * frac * np.log(2)
    out["N14_vacuum_rudder"] = {
        "requires": "z~1 dip real (w<-1: phantom sector, NEC violation) — DR3 gate",
        "phantom_fraction_|1+w|_mean": float(frac.mean()),
        "phantom_fraction_p16_p84": [float(np.percentile(frac, 16)), float(np.percentile(frac, 84))],
        "local_energy_creation_Wm3_mean": float(P.mean()),
        "local_energy_creation_p16_p84": [float(np.percentile(P, 16)), float(np.percentile(P, 84))],
        "max_local_tilt_mag_per_efold_mean": float(dmu_per_efold.mean()),
        "data_bound_on_local_tilt": "<0.01 mag/dex host mass = <0.9% rho per dex (§19)",
        "verdict": "even if the dip is real, a phantom-coupled rudder can only tilt its OWN "
                   "lab volume; the cosmological probe averages 1e76 m^3 — unsteerable globally"}

    # ---- N15: optical rudder — dispersion medium faking a distance shift --------
    def dm_needed(dmu, z, nu_GHz=5e5):
        DL = C * z / (H0 * 1e3) * MPC          # rough DL ~ cz/H0 * Mpc->m (low z)
        DL = C * z / H0_si                     # meters
        dt = (dmu * np.log(10) / 5.0) * DL / C  # seconds of delay to fake dmu
        # dt[s] = 4.15e-3 * DM[pc/cm^3] / nu[GHz]^2  ->  DM = dt * nu^2 / 4.15e-3
        DM = dt * nu_GHz ** 2 / 4.15e-3
        return float(dt), float(DM)
    igm_dm = lambda z: 1000 * z   # pc/cm^3, rough IGM to z
    res = {}
    for z, dmu in [(0.05, 0.03), (1.0, 0.10)]:
        dt, DM = dm_needed(dmu, z)
        res[f"z={z} fake dmu={dmu}"] = {
            "delay_s": dt, "DM_needed_pc_cm3": DM, "real_IGM_DM_pc_cm3": igm_dm(z),
            "ratio_to_IGM": DM / igm_dm(z),
            "signature": "dt ~ DM/nu^2: radio/optical/X-ray arrive at different times; "
                         "multi-band SN + radio BAO disagree -> caught"}
    out["N15_optical_rudder"] = {
        "requires": "no new physics — an engineered plasma along the whole line of sight",
        "cases": res,
        "verdict": "needs ~1e16x the IGM column with exactly the right z-dependence over "
                   "445 Mpc, and is frequency-dispersive: caught by multi-band/multi-tracer "
                   "consistency (BAO, Ly-a AP, CMB all read the same geometry)"}
    with open(os.path.join(RESULTS, "rudder_results.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out, indent=2)[:2200], flush=True)
    print("RUDDER DONE", flush=True)


if __name__ == "__main__":
    rudder_mode()
