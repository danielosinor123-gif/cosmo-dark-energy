"""session11.py — harnessing feasibility: cosmological bounds -> engineering numbers.

For each candidate mechanism that could release/redirect vacuum energy, take the bound
established in our fits (RESULTS.md) and convert it into:
  - max releasable energy density  u_max = f * rho_vac * c^2        [J/m^3]
  - energy in a galaxy-scale volume V = (10 kpc)^3
  - power density if released on timescale tau (Hubble time, Gyr, yr)
  - comparison vs solar constant 1361 W/m^2 and vs global power use 6e20 J/yr
Outputs results/harness_budget.json + prints a markdown table.

Usage: python experiments/session11.py budget
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

C = 299792458.0            # m/s
RHO_CRIT = 8.6e-27         # kg/m^3 (H0=70)
OMEGA_L = 0.68
RHO_VAC = OMEGA_L * RHO_CRIT          # kg/m^3
U_VAC = RHO_VAC * C ** 2              # J/m^3
SOLAR = 1361.0            # W/m^2
YEAR = 3.156e7            # s
HUBBLE_T = 1.4e10 * YEAR  # s
GAL_V = (10 * 3.086e19) ** 3          # m^3  (10 kpc cube)
EARTH_V = (6.371e6) ** 3 * 4 / 3      # m^3
GLOBAL_USE = 6.0e20       # J/yr humanity
SUN_LUM = 3.828e26        # W

MECH = [
    # name, f_bound (fraction of rho_vac), timescale_s, note
    ("N5 vacuum battery (LDPT release)", 0.3, HUBBLE_T, "f<=0.3, z_c<=0.2 from BF 0.11 (§17)"),
    ("N5 optimistic trigger", 1.0, 1e9 * YEAR, "f=1 hypothetical, tau=1 Gyr"),
    ("N4 screenable tilt (local)", 0.01, HUBBLE_T, "<0.01 mag/dex host mass (§19); tilt fraction ~1%"),
    ("RVM channel", 0.008 * 1.0, HUBBLE_T, "nu=0.008+-0.005 (§22) -> releasable ~ nu*rho_vac"),
    ("Casimir cavity (lab, not cosmo)", None, None, "u = -pi^2 hbar c / (240 a^4), a=1um gap"),
    ("Dynamical Casimir (lab)", None, None, "P ~ hbar * omega^2 * ... order: 1 W/m^2 at GHz + mm gaps (estimate)"),
]

HBAR = 1.0546e-34


def casimir(a):
    return -np_pi2() * HBAR * C / (240.0 * a ** 4)


def np_pi2():
    import math
    return math.pi ** 2


rows = []
for name, f, tau, note in MECH:
    if name.startswith("Casimir cavity"):
        a = 1e-6
        u = casimir(a)
        rows.append({"mechanism": name, "u_J_m3": u, "note": note + f" (a={a:g} m)",
                     "kind": "lab-casimir"})
        continue
    if name.startswith("Dynamical"):
        rows.append({"mechanism": name, "u_J_m3": None,
                     "note": note + "; lab-scale, not cosmological", "kind": "lab-dce"})
        continue
    if name.startswith("RVM"):
        u = f * U_VAC
    else:
        u = f * U_VAC
    e_gal = u * GAL_V
    e_earth = u * EARTH_V
    p = u / tau if tau else None
    rows.append({"mechanism": name, "f_bound": f, "u_J_m3": float(u),
                 "energy_galaxy_vol_J": float(e_gal),
                 "energy_earth_vol_J": float(e_earth),
                 "power_density_W_m3": float(p) if p else None,
                 "power_density_vs_solar": (float(p) * 1.0 / SOLAR) if p else None,
                 "years_of_humanity_per_m3": (GLOBAL_USE * YEAR / (p) if p else None),
                 "tau_s": tau, "note": note})

out = {"rho_vac_kg_m3": float(RHO_VAC), "u_vac_J_m3": float(U_VAC),
       "solar_W_m2": SOLAR, "global_use_J_yr": GLOBAL_USE, "galaxy_vol_m3": float(GAL_V),
       "rows": rows}
with open(os.path.join(RESULTS, "harness_budget.json"), "w") as f:
    json.dump(out, f, indent=2)

print(f"rho_vac = {RHO_VAC:.3e} kg/m^3, u_vac = {U_VAC:.3e} J/m^3")
for r in rows:
    if r.get("u_J_m3") is not None and r.get("f_bound") is not None:
        p = r["power_density_W_m3"]
        yrs = r["years_of_humanity_per_m3"]
        print(f"  {r['mechanism']:34s} f={r['f_bound']:.3g} u={r['u_J_m3']:.3e} J/m3 "
              f"E_gal={r['energy_galaxy_vol_J']:.3e} J  p={p:.3e} W/m3 "
              f"(1 m3 powers humanity for {yrs/3.156e7:.2e} yr)  [{r['note']}]")
    else:
        print(f"  {r['mechanism']:34s} {r.get('u_J_m3')}  [{r['note']}]")
print("BUDGET DONE")
