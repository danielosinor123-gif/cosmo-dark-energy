# E3 "Λ-flow map + harvest router" and S2-Opt implementation spec (concept level)

Plain-English rule up front: **dark energy does not flow from place to place.** In the
smooth universe its density is the same everywhere; what varies is the *distance
reading* (H, μ residuals, peculiar velocities, host mass). So the tool below **maps
where the universe reads faster/slower**, but it cannot find a “pocket of extra dark
energy,” because our own tests say there are none (sky dipole/quadrupole null §26,
host-mass tilt <0.01 mag/dex §19, vpec-robust §28). And even a real gradient could
not be piped to us — Lovelock/EP/passivity/QI plus the material-cost theorem
(`REDISTRIBUTION_THEOREM.md`) forbid concentrating or redirecting Λ.

## E3-A — conceptual tool: the "expansion-reading sky map" (Λ-flow mapper)

**Inputs**
- DESI BAO anisotropic (AP) per sky region; DR3 full-shape when public.
- SN distance residuals per sky pixel (Pantheon+ RA/DEC; DESY5/Dovekie cross-check),
  with the 77 Cepheid calibrators flagged (§27) and a free δM step at z<0.1 (§21).
- Peculiar-velocity field (e.g. CosmicFlows/2M++ reconstruction) to separate local
  matter motion from expansion (§28).
- Host stellar mass per SN to test environment tilts (§19).
- CMB (R, ω_b) and Lyα AP at z=2.33 as high-z anchors (§33).

**Outputs (sky pixels × z bins)**
- H(x,z): apparent local expansion rate.
- w(x,z) / 5-bin w(z) per region: effective dark-energy equation of state on the sky.
- δM(x), dipole/quadrupole coefficients of μ: calibration/void/systematics map.
- A "harvest score" per region: the only physically meaningful one is the enclosed
  comoving volume V(x) and its expansion power **3Hρ_Λc²V** — not a flow.

**Expected result from existing data:** no sky dipole/quadrupole in DE (Δχ²≈−0.05…−0.27
for 8 params, §26); the only strong sky signal is the **calibrator step +0.203±0.030 mag**
(§27). If DR3 changes that, the dial (`dr3_plan.json`) fires.

**Why you can't harvest toward us even if the map shows a high-H region**
- Λ is uniform by the equivalence principle; a "high-H" pixel is a distance-ladder or
  local-void effect, not extra Λ.
- To extract expansion work you must enclose a comoving volume with a boundary. Power =
  3Hρ_Λc²V. Locally that is 3.4×10⁻²⁷ W/m³; even enclosing the whole observable
  universe (V~4×10⁸⁰ m³) would give ~10⁵⁴ W **only with a piston shell at the horizon** —
  impossible and causally unavailable; per cubic meter it is still 10⁻²⁷ W.
- Redirection needs material boundaries whose rest energy exceeds the shifted Λ energy
  by R ≥ 10⁹–10¹⁴ (`REDISTRIBUTION_THEOREM.md`).

## S2-Opt — implementation spec for the optimized vacuum levitator/spring

Goal: a stable, manufacturable micro-levitator/actuator that uses squeezed-vacuum
force, with dark energy only as the 40 μm calibration mark.

**Design**
1. **Repulsive Lifshitz stack**: lower plate / fluid / upper plate, choosing the fluid
   permittivity between the plate materials so Casimir force is **repulsive** (stable).
   Candidate geometry: silica and metal films with an intermediate dielectric liquid;
   verify sign by F(d) measurement first.
2. **MEMS array**: N small cells (pitch ~100 μm) in parallel, total area A = N·a. Target
   operating gap d₀ = 0.5–1 μm for robustness (P_rep order 10⁻⁴–10⁻³ Pa). Loads:
   1 mg mirror needs ~10–100 mm² total area; a 1 kg load would need ~0.75 m² at 100 nm
   or ~7500 m² at 1 μm — so keep it micro-scale.
3. **Active gap lock**: capacitive gap sensing + piezo or fluid temperature (ε(T))
   feedback to hold d₀; belt-and-suspenders against pull-in if any attractive branch
   appears.
4. **Vacuum chamber** (≤10⁻³ mbar) to remove squeeze-film damping; optional cryo to cut
   thermal noise.
5. **Patch-charge control**: conductive coatings + Kelvin probe; null to < μV.
6. **Metrology**: white-light interferometry for d (target <1 nm), force balance for
   F(d), measure hysteresis and stable branch.
7. **Materials**: ULE glass / silicon / ALD coatings; graphene or 2D conductive layers
   for strength + charge control; high-index dielectric or metamaterial coatings to
   boost repulsive pressure per area.

**Acceptance tests**
- F(d) shows a repulsive branch with positive effective stiffness at d₀.
- Levitate a 1 mg test mass at d₀ ≥ 0.5 μm for >1 h, drift <1 nm, force margin ≥2×.
- Cycle test: work out − work in ≈ 0 within error (passivity), proving no free energy.
- Calibration: at d\* = 39.7 μm the pressure equals ρ_Λc² = 5.25×10⁻¹⁰ Pa; load there
  is only ~5×10⁻¹¹ kg/m², so d\* is a sensor point, not a levitation point.

**Honest limits**
- Not dark-energy powered: useful force lives at nm–μm gaps where QED Casimir ≫ Λ.
- 100 nm gaps can in principle carry ~10⁴ kg/m² but are fragile (roughness, stiction,
  pull-in); 0.5–1 μm is the manufacturable optimum.
- Zero-sum per cycle; material mass always dwarfs handled energy.

**Uses if built**: frictionless micro-bearings, MEMS actuators/mirrors, force standards,
position sensors, and a precise 40 μm Λ-pressure gauge (G1). Not a power source.
