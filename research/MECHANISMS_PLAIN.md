# What each conceptual mechanism actually does (plain language)

One-line version: **G1 reads Λ, S2 borrows boundary force, E3/P5 tap expansion or phantom growth at ~10⁻²⁷ W/m³, N20 is a dead light-bender.**

| mechanism | what it physically does | input → output | honest limit |
|---|---|---|---|
| **G1 Λ-pressure gauge** | Holds two plates at the matching gap d\* = 39.7 μm and measures the tiny Casimir pressure there, which equals ρ_Λc² = 5.25×10⁻¹⁰ Pa | laser gap control + force sensor → a ppm readout of the vacuum density ρ_Λ (1 m² + 1 fN ⇒ ~2×10⁻⁶) | measures normalization only; yearly drift ~2×10⁻¹¹ is invisible; not a power source |
| **S2 vacuum spring / Casimir levitator** | Uses Casimir attraction (or repulsive Lifshitz in a fluid) between plates as a force: at nm gaps it can hold kg on cm²; at d\* it equals the cosmic pressure but is only ~5×10⁻¹¹ kg/m² | plates + fluid/actuator → levitation/actuation force P·A | attractive case is unstable; reset is zero-sum; real load uses nm gaps where Casimir ≫ Λ, so dark energy is only the 40 μm calibration point |
| **E3 expansion harvester / tether** | Lets cosmic expansion do work on a comoving piston, or holds a tether at fixed proper separation | comoving volume → 3Hρ_Λc² = 3.44×10⁻²⁷ W/m³; tether force on 1 kg ≤ 10⁻¹³ N even at 1 Mpc | 1 W needs 2.9×10²⁶ m³ ≈ 13 million Moon volumes; no ratchet beats this (2nd law, T_dS=2.7×10⁻³⁰ K) |
| **P5 phantom-dip rudder** | If DR3 shows the z~1 dip at ≥5σ, couples a device to the phantom dark field to locally modulate Λ_eff | portal coupling gφF²/ξRφ² → local Λ_eff dial; tap 3H|1+w|ρ_Λc² = 8×10⁻²⁸ W/m³ | gated on DR3 ≥5σ + surviving z>0.1/anchors/δM/2pt; even then 1 W needs 1.2×10²⁷ m³ and host-mass tilt is bounded <0.01 mag/dex |
| **N20 vacuum prism** | A slab that cancels Λ would act as a negative-index lens for light | would deflect light by θ ~ (8πG/c⁴)(ρ_Λc²)L | θ ~ 10⁻⁵² rad for L=1 m — dead by 52 orders; material cost R≥10⁹ |

**Bottom line:** none of them moves dark energy across the sky or powers a device.
G1 is a sensor of Λ; S2 is an ordinary Casimir actuator whose 40 μm point is Λ-calibrated;
E3/P5 can only harvest expansion/phantom growth at ~10⁻²⁷ W/m³; N20 does nothing.
The only thing that could change this picture is a ≥5σ DR3 dip (`dr3_plan.json`).
