# Concept designs for the mechanisms — crazy but internally consistent, each with a computed kill switch
(session 31, `experiments/session31.py concepts` → `results/concept_designs.json`)

Everything below is a **design envelope**, not a claim of a working device. The rule we
followed: invent freely, but every concept must carry its own number that kills it if
physics says no.

## G1 — the Λ-pressure gauge (real instrument concept)

**Design.** Twin parallel-plate Casimir cavities, area A, gap servoed around
**d\* = 39.7 μm** where P(d\*) = K/d\*⁴ = ρ_Λc² = **5.253×10⁻¹⁰ Pa**. Differential force
vs a reference gap, Kelvin-probe patch compensation, interferometric gap metrology.
Output: d\* = (K/ρ_Λc²)¼ is a direct readout of the cosmological vacuum density **if
the QED Casimir coefficient K is trusted**.

**Sensitivity (`concept_designs.json`).** σ_P = F_noise/A, so σ_ρ/ρ = σ_P/ρ_Λc²:

| area | force noise | σ_P | σ_ρ/ρ | gap resolution | patch budget |
|---|---|---|---|---|---|
| 1 mm² | 1 fN | 1e-9 Pa | 190% | 19 μm | 0.6 mV |
| 1 cm² | 1 fN | 1e-11 Pa | **1.9%** | 190 nm | 60 μV |
| **1 m²** | **1 fN** | **1e-15 Pa** | **1.9×10⁻⁶ (ppm)** | **19 nm** | **0.6 μV** |
| 1 km² | 1 nN | 1e-15 Pa | 1.9×10⁻⁶ | 19 nm | 0.6 μV |

So a 1 m² cavity with fN force noise and sub-μV patch potentials measures the
normalization of ρ_Λ to **~2 ppm**. **Kill switch:** cosmological drift is invisible on
lab timescales — δρ/ρ per year = 3(1+w)H₀Δt ≈ **2×10⁻¹¹** (w=−0.9) to −6×10⁻¹¹
(w=−1.27) — far below ppm, so this is a one-shot gauge, not a drift monitor.
**Falsification:** if d\* disagrees with the cosmology-inferred value by >ppm after
systematics, “Casimir QED normalization = cosmological Λ” fails → new physics.

## S2 — vacuum spring / Casimir levitator (pull/levitation)

Load capacity per area for Casimir pressure P=K/d⁴ is m/A = P/g:

| gap | pressure | load | area to hold 1 kg |
|---|---|---|---|
| 10 nm | 1.3×10⁵ Pa | **1.33×10⁴ kg/m²** | **0.75 cm²** |
| 100 nm | 13 Pa | 1.33 kg/m² | 0.75 m² |
| 1 μm | 1.3×10⁻³ Pa | 1.3×10⁻⁴ kg/m² | 7.5×10³ m² |
| **39.7 μm = d\*** | **5.25×10⁻¹⁰ Pa = ρ_Λ** | **5.4×10⁻¹¹ kg/m²** | **1.9×10¹⁰ m² (≈19,000 km²)** |
| 100 μm | 1.3×10⁻¹¹ Pa | 1.3×10⁻¹² kg/m² | 7.5×10¹¹ m² |

Attractive Casimir is **unstable** (force falls as the gap opens); a real levitator
needs repulsive Lifshitz (fluid with ε₁<ε_fluid<ε₂) or a metamaterial/active feedback.
**Dark-energy connection:** only at d\* does the pressure equal ρ_Λ — there the load is
~10⁻¹⁰ kg/m², useless for levitation. Real levitation uses nm gaps where Casimir is
≫ Λ, so the force is QED boundary physics, not dark energy. Kills: roughness, patch
potentials, squeeze-film damping, stability.

## E3 — expansion tether / harvester (the only net-positive term)

Holding a mass m at fixed proper separation r in de Sitter needs
F = m·(8πGρ_Λ/3c²)·r:

| r | F on 1 kg |
|---|---|
| 1 m | 3.3×10⁻³⁶ N |
| 1 km | 3.3×10⁻³³ N |
| 1 AU | 4.9×10⁻²⁵ N |
| 1 Mpc | 1.0×10⁻¹³ N |

A “comoving ratchet” cannot beat the expansion cap 3Hρ_Λc² = **3.44×10⁻²⁷ W/m³**: 1 W
needs **2.9×10²⁶ m³** (cube side 6.6×10⁸ m ≈ **1.3×10⁷ Moon volumes**). The de
Sitter horizon heat engine is colder still: T_dS = 2.7×10⁻³⁰ K → total 6.7×10⁻⁷³ W.

## P5 — phantom-dip rudder (only if DR3 pushes the dip past 5σ)

If the z∈[0.75,1.5] dip reaches **|w+1|/σ ≥ 5** and survives z>0.1, anchors, δM, 2pt,
then w<−1 there and a portal-coupled device (gφF², ξRφ²) could locally modulate
Λ_eff. Posterior power 3H|1+w|ρ_Λc² = **8.08×10⁻²⁸ W/m³**; 1 W needs **1.24×10²⁷ m³**
(side 1.07×10⁹ m) — *worse* than plain expansion today. So even if real, it is a
**probe/signature**, not a power source; its falsifiers are the pre-registered DR3 cuts
plus the host-mass tilt bound <0.01 mag/dex (§19).

## N20 — the crazy one we killed: “vacuum prism / Λ sky lens”

Idea: a slab that locally cancels Λ should bend light like a lens with index from
δρ = −ρ_Λ. Deflection over thickness L: θ ~ (8πG/c⁴)(ρ_Λc²)L. For L = 1 m,
θ ~ **10⁻⁵² rad** — 52 orders below anything. Material cost from
`REDISTRIBUTION_THEOREM.md` (R ≥ 10⁹–10¹⁴) kills it anyway. Dead on arrival; kept here
so the crazy ideas are on record.

## Bottom line

The designs reduce to: **G1 is a real ppm-level Λ-normalization gauge at 40 μm; S2 can
levitate only at nm gaps where Casimir ≫ Λ (the d\* point holds 10⁻¹⁰ kg/m²); E3 is
capped at 3.4×10⁻²⁷ W/m³; P5 is gated on a ≥5σ DR3 dip and is still ≤10⁻²⁷ W/m³.**
Dark energy enters only as the calibration point of the gauge and the possible phantom
signature — never as a force or fuel.
