# Why Λ cannot be steered physically — the full argument (session 25b)

"Manipulating the reading" = changing the local effective cosmological term
Λ_eff = Λ + δ(x). Here is every knob that has ever been proposed, and exactly where
each one dies. Numbers from `results/manip_results.json`, `harness_caps.json`, §19/§28/§33.

## 0. What a knob would have to be, mathematically

Einstein's equation: G_μν + Λ g_μν = 8πG T_μν. A "local Λ" means a contribution to
T_μν proportional to g_μν, i.e. a fluid with p = −ρ exactly. Three theorems close
the door:

1. **Lovelock (1917):** in 4D the only divergence-free, metric-only, second-order
   tensor is G_μν + Λg_μν with Λ a *constant of nature*. There is no "material" whose
   stress tensor is an adjustable multiple of g_μν. To make Λ_eff vary in space you
   need extra *fields*, and every field obeys:
2. **NEC / passivity:** for a scalar, ρ + p = φ̇² ≥ 0, so any dynamical stand-in has
   w ≥ −1 — it can only *approach* Λ from above, never be pumped past it, and the
   ground state (the actual vacuum) is **passive**: no cyclic process extracts work
   from it (2nd law). Turning w past −1 (phantom) needs ρ + p < 0 — allowed only if
   the z~1 dip is real, and even then |w+1| = 0.27±0.16 in one bin (§34).
3. **QI (Ford–Roman):** sustained negative energy density in free fields is bounded
   (|ρ_neg| τ² ≲ ħ/c³); the only sustained source is material boundaries (Casimir),
   and those are zero-sum — the plates hold the energy.

## 1. The knobs, one by one

| knob | what it would need | verdict (number) |
|---|---|---|
| **Concentrate Λ astrophysically** (screening, chameleon, N4) | δρ/ρ ~ 2.8% inside z<0.1 sphere to fake the step | excluded: host-mass tilt <0.9%/dex (<0.01 mag/dex, §19, <1.6σ) |
| **Backreaction / local void** | δH/H ~ 2.5% near us → δμ ≈ 0.054 mag at z~0 | REAL but only at z≈0 — it can steer the *anchor* (this is the local-void H₀ effect, and why §28 vpec test mattered); at z~1 peculiar velocities give δμ ~ 0.002 mag — cannot touch the dip |
| **Lensing** | stochastic σ_μ ≈ 0.02 mag at z~1 | dip signal ≈ 0.10 mag coherent across 3 SN samples → 5× above lensing noise; lensing cannot *fake* a coherent dip |
| **Casimir redistribution at scale (N12)** | plates over the z<0.1 sphere | area needed = 1.2×10²⁷ × observable-universe area; lab gaps are 40 μm, zero-sum |
| **Phantom field local excitation (if dip real, N13)** | couple to ∂_tρ_Λ | power ≤ 8×10⁻²⁸ W/m³ (posterior, §37); local tilt would show in host-mass/dex — excluded at <1%/dex |
| **Bubble nucleation (N7/N11)** | catalyze vacuum decay | EW rate ~10⁻⁸⁵⁸ GeV⁴; cosmological LDPT BF 0.11, no sky dipole (<1σ, §26) |
| **Global change of Λ** | unimodular: Λ is the boundary condition of the whole spacetime | not a local operation; shifts *all* redshifts — excluded by z>0.1 agreeing with Λ (§8) and Lyα AP at z=2.33 (χ²=0.97/2, §33) |
| **Quantum observation (Unruh)** | detector sees vacuum as a bath | energy billed to the accelerating agent (2nd law); global reading unmoved (§37) |
| **Horizon engines** | harvest de Sitter horizon | T_dS = 2.7×10⁻³⁰ K → total 6.7×10⁻⁷³ W (§37) |

## 2. Why the *reading* itself is self-protecting

Λ enters **every** distance and **every** redshift, so any manipulation must be
coherent across ~10⁷⁶ m³ — and every independent tracer must agree, or the
manipulation shows up as a *disagreement*. Our multi-tracer checks are exactly this
vigilance:

- BAO ratios (DM/rd, DH/rd) are insensitive to local masses — they read only the
  background; they agree with Λ at z>1.6 and with the low-z SN story (§29/§33).
- Three SN samples (Dovekie / DESY5 / Pantheon+) give the *same* dip to 0.016 (§34)
  — a local manipulation would have to be sample-independent, i.e. cosmological.
- DESI's own published chains reproduce our w(0.25)>−1, w(1)≈−1.1…−1.4 (§37).
- CMB (R, ω_b) and Lyα AP pin the high-z anchor independently of the SN systematics.

There is no slack: the readings agree with each other, so there is nothing to steer.

## 3. The one theoretical crack — and its size

If the z~1 dip survives DR3, then w < −1 there: the vacuum is a *dynamical phantom
field* (N2), NEC is violated in that epoch, and a local coupling (∂_μρ_Λ) becomes
conceivable in principle. But even then: the field is homogeneous on large scales,
any local gradient costs kinetic energy and shows up as a tilt (bounded <1%/dex),
and the extractable power is the growth rate itself, 3H|1+w|ρ_Λc² ≤ 8×10⁻²⁸ W/m³
(§37) — 26 orders below sunlight. So the crack exists only as a *measurement*, not a
knob: DR3 decides whether the vacuum is a moving field, and if it is, we read it —
we still can't turn it.

## Bottom line

Λ is not a substance in space; it is a property of spacetime (or of the fields'
ground state). Locally you can manipulate *curvature* — that's just mass, and it
steers only the z≈0 anchor (local void, ~0.05 mag) — but the constant part cannot be
varied locally without violating NEC/Lovelock/passivity, and globally it is a
boundary condition of the universe, visible at all z and pinned by the data. The
probe is self-protecting precisely because every tracer reads the same Λ. What
remains steerable is *what we choose to measure and what we pre-declare as a
detection* — which is why the DR3 plan is pre-registered (`dr3_plan.json`).
