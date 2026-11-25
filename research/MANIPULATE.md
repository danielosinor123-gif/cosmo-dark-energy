# Can we manipulate the dark-energy reading? (session 25, `results/manip_results.json`)

Two different questions hide in "steer the probe": **steer the analysis** vs **steer
the physics**. Only the first is possible.

## 1. Steering the analysis — YES, and we proved it happened (§15, §21, §27)

The "3–4σ evolving dark energy" reading is steered entirely by choices that change no
physics at all:

| choice | effect on the "detection" |
|---|---|
| H₀ anchor: Planck vs Cepheid vs TRGB vs JAGB | 4.36σ → 0.42σ / 1.65σ / 1.80σ |
| z>0.1 cut on SNe | Δχ² −6.6 → −0.0 |
| free δM step at z<0.1 | 1.72σ → 0.00σ |
| SN sample (Dovekie / DESY5 / Pantheon+ / Union3) | signal strength varies (§9, §20) |

The +0.03 mag low-z step *is* the Cepheid calibrators (§27). So the reading can be
moved by ~4σ with analysis choices — which is exactly why `dr3_plan.json`
**pre-registers** the cuts before DR3 lands. That is the honest way to "aim" the probe.

## 2. Steering the physics — NO, quantified

**(a) Fake the low-z step physically:** a +0.03 mag step at z~0.05 needs
δDL/DL = 1.4% → **δρ_Λ/ρ_Λ ≈ 2.8% uniformly inside the z<0.1 sphere**
(R ≈ 445 Mpc, V ≈ 1.1×10⁷⁶ m³, injected energy ≈ 1.6×10⁶⁵ J) — with no mechanism:
the equivalence principle says Λ cannot be concentrated or screened astrophysically,
and our own N4 test bounds environmental tilts at **<0.9% per dex of host mass**
(<0.01 mag/dex, §19) — far below what faking the step would need, and in the wrong
place (galaxies, not the 445-Mpc sphere).

**(b) Fake the w0wa hint physically:** wa=−0.5 means w≈−1.5 at z~1, i.e. the vacuum
density **grows ×2.8 from z=1 to today** — a *global* field change, visible at every
redshift, not just where we look. Our data cap |w+1| in the dip bin at 0.27±0.16
(1.7σ). A global change of Λ is a boundary condition of the whole universe
(unimodular gravity) — not a local operation; and the z>0.1 data already agree with
plain Λ (§8), so a global change is excluded, not steerable.

**(c) Casimir-shield the reading (N12 at cosmic scale):** cancelling Λ over the z<0.1
sphere with 40-μm cavities would need plate area **1.2×10²⁷× the area of the
observable universe**. Lab cavities affect only their gap volume, zero-sum.

**(d) Quantum "observation changes the vacuum":** an Unruh detector makes the vacuum
radiate, but the energy comes from the detector's agent (2nd law, §37) — the global
reading does not move.

## Bottom line

You can aim the probe by *pre-registered analysis choices* (that is legitimate, and it
is how the 2025–26 hint was manufactured — and dismantled). You cannot aim it by
*doing things to the vacuum*: Λ is uniform by the equivalence principle, our data
bound any local tilt at <1%/dex, Casimir redistribution is zero-sum and microscopic,
and a global change of Λ is a change of the universe's boundary conditions, visible
everywhere at once and excluded by the z>0.1 data. The probe reads the universe as
it is — the only steering available is choosing honestly what to look at, which is
why the DR3 plan is pre-registered.
