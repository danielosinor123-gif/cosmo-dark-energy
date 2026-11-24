# New-physics proposals for dark energy (session 6) — motivated by our own measurements

Everything below is constrained by what *our* pipeline measured on public data
(RESULTS.md §1–§15):

- binned w(z): w(0) ≈ −0.93, dip to ≈ −1.2 at z ∈ [0.75, 1.5], back to ≈ −0.9 at z > 1.5
  (probe-consistent, max tension 0.51σ; nonparametric BF vs constant w = 17–135)
- w0wa: w0 ≈ −0.87, wa ≈ −0.35, but carried by z < 0.1 SNe and degenerate with H₀
  (wa ≈ −0.17 + 0.32·(H₀ − 70)); Cepheid anchor → 0.4σ, Planck anchor → 4.4σ
- growth (24-pt fσ8): consistent with *constant* w ≈ −0.9; w0wa adds nothing once fσ8
  is in the joint (Δχ² ≈ 0 vs wCDM)
- one coherent low-z SN offset +0.03…+0.05 mag common to Dovekie/Pantheon+/DESY5

Design targets for a viable model: (i) w(0) slightly above −1, (ii) a phantom-ish dip
near z ~ 1, (iii) growth close to ΛCDM (no big fσ8 distortion), (iv) optionally absorbs
part of the H₀ tension.

---

## N1 — LDPT: Late Dark-sector Phase Transition ("late vacuum metamorphosis")

**Mechanism.** The dark sector undergoes a first-order phase transition at z = z_c
(~1–2): a fraction f of the vacuum energy ρ_DE is released and converted into dark
matter (or radiation) — a *late-time* analogue of vacuum metamorphosis (Wetterich 1988)
and of early-dark-energy, but at z ~ 1 instead of z ~ 3000. Parametrisation:

ρ_DE(z) = ρ_DE0 · [1 − f·σ((z−z_c)/Δ)] / [1 − f·σ(−z_c/Δ)],  σ = logistic,

so ρ_DE(0) is unchanged, the density steps down by f across z_c, and
w_eff(z) = −1 − (1/3) d ln ρ_DE / d ln a develops a dip below −1 inside the transition.

**Why it fits our constraints.** A step at z_c ~ 1 naturally produces the binned dip at
z ∈ [0.75, 1.5]; the released energy adds to matter at z > z_c, which *raises* the
high-z matter density relative to ΛCDM and shifts the BAO/CMB-inferred H₀ — potentially
absorbing part of the H₀ tension (our §15 map shows ΛCDM at z>0.1 wants H₀ ≈ 70.5,
SH0ES says 73).

**Predictions (falsifiable).**
1. w_eff(z) has a *localised* dip below −1 confined to z_c ± 2Δ — DR3/full-shape or a
   5-bin w(z) reconstruction should show a step-like feature, not a smooth CPL ramp.
2. Ω_m(z > z_c) is enhanced by ≈ f·Ω_DE relative to ΛCDM → fσ8 at z > 1.5 is shifted
   (we test against the public 24-pt fσ8 benchmark in the joint fit).
3. χ²(H₀) flattens and its minimum moves toward H₀ ≈ 73 (h0map_s6_ldpt).
4. CMB damping tail / lensing: a late transition changes the late-time ISW and the
   matter density at recombination is unchanged, so Planck alone barely sees it — the
   signature is in late-time distances + growth, exactly where our data look.

**Status:** implemented in `src/models.py` (`lnF_ldpt`, `BackgroundFn`); fit in
`experiments/session6.py`.

χ²(H₀) profile on BAO+Dovekie(full)+CMB-compressed (`h0map_s6_ldpt.json`):

| H₀ | χ² ΛCDM | χ² LDPT | Δχ² | LDPT best (f, z_c, Δ) |
|---|---|---|---|---|
| 66 | 2136.7 | 2133.3 | −3 | (0.73, 2.2, 0.03) |
| 70 | 1647.4 | 1646.6 | −1 | (0.09, 0.9, 0.03) |
| 72 | 1726.9 | 1690.4 | **−36** | (0.40, 0.2, 0.56) |
| 73 | 1847.6 | 1742.1 | **−105** | (0.46, 0.2, 0.40) |
| 74 | 2021.5 | 1796.8 | **−225** | (0.31, 0.2, 0.12) |

**Reading:** LDPT eats a large and growing chunk of the high-H₀ misfit — at H₀=73 it
removes Δχ² ≈ 105, at H₀=74 ≈ 225 — so a late dark-sector transition *can* absorb most
of the SH0ES-vs-BAO/CMB H₀ gap. Caveat (important): at high H₀ the fit pushes z_c to the
prior boundary z_c ≈ 0.2, i.e. the "transition" the data want sits at z ≲ 0.2 — exactly
where the low-z SN calibration step lives (§8, §13, §16). Physics-at-z~0.2 and a low-z
calibration systematic are degenerate in current data; DR3/full-shape or a non-Cepheid
anchor breaks the degeneracy. **MCMC verdict (RESULTS.md §17): LDPT does not improve the
full-sample fit at all (Δχ² = −0.34 no growth, −0.91 with fσ8) and is disfavoured by the
Laplace evidence (BF = 0.11 with growth) — killed by Occam.**

---

## N2 — QUINTOM-DIP: localized phantom-divide crossing (two-field dark sector)

**Mechanism.** A two-field dark sector (quintom: one quintessence + one phantom
component, or a single field with a non-canonical kinetic term that permits a
transient crossing) whose *effective* equation of state crosses w = −1 near z ~ 1:

w(z) = w_b + A·exp(−((z−z_c)/Δ)²),   A < 0,

with w(0) = w_b + A·e^{−(z_c/Δ)²} > −1, a phantom dip w < −1 around z_c, and w → w_b
at high z. Single-field canonical scalar models *cannot* cross w = −1 (would require a
ghost), so observing a crossing is a smoking gun for multi-field/non-canonical dark
energy — this is the physics claim.

**Why it fits our constraints.** It reproduces the measured binned shape directly
(shallow above −1 at z=0, dip at z~1, back above −1 at high z) with 4 parameters, and
because the phantom phase is confined to z ~ 0.6–1.6, the integrated effect on growth
(fσ8 at z ≲ 2.5 mostly weights z < 1.5) can stay close to ΛCDM — consistent with §7.

**Predictions (falsifiable).**
1. w(z) crosses −1 *twice* (at z ≈ z_c ± Δ·arcosh(√…)) — a 5-bin reconstruction or
   w0wa-vs-quintom BF (session6) discriminates: if the data prefer the dip over CPL,
   crossing is supported.
2. Growth: phantom DE with no DE perturbations gives *enhanced* growth inside the
   phantom window (f slightly above ΛCDM at z ~ 1) — testable with the fσ8 benchmark
   and with RSD at z ~ 1–2 (DESI DR3).
3. No H₀ absorption (pure background shape) — the χ²(H₀) profile stays ΛCDM-like;
   if the real driver is H₀ systematics, N2 alone won't fix it (contrast with N1).

**Status:** implemented (`w_quintom`, `make_lnF_quintom`); fitted (`session6_results.json`,
RESULTS.md §17). Best fits: no-growth w_b=−0.92, A=−0.68, z_q≈1.0, Δ=0.14 (narrow deep
phantom spike, w(1)≈−1.6 at MAP); with fσ8 z_q=1.93±0.51, A=−0.52±0.24 (broader dip).
Δχ² vs ΛCDM: −8.06 (no growth, k=4), −7.70 (with growth); Laplace BF = 2.2 with growth —
**consistent with the data but not required**: w0wa (BF 8–9, k=2) is the more
parsimonious description. A confirmed double crossing would still be a discovery, but
current public data do not demand it.

---

## N3 — (screened) IDE-threshold: interacting DE switching on at late times

Q = 3βHρ_DE for z < z_on, else 0 (dark energy decays into dark matter only after z_on).
Changes both E(z) and growth; needs the perturbed continuity equations (not the Heath
integral), so deferred — flagged here because it is the natural model if fσ8 at low z
turns out to deviate from ΛCDM (our current data say it doesn't).

---

## What would count as a discovery

- BF(N1 or N2 / ΛCDM) > 10 on BAO+SN+CMB **and** survival under the §8 z>0.1 cut,
  §12 anchor swap, and §7 growth joint — currently the bar none of w0wa/wCDM/N-model
  candidates has cleared robustly.
- A confirmed double crossing of w = −1 (N2) at >3σ with 5+ bin reconstruction.
- A measured χ²(H₀) flattening of N1 with minimum at H₀ ≈ 73 (N1 absorbing the H₀
  tension) — that would be the first late-time dark-sector transition candidate
  consistent with both SH0ES and BAO/CMB.

## Final verdict (session 6, RESULTS.md §16–§17)

- No new-physics model beats plain w0wa on evidence; LDPT is rejected (BF 0.11),
  quintom-dip is allowed but not required (BF 2.2 vs w0wa's 8.0).
- The single most quotable number (§16): **a δM ≈ +0.02–0.03 mag step at z<0.1 absorbs
  the entire evolving-DE preference (1.72σ → 0.00σ)** — same size as the coherent
  low-z offset common to all three SN samples (§13).
- Discovery bar (BF>10 AND survival of z>0.1 cut + anchor swap + δM step): met by
  nothing. Predictions for DR3/non-Cepheid anchors stand as written above.

## N4 — "Screenable vacuum energy" (the harnessable-vacuum hypothesis)

**Mechanism (novel extrapolation).** Treat ρ_vac not as a constant but as the energy of an
order parameter (a condensate / false-vacuum field) that local matter density can tilt —
unscreened in voids, shifted in massive halos. If so, the vacuum energy is *locally
modifiable in principle* — the only route to "harnessing" it that doesn't violate energy
conditions — and H₀/μ measured in different environments should differ.

**Prediction.** Low-z SN distance residuals correlate with host-galaxy stellar mass.

**Test + verdict (RESULTS.md §19, `env_results.json`).** Pantheon+ host masses: host-split
Δμ = +0.0008 ± 0.0096 mag (0.08σ) at z<0.3; slope dr/dlogM = −0.0055 ± 0.0035 mag/dex
(1.57σ, sign opposite to the z>0.3 control). **No environmental coupling detected.**
Bound: matter–vacuum coupling < ~0.01 mag/dex in host mass at z<0.3. So: either the vacuum
is not locally tiltable at galaxy scales, or the coupling is below this bound. The common
low-z offset is global, not environmental.

## N5 — "Metastable vacuum battery" (decay-triggered release)

**Mechanism.** The dark sector sits in a false vacuum; a first-order transition (natural
at z_c, or nucleated) converts a fraction f of ρ_vac into radiation/DM — a global "battery
discharge". Cosmological signature = the LDPT step (N1).

**Energy budget if harnessable.** ρ_vac c² ≈ 6×10⁻¹⁰ J/m³. With the data bound f ≲ 0.3
(§17: LDPT disfavored, z_c ≲ 0.2), the releasable density is ≲ 2×10⁻¹⁰ J/m³; a
galaxy-scale volume (~10⁶³ m³) holds ≲ 2×10⁵³ J ≈ 0.1 M☉c² — modest per volume, but a
*global* transition would also change the expansion history, which our data already
constrain (BF 0.11). Verdict: even as a "battery," current data cap f and pin z_c to the
low-z SN step — indistinguishable from calibration (§16).

## Session 7–8 harness summary

- H₀-tension hypotheses tested quantitatively (RESULTS.md §18): early-universe fixes
  (ΔN_eff, EDE) absorb **nothing** in public compressed data; only late-time freedom does
  (w0wa flat in H₀; LDPT Δχ²≈105 at H₀=73, degenerate with the δM step).
- "Harnessable vacuum" (local tilt, N4) excluded at <0.01 mag/dex host-mass level (§19);
  "vacuum battery" (N5) capped by the same f ≲ 0.3, z_c ≲ 0.2 bounds as N1.

All numbers in `results/session6_results.json`, `results/step_results.json`,
`results/h0tension_results.json`, `results/env_results.json`,
`results/h0map_s6_ldpt.json`, `results/predict_s6_*.json`.

---

## N6–N10 — harnessing proposals (S10 subagent)

Full write-up in **`research/HARNESS_N6_N10.md`** (mechanism Lagrangians, no-go evasion,
energy budgets, public-data tests, 5–15 yr experiments, comparison table, references).
One-liners:

- **N6** vacuum-condensate phonon tap (ghost-condensate w=−1 with real phonon + F²φ/4
  channel): continuous cap H0·ρ_Λc² ≈ 1.2×10⁻²⁷ W/m³; decisive test = fσ8 in the
  z∈[0.75,1.5] dip bin (24-pt, Δχ²≈0).
- **N7** local bubble-nucleation catalysis of the dark false vacuum: one-shot cap
  f·ρ_Λc² ≲ 1.6×10⁻¹⁰ J/m³ (LDPT BF 0.11, z_c≤0.2); dead as power (0.17 W/m² of wall);
  decisive test = angular spectrum of low-z SN residuals vs bubble template (§13: none).
- **N8** Casimir negative-energy cavity "vacuum pump": QI caps sustained ρ_neg>ρ_Λ to
  τ≲0.1 ms (drive must pay); Λ shows up as a constant −5.3×10⁻¹⁰ Pa offset at gaps
  ≳40 μm — measurable; decisive test = low-z SN residual vs host dust A_V (cap
  <0.01 mag/dex, §19).
- **N9** axion-condensate stirring: resonant power ~10⁻⁵² W/m³ (dead at m~H0); explains
  the z~1 dip shape (m~H(z~1)); decisive test = cert5 5-bin w(z)+fσ8 vs cosine template.
- **N10** time-crystal vacuum: T_TC≳3 Gyr from smooth 5-bin w(z); A≲0.3; decisive test =
  cert5 alternating-bin residuals + clock/21-cm modulation searches.

Common caps: one-shot ≤1.6×10⁻¹⁰ J/m³ (degenerate with δM step at z<0.1), continuous
≤1.2×10⁻²⁷ W/m³, static local tilt <0.01 mag/dex (§19). None solves the 10^120 problem;
all are falsifiable with existing public data + one lab measurement each.
