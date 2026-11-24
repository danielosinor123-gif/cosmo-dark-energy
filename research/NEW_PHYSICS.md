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
anchor breaks the degeneracy. MCMC fits + Bayes factors vs ΛCDM: TBD (running).

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

**Status:** implemented (`w_quintom`, `make_lnF_quintom`); fits running (TBD).

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

All numbers to be filled from `results/session6_results.json`, `h0map_s6_ldpt.json`,
`predict_s6_*.json`.
