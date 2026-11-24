# VERDICT.md — what 19 sessions of ruthless re-analysis say about dark energy (2026-11-24)

All numbers: our independent pipeline on public data — `results/RESULTS.md` §1–§30.
Code: `src/cosmo.py`, `src/models.py`, `experiments/session*.py`. Repo: cosmo-dark-energy.

## The 2025–26 "evolving dark energy" hint

- It is real in the data (our w0wa: Δχ²=−6.3, 1.7σ, BF 9–12; DESI's own MAP 3.05σ) — §2.
- It lives **entirely in the z<0.1 supernova tail** (z-cut: −6.6→−0.0; DESI's own
  zmin=0.1 chain: 2.09σ→1.06σ) — §8, §12.
- It is **SN-sample-dependent, not geometry-driven**: DESI's own BAO+2pt = BAO-only;
  BAO+SN+2pt = BAO+SN; Lyα alone has zero preference (Δχ²=+0.09) — §20, §29.
- It is quantitatively **equivalent to a +0.03 mag sharp step at z=0.1** (δM=+0.034±0.012
  absorbs everything; capstone Δχ²=−0.26, wa=−0.07±0.13) — §16, §21, §25.
- That step is **the Cepheid-calibrator hosts**: 77 calibrator SNe sit +0.20±0.03 mag
  (6.7σ) off BAO+CMB geometry → implied H₀=74.3±0.2 (SH0ES); the other 664 low-z SNe
  agree at −0.012±0.007 (H₀=67.7). Not a v_pec artifact (stable across z columns) — §27, §28.
- H₀↔wₐ map: wₐ ≈ −0.17 + 0.32(H₀−70); Planck-vs-SH0ES ≡ Δwₐ≈1.8 — §15.
- **Non-constant w is still Bayesianly robust** (3-bin BF 17–135; 5-bin dip
  w=−1.26±0.16 at z∈[0.75,1.5], w(0)=−0.93) — §14, §24. Allowed, not required.
- **Two separable features:** the z<0.1 step is entirely the Cepheid calibrators (§27),
  while the z~1 dip is **identical with and without them** (§31: −1.264±0.157 vs
  −1.257±0.159) — the mid-z dip is anchor-independent and is the one thing DR3 must
  confirm or kill. On the anchor-clean sample the quintom crossing (w=−1.19 at z=1,
  back to −0.93 at z=2) is still the best raw fit (Δχ²=−7.4, §32) — ~1.6σ, below the
  discovery bar but the strongest surviving candidate.

**Verdict:** Λ is not dead, but the 3–4σ headline is, in our quantification, the H₀
anchor (localized to the Cepheid calibrators) plus a +0.03 mag low-z SN step wearing a
w(z) costume. A real z~1 dip remains possible and is pre-registered for DR3.

## New physics we invented and how the data judged it

| model | verdict | key number |
|---|---|---|
| N1/N5 LDPT (late phase transition) | **killed** | BF 0.11, z_c≤0.2, f≲0.3 — §17 |
| N2 quintom phantom dip | **allowed, not required** | BF 2.2 < w0wa 8.0 — §17 |
| RVM (Λ∝ρ_m channel) | **null** | ν=0.008±0.005 — §22 |
| N4 screenable vacuum | **excluded** | <0.01 mag/dex host mass; no dust coupling (0.8–0.9σ) — §19, §26 |
| N6–N10 harnessing mechanisms | **all capped** | one-shot ≤1.6×10⁻¹⁰ J/m³; continuous ≤1.2×10⁻²⁷ W/m³; no angular/alternation signatures — §26, HARNESS_N6_N10.md |

## Can dark energy be harnessed? — the quantitative answer

- ρ_Λc² = 5.3×10⁻¹⁰ J/m³. Even a **full** release gives ~10⁻²⁶ W/m³: 1 m³ of vacuum
  powers humanity for ~10⁴⁶ yr. Dilution kills it ~26 orders of magnitude below solar
  flux — §23. Literature agrees: no net vacuum extraction without external work
  (HARNESS_LITERATURE.md; Maclay, Lafleur, Ford–Roman no-go canon).
- The only denser vacuum is the **EW false vacuum** (7.6×10⁴⁶ J/m³ = 1.5×10⁵⁷ ρ_Λ),
  bounded by Higgs stability, not by cosmology — §23.
- QI cap: sustained negative density above ρ_Λ lasts ≲0.1 ms (Ford–Roman) — HARNESS_N6_N10 §0.2.

**Verdict: dark energy is a probe, not a resource.** The lever we actually have is the
H₀ anchor decision + DR3 full-shape.

## Pre-registered DR3 falsification plan (`dr3_plan.json`, §30)

1. w0wa must reach BF>10 **with z<0.1 SNe excluded**.
2. Must survive TRGB/JAGB anchors and a free δM step.
3. 5-bin w(z): dip at z~1 with σ<0.1 and double crossing of w=−1 → N2-class physics.
4. DR3 full-shape must *move* wₐ (Y3 2pt did not, §20).
5. Precision needed to break the w0wa↔δM degeneracy at 3σ: **~3× BAO** (§30).

If 1–4 fail: the 2025–26 hint was calibration, Λ stands. If they pass: multi-field dark
energy (N2) is real, and *then* the harnessing question reopens.
