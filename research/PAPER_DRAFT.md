# PAPER DRAFT — "The 2025–26 evolving-dark-energy hint is a +0.03 mag low-z supernova step plus the H₀ anchor choice"

Working title. All numbers are from our own pipeline on public data
(`results/RESULTS.md` §1–§26; code in `src/`, `experiments/`; repo: cosmo-dark-energy).

## Abstract (draft)

DESI DR2 BAO + supernovae + CMB report a 3–4σ preference for evolving dark energy
(w₀wₐ, wₐ<0). We re-analyze the public data with an independent pipeline and show the
preference is (i) entirely in the z<0.1 supernova tail, (ii) quantitatively equivalent to
a +0.03 mag step in the low-z SN distance scale — the same offset present, mutually
consistent, in all three public SN samples — and to the local H₀ anchor
(1 km/s/Mpc ≡ Δwₐ ≈ 0.32), and (iii) absent from SN-independent geometry: DESI's own
BAO+full-shape chains are statistically identical to BAO-only, and adding 2pt to BAO+SN
does not move wₐ. With a free low-z step, a Cepheid H₀ prior, and growth all marginalized,
w₀wₐ has zero preference over ΛCDM (Δχ² = −0.26). A nonparametric 5-bin w(z) still favors
non-constant w (BF 17–135) with a localized phantom dip w = −1.26±0.16 at z∈[0.75,1.5] —
allowed, not required. We quantify the harnessing caps (one-shot ≤1.6×10⁻¹⁰ J/m³,
continuous ≤1.2×10⁻²⁷ W/m³) and pre-register the DR3 tests that would distinguish a real
z~1 transition from calibration.

## 1. Introduction

Dark energy ~68% of the universe; Λ is the simplest fit but has the 10¹²⁰ problem.
2025–26 DESI Y3: evolving-DE hint at 3–4σ in w₀wₐ (BAO+DESY5 SN). Question: is the
hint robust to systematics, or is it the known H₀ tension + low-z SN calibration wearing
a w(z) costume?

## 2. Data & pipeline

- DESI DR2 BAO "ALL" (13 rows, full cov) — data.desi.lbl.gov
- SN: DES-Dovekie HD (1820), Pantheon+ (1701, w/ host mass, RA/DEC, MWEBV), DESY5
- CMB: Planck 2018 PR3 compressed (R, ω_b) diagonal Gaussian (ℓ_A omitted — documented)
- Growth: public 24-pt fσ8 compilation, s8 free
- Models: ΛCDM, wCDM, w₀wₐ, 3/5-bin w(z), LDPT, quintom, RVM; Metropolis + Laplace BF

## 3. Results

1. Our BAO+SN+CMB w₀wₐ: Δχ²=−6.34 (1.73σ), BF 9–12 (§2). DESI's own MAP: −13.58 (3.05σ).
2. z>0.1 cut kills it: ours −0.01; DESI's own zmin=0.1 chain: 2.09σ→1.06σ (§8, §12).
3. All three SN samples show the same +0.03…+0.05 mag low-z excess vs BAO+CMB geometry,
   pairwise <1σ (§13) — one common offset ≈ ΔH₀ ≈ 1.1 km/s/Mpc.
4. H₀↔wₐ map: wₐ ≈ −0.17 + 0.32(H₀−70); Planck-vs-SH0ES ≡ Δwₐ≈1.8 (§15).
5. δM step: +0.022±0.020 mag at z<0.1 → 1.72σ→0.00σ (§16). Edge test: sharpest step
   preferred, Δχ²=−7.3 (§25).
6. DESI's own chains: BAO+2pt = BAO-only; BAO+SN+2pt = BAO+SN (wₐ unmoved, §20) —
   the signal is SN-sample-dependent, not geometry-driven.
7. Nonparametric: 3-bin BF 17–135 (§14); 5-bin: phantom dip w=−1.26±0.16 at z∈[0.75,1.5],
   w(0)=−0.93, growth Λ-like (§24).
8. Capstone (all systematics free): w₀wₐ vs ΛCDM Δχ²=−0.26, wₐ=−0.07±0.13, δM=+0.034±0.012
   (§21).
9. New physics judged by data: LDPT BF 0.11 (killed); quintom BF 2.2 (allowed, not
   required); RVM ν=0.008±0.005 (null) (§17, §22).
9b. The H₀ tension is localized to the 77 Cepheid-calibrator SNe: they sit
   +0.20±0.03 mag (6.7σ) above BAO+CMB geometry — implying H₀=74.3±0.2 — while the 664
   other z<0.1 SNe agree at −0.012±0.007 mag (H₀=67.7). The full-sample "+0.03 mag
   low-z offset" is the diluted signature of these calibrators (§27, v_pec-robust §28).
9c. Two separable features: the z<0.1 step is entirely the calibrators, but the z~1
   phantom dip is identical without them (−1.264±0.157 vs −1.257±0.159, §31) —
   anchor-independent; the Y3 Lyα BAO shows no preference (§29); the 2026 DR2 Lyα
   full-shape AP (arXiv:2607.27410, 1% at z=2.33) is ΛCDM-consistent alone (χ²=0.97/2,
   Lyα+BBN H₀=66.5±1.3) and weakens w0wa when added (−6.6→−5.6, §33) — the high-z
   anchor moved toward ΛCDM; DR3 needs ~3× BAO precision to break the w0wa↔δM
   degeneracy at 3σ (§30).
10. Harnessing caps: one-shot ≤1.6×10⁻¹⁰ J/m³, continuous ≤H₀ρ_Λc²≈1.2×10⁻²⁷ W/m³,
    tilt <0.01 mag/dex; EW vacuum 57 orders denser than Λ (§23, §26,
    HARNESS_N6_N10.md, HARNESS_LITERATURE.md).

## 4. Discussion

Λ survives as most parsimonious. The hint = low-z SN step + H₀ anchor. Growth vs
geometry mild tension (constant w≈−0.9 from fσ8).

## 5. Falsification / pre-registered DR3 tests

(a) w₀wₐ with z>0.1 SNe only: discovery requires BF>10 there. (b) Anchor swap to
TRGB/JAGB: signal must survive. (c) 5-bin w(z): dip at z~1 with σ<0.1 and double
crossing of w=−1 → N2-class physics; collapse → calibration. (d) DR3 full-shape must
move wₐ (Y3 2pt did not). (e) Non-Cepheid low-z anchor (e.g. SNe-only/BAO-anchored
calibration) removing the z<0.1 step.

## Caveats

Compressed CMB prior (R,ω_b only; ℓ_A off ~2% in our r_s); fσ8 24-pt public
compilation; BAO per-tracer covs diagonal where 1×1; single low-z offset could still be
a real z≲0.1 transition (degenerate with calibration — the sharp edge, §25, keeps that
possibility open until (e)).
