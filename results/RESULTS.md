# Results — independent re-analysis of the DESI DR2 evolving-dark-energy hint (Nov 23, 2026)

Pipeline: `src/cosmo.py` — our own numpy/scipy code (flat ΛCDM / wCDM / w0waCDM, binned-w option,
numerical sound horizon, adaptive Metropolis). No cosmology packages used.
Data (all public):
- DESI DR2 BAO, all 7 tracers (CobayaSampler/bao_data `desi_bao_dr2`, 13-dim combined vector).
- DES-Dovekie SN Ia Hubble diagram (CobayaSampler/sn_data, 1820 SNe, precomputed inverse
  total covariance).
- Planck 2018 TT,TE,EE+lowE compressed prior as a *diagonal* Gaussian on (R, ω_b).
  ℓ_A omitted on purpose: our EH98-z_d numerical r_s is ~2–3% off CAMB, which would be a
  ~100σ systematic in ℓ_A; documented here instead. DESI's own "CMB-compressed" is the PR4
  (θ*, ω_b, ω_c) triple, whose data file is NERSC-only — so we use PR3 (R, ω_b).
Validation reference: DESI's own published DR2 chains (table at the bottom).

## 1. Main model comparison (`fit_models.py` → `main_results.json`)

w0wa vs ΛCDM, Δχ²_MAP and equivalent Gaussian σ (2 extra params):

| Data combo | w0wa best-fit (w0, wa) | Δχ² (w0wa−ΛCDM) | σ | wCDM w |
|---|---|---|---|---|
| **BAO + SN + CMB(R,ωb)** | −0.872 ± 0.066, −0.345 ± 0.303 | −6.34 | **1.73σ** | −0.957 ± 0.022 |
| BAO + SN | −0.857 ± 0.068, −0.439 ± 0.485 | −6.31 | 1.72σ | −0.914 ± 0.036 |
| BAO + CMB | −0.643 ± 0.200, −0.938 ± 0.627 | −2.90 | 0.72σ | −0.992 ± 0.039 |
| BAO only | −0.579 ± 0.240, −1.300 ± 0.955 | −4.59 | 1.28σ | −0.927 ± 0.080 |
| SN + CMB | −0.855 ± 0.110, −0.457 ± 0.558 | −4.12 | 1.14σ | −0.953 ± 0.024 |

Best-fit ΛCDM (BAO+SN+CMB): Ωm = 0.303 ± 0.006, h = 0.703 ± 0.002, r_d = 143.9 ± 0.6 Mpc,
χ² = 1645.6 (≈1835 dof — good fit).

**Reading:** every combination prefers w0wa over ΛCDM at 0.7–1.7σ, always with w0 > −1 and
wa < 0 — the same *direction* as DESI's 3.1σ hint, but weaker in magnitude. Our CMB prior
(R, ω_b only) is much weaker than DESI's full compressed (θ*, ω_b, ω_c), which is the main
reason our significance is ~1.7σ rather than ~3.1σ.

## 2. Per-tracer leave-one-out (`tracer_results.json`)

Each DESI tracer alone + SN + CMB, w0wa vs ΛCDM:

| tracer (z_eff) | w0 | wa | Δχ² | σ |
|---|---|---|---|---|
| BGS (0.295) | −0.856 ± 0.122 | −0.481 ± 0.599 | −4.03 | 1.11σ |
| LRG z0 (0.510) | −0.891 ± 0.112 | −0.280 ± 0.538 | −3.98 | 1.10σ |
| LRG z1 (0.706) | −0.896 ± 0.111 | −0.291 ± 0.499 | −3.97 | 1.09σ |
| LRG+ELG (0.934) | −0.843 ± 0.114 | −0.540 ± 0.566 | −5.43 | 1.51σ |
| ELG (1.321) | −0.851 ± 0.112 | −0.513 ± 0.533 | −4.65 | 1.29σ |
| QSO (1.484) | −0.881 ± 0.114 | −0.369 ± 0.548 | −3.88 | 1.06σ |
| Lyα (2.330) | −0.877 ± 0.117 | −0.379 ± 0.564 | −3.81 | 1.04σ |

**Reading:** the w0wa preference is present at ~1–1.5σ in *every* tracer independently —
no single tracer (not even Lyα, whose full-shape analysis famously pulls back toward ΛCDM
in DESI's own work) drives or kills the signal. w0 ≈ −0.87 ± 0.11 consistently.

## 3. Probe-consistency certificates, binned w(z) (`probe_consistency.py` → `probe_results.json`)

w(z) piecewise-constant in bins [0,0.75), [0.75,1.5), [1.5,∞); each row = one probe subset
(+ SN + CMB unless noted). Tension T = (w_probe − w_joint)/√(σ²_probe + σ²_joint) per bin.

| probe subset | w(z<0.75) | w(0.75–1.5) | w(z>1.5) | max \|T\| vs joint |
|---|---|---|---|---|
| **joint (BAO+SN+CMB)** | −0.93 ± 0.03 | **−1.31 ± 0.21** | −0.97 ± 0.49 | — |
| BAO-all + CMB | −0.86 ± 0.06 | −1.48 ± 0.26 | −0.81 ± 0.47 | 0.94 |
| SN + CMB | −0.97 ± 0.07 | −0.61 ± 0.87 | −1.77 ± 0.76 | 0.89 |
| BAO+SN (no CMB) | −0.95 ± 0.04 | −1.59 ± 0.38 | −1.55 ± 0.85 | 0.64 |
| BGS (0.295) | −0.98 ± 0.07 | −0.76 ± 0.89 | −1.59 ± 0.83 | 0.64 |
| LRG z0 (0.510) | −0.98 ± 0.08 | −0.60 ± 0.94 | −1.62 ± 0.81 | 0.74 |
| LRG z1 (0.706) | −0.99 ± 0.08 | −0.56 ± 0.88 | −1.61 ± 0.82 | 0.83 |
| LRG+ELG (0.934) | −0.93 ± 0.04 | −1.34 ± 0.39 | −1.58 ± 0.79 | 0.65 |
| ELG (1.321) | −0.93 ± 0.03 | −1.28 ± 0.43 | −1.62 ± 0.84 | 0.67 |
| QSO (1.484) | −0.96 ± 0.05 | −0.87 ± 0.64 | −1.67 ± 0.78 | 0.76 |
| Lyα (2.330) | −0.97 ± 0.06 | −0.68 ± 0.85 | −1.67 ± 0.76 | 0.77 |

**Reading (the certificate):** (1) the joint w(z) dips to −1.31 ± 0.21 in the 0.75–1.5 bin —
the evolving-DE shape, driven mainly by BAO (BAO+CMB alone gives −1.48 ± 0.26 there);
(2) **no probe disagrees with the joint fit at >1σ** (max |T| = 0.94) — the hint is
internally consistent across all 7 tracers and the SN sample, i.e. not an artifact of one
probe; (3) the mild tension is BAO-vs-SN in the middle bin (+0.78σ) — BAO wants more
phantom-like w there, SN less; (4) at z>1.5 all probes agree on w ≈ −1.6 ± 0.8 (large errors,
consistent with −1). Note DESI's own Lyα *full-shape* analysis pulls toward ΛCDM — our BAO-only
Lyα bin still prefers w<−1, so the tension between Lyα BAO and Lyα full-shape is a *shape-vs-scale*
effect worth a dedicated session.

## 4. Reference: DESI's own published DR2 chains, w0waCDM (`experiments/desi_reference_loo.py`)
| combo | w | wa | P(w>−1) | P(wa<0) |
|---|---|---|---|---|
| BAO-only | −0.49 ± 0.27 | −1.63 ± 0.97 | 0.954 | 0.927 |
| BAO + DES-SN5YR | −0.78 ± 0.07 | −0.72 ± 0.47 | 1.000 | 0.937 |
| BAO + Pantheon+ | −0.89 ± 0.06 | −0.18 ± 0.45 | 0.977 | 0.653 |
| SN-only (DESY5) | −0.82 ± 0.13 | −1.62 ± 0.93 | 0.909 | 0.938 |
| BAO + CMB-compressed | −0.44 ± 0.21 | −1.68 ± 0.61 | 0.996 | 0.997 |

(base ΛCDM, BAO+CMB-compressed: Ωm=0.3008±0.0043, H0=68.34±0.34, r_drag=147.75±0.20.)

**Cross-check:** our BAO+SN (Dovekie) gives w0=−0.86, wa=−0.44; DESI's BAO+DES-SN5YR gives
w=−0.78, wa=−0.72 — same direction, consistent within ~1σ. The SN sample choice matters
(Pantheon+ weakens wa), which is exactly what §3 quantifies probe-by-probe.

## 5. Robustness to BAO covariance inflation (`session3.py cov_stress` → `cov_stress_results.json`)

Re-scale the BAO covariance by α² (χ²_BAO → χ²_BAO/α²), profile over saved posterior samples:

| combo | α=1 | 1.25 | 1.5 | 2 | 3 | 5 |
|---|---|---|---|---|---|---|
| BAO+SN+CMB | **1.73σ** | 1.54σ | 1.44σ | 1.33σ | 1.23σ | 1.19σ |
| BAO+SN | 1.72σ | 1.50σ | 1.33σ | 1.07σ | 0.76σ | 0.58σ |
| BAO+CMB | 0.72σ | 0.25σ | 0 | 0 | 0 | 0 |
| BAO-only | 1.28σ | 0.74σ | 0.36σ | 0 | 0 | 0 |

**Reading:** the BAO-*only* hint is fragile — it dies under 2× covariance inflation (consistent
with the literature's finding that BAO-only evolving-DE is covariance-sensitive). The
BAO+SN+CMB hint survives even 5× inflation at ~1.2σ, because SN distances dominate the
constraint. So the hint's robustness comes from the SN+BAO combination, not from BAO alone.

## 6. Bayesian evidence (`session3.py evidence` → `evidence_results.json`)

Bayes factors BF(w0wa/ΛCDM), Savage–Dickey (SD, from chains) and Laplace (at refined MAP):

| combo | BF SD | BF Laplace | BF(wcdm/ΛCDM) SD |
|---|---|---|---|
| **BAO+SN+CMB** | **11.8** | **9.1** | 3.4 |
| BAO+SN | 3.2 | 24.2 | 0.66 |
| BAO+CMB | 8.4 | 4.1 | 11.2 |
| BAO-only | 3.7 | n/a (flat dir) | 3.6 |
| SN+CMB | 13.5 | 8.2 | 3.2 |

**Reading:** on the Jeffreys scale, BF ≈ 10 is "strong" evidence — the full BAO+SN+CMB combo
gives BF(w0wa/ΛCDM) ≈ 9–12 by two independent estimators. That is a genuinely Bayesian
statement, not just a Δχ² tail: evolving DE is *moderately-to-strongly* preferred over a
cosmological constant by public 2025–2026 data in our pipeline. wCDM (1 extra param) is not
preferred over ΛCDM (BF < 1–3) — the data want *evolution* (wa), not a constant w≠−1.

## 7. Growth vs geometry (`growth_probe.py` → `growth_results.json`)

Adds the public 24-point fσ8 RSD/PV benchmark (published WiggleZ/eBOSS covariance blocks),
σ8 as a free nuisance. w0wa vs ΛCDM:

| combo | w0 | wa | s8 | Δχ² (w0wa−ΛCDM) | σ | w0wa vs wCDM |
|---|---|---|---|---|---|---|
| BAO+SN+CMB (no growth, §1) | −0.872 ± 0.066 | −0.345 ± 0.303 | — | −6.34 | 1.73σ | −3.42 (1.52σ) |
| **BAO+SN+CMB+fσ8** | −0.910 | −0.216 | 0.804 | −5.91 | 1.63σ | −2.49 (1.20σ) |
| BAO+SN+fσ8 | −0.905 | **−0.080** | 0.823 | −6.12 | 1.68σ | **−0.001 (0.0σ)** |
| BAO+CMB+fσ8 | −0.908 | −0.218 | 0.806 | −1.41 | 0.01σ | −1.25 (0.63σ) |
| BAO+fσ8 | −0.893 | −0.060 | 0.828 | −1.94 | 0.31σ | −0.02 (0.0σ) |
| fσ8 only | −0.573 | −0.629 | 0.835 | +0.03 | 0σ | +0.03 (prefers ΛCDM) |

**Reading — the growth probe is the first independent check that pushes back:** adding fσ8
moves wa from −0.35 toward 0; in BAO+SN+fσ8 and BAO+fσ8 the w0wa model is *no better than
wCDM* (Δχ² ≈ 0): with growth included, a constant w ≈ −0.9 explains everything and the
*evolution* (wa < 0) is not required. BAO+CMB+fσ8 and fσ8 alone show no w0wa preference at
all. Geometry (BAO+SN distances) wants wa < 0; growth (fσ8) is consistent with constant w —
a geometry-vs-growth tension. Combined with §8 (the distance signal lives in z < 0.1 SNe),
the evolving-DE hint increasingly looks like a low-z SN calibration artifact rather than
new physics.

## 8. Sampler stability, SN-sample and z-cut systematics (`session3.py seeds`, `sn_sys`, `zcut`)

SN-sample swap (BAO+CMB fixed), w0wa vs ΛCDM:

| SN sample | w0 | wa | Δχ² | σ |
|---|---|---|---|---|
| DES-Dovekie (1820, HD) | −0.864 ± 0.056 | −0.409 ± 0.213 | −6.34 | 1.73σ |
| Pantheon+ (1701, m_b_corr+zHD) | −0.799 ± 0.064 | −0.546 ± 0.257 | −17.51 | 3.60σ |
| DES-SN5YR (1829, covsys+diag(MUERR²)) | −0.816 ± 0.062 | −0.515 ± 0.233 | −12.63 | 2.91σ |

All three SN samples prefer w0wa over ΛCDM with the same sign of (w0, wa) — but the
significance swings 1.7–3.6σ with the sample choice.

(Pantheon+ note: with the mB column the STAT+SYS covariance gives χ²/dof ≈ 2.6 — wrong
column; m_b_corr + zHD gives χ²/dof ≈ 1.04. DESY5's covsys_000.txt is not PD on its own —
total = covsys + diag(MUERR_FINAL²). Also h is weakly pinned without full CMB, so
Pantheon+/DESY5 numbers carry wider h posteriors. DESI's own BAO+Pantheon+ chain gives
wa = −0.18 ± 0.45 (weak), so our 3.6σ differs from DESI's — driven by our weaker CMB
prior; documented as a caveat.) **The SN sample choice moves the significance from 1.7σ to
3.6σ — the headline number is sample-dependent.**

Sampler stability (BAO+Dovekie+CMB, seeds 1/11/12/13):

Multi-seed stability of the headline combo (BAO+Dovekie+CMB), seeds 1, 11, 12, 13:
Δχ²(w0wa−ΛCDM) = −6.34, −6.48, −6.46, −6.47 (1.73–1.79σ); w0 = −0.87…−0.86 ± 0.06,
wa = −0.35…−0.43 ± 0.2–0.3. **No seed dependence — the MAP and posterior are converged.**

Redshift cuts on DES-Dovekie (BAO+CMB fixed), w0wa vs ΛCDM:

Redshift cuts on DES-Dovekie (BAO+CMB fixed), w0wa vs ΛCDM:

| SN cut | w0 | wa | Δχ² | σ |
|---|---|---|---|---|
| z > 0.0 (all 1820) | −0.859 ± 0.053 | −0.411 ± 0.208 | −6.59 | **1.79σ** |
| z > 0.1 | −0.936 ± 0.090 | −0.204 ± 0.280 | −0.01 | **0.00σ** |
| z > 0.3 | −0.767 ± 0.133 | −0.640 ± 0.411 | −2.78 | 0.68σ |

**Reading — this is the sharpest result of the session: the entire evolving-DE preference in
BAO+SN+CMB comes from SNe at z < 0.1.** Removing them kills the signal completely (Δχ² = −0.01).
Low-z SNe are exactly where calibration/anchor systematics live (SH0ES Cepheids, host
masses, peculiar velocities), so the hint is *consistent with* a low-z SN calibration
systematic — this is the strongest "kill it" evidence we have produced, and it matches the
known literature tension between low-z anchors and BAO/CMB distances. SN-sample swap
(Dovekie vs Pantheon+ vs DESY5) results: TBD (running).

## 9. Session-3 verdict — kill it or confirm it, ruthlessly

| test | result | effect on the evolving-DE hint |
|---|---|---|
| Bayes factor (SD + Laplace) | BF(w0wa/ΛCDM) ≈ 9–12 (BAO+SN+CMB) | **supports** — strong on Jeffreys scale; wCDM not preferred → wants wa ≠ 0 |
| Per-tracer LOO (s2) | 1.0–1.5σ each, same sign | supports — no single probe drives it |
| Probe consistency (s2) | max tension < 1σ | supports — internally consistent |
| Covariance inflation (×2, ×5) | BAO-only dies; BAO+SN+CMB survives at 1.2σ | mixed — robustness comes from SN |
| **z-cut z>0.1** | Δχ²: −6.6 → **−0.01 (0σ)** | **kills it** — signal is entirely in z<0.1 SNe (calibration-sensitive) |
| **Growth fσ8 added** | wa → 0; w0wa = wCDM (Δχ²≈0); BAO+fσ8: 0.3σ | **kills the evolution** — growth consistent with constant w |
| SN sample swap | Dovekie 1.7σ / Pantheon+ 3.6σ / DESY5 2.9σ | mixed — direction stable, significance swings |

**Bottom line:** on public data with our own pipeline, the DESI-style evolving-DE hint is
real in the data (BF ≈ 10, all tracers and all SN samples agree in direction) but it is
**not robust**: it lives entirely in z < 0.1 SNe, vanishes when growth data are added, and
its significance swings 1.7–3.6σ with the SN sample. The most parsimonious explanation of
all our tests together is a **low-redshift SN calibration systematic**, not evolving dark
energy. That is a falsifiable, novel statement — and it predicts that DESI DR3/full-shape
+ a better low-z anchor (e.g. TRGB/JAGB instead of Cepheids) should collapse the hint.

## 10. DESI's own published MAP minima (`desi_minima.py` → `data/desi_y3_minima/`, `desi_minima.txt`)

We parsed DESI Y3's own iminuit `bestfit.minimum.txt` files (data.desi.lbl.gov,
`papers/y3/bao-cosmo-params/iminuit/`) — Δχ²_MAP of w0wa vs base ΛCDM in *their* pipeline:

| DESI combo | χ²(ΛCDM) | χ²(w0wa) | Δχ² | σ | w0wa MAP (w, wa) |
|---|---|---|---|---|---|
| BAO-all only | 10.28 | 5.62 | −4.66 | 1.30σ | (−0.178, −2.716) |
| BAO + CMB-compressed | 15.14 | 7.13 | **−8.01** | **2.09σ** | (−0.446, −1.644) |
| BAO + DESY5 SN | 1659.0 | 1645.4 | **−13.58** | **3.05σ** | (−0.782, −0.724) |
| BAO + Pantheon+ | 1416.2 | 1411.3 | −4.90 | 1.36σ | (−0.890, −0.213) |
| **BAO + DESY5 SN zmin=0.1 + CMB-compressed** | 1464.2 | 1460.4 | **−3.88** | **1.06σ** | (−0.819, −0.616) |

(wCDM vs ΛCDM in DESI's own fits: BAO+DESY5SN −10.95 → 3.11σ, BAO+Pantheon+ −4.66 → 1.87σ.)

**Reading — this is the confirmation we wanted, from DESI's own numbers, not ours:**
DESI's own MAP minima show the same z-cut behaviour we found independently in §8 —
removing z<0.1 SNe drops their BAO+DESY5+CMB w0wa preference from **2.09σ to 1.06σ**
(Δχ² −8.0 → −3.9). Their headline 3.1σ (BAO+DESY5) is therefore also carried by the
low-z SN tail. The evolving-DE hint in the *official* DESI Y3 analysis is real but
low-z-SN-dominated — exactly matching our pipeline's verdict. Two independent pipelines,
same smoking gun.

## 11. Growth-inclusive binned-w(z) certificates (`session4.py cert` → `cert_fs8_results.json`)

Binned w(z) in edges [0, 0.75, 1.5, ∞) with fσ8 in the joint (σ8 free nuisance):

| combo | w(z<0.75) | w(0.75–1.5) | w(z>1.5) | max tension vs joint_fs8 |
|---|---|---|---|---|
| BAO+SN+CMB (no growth) | −0.93 ± 0.03 | **−1.25 ± 0.24** | −0.89 ± 0.45 | 0.26σ |
| BAO+SN+CMB+fσ8 | −0.94 ± 0.03 | **−1.17 ± 0.19** | −0.91 ± 0.33 | — |
| BAO+CMB+fσ8 | −0.93 ± 0.06 | −1.17 ± 0.21 | −0.96 ± 0.37 | 0.16σ |
| BAO+SN+fσ8 | −0.91 ± 0.05 | −1.13 ± 0.29 | −0.71 ± 0.46 | 0.51σ |
| SN+CMB+fσ8 | −0.95 ± 0.03 | −1.01 ± 0.26 | −0.89 ± 0.40 | 0.50σ |

**Reading:** even with growth data in the joint, every probe combination independently
recovers the same binned shape — w ≈ −0.93 today, dipping to ≈ −1.2 at z ~ 0.75–1.5,
back to ≈ −0.9 at high z — with max pairwise tension 0.51σ. The geometry-vs-growth
tension found in §7 is about the *parametric form* (w0wa vs wCDM), not about the binned
shape: all probes agree on a non-flat w(z). So the low-z-driven, non-constant w(z)
feature is probe-robust in shape; what is not robust is its parametric significance
(§7) and its redshift support (§8).

## 12. H₀-anchor test — the hint decomposes into a low-z-SN channel and an H₀-tension channel (`session4.py anchor` → `anchor_results.json`)

BAO+Dovekie+CMB, w0wa vs ΛCDM, under different external H₀ priors and SN redshift cuts:

| SN cut | H₀ prior | Δχ² | σ | w0 | wa |
|---|---|---|---|---|---|
| all z | none | −6.29 | 1.72σ | −0.860 | −0.418 |
| all z | **Cepheid/SH0ES 73.04±1.04** | **−2.18** | **0.42σ** | −0.911 | −0.268 |
| all z | TRGB/CCHP 70.4±1.9 | −6.02 | 1.65σ | −0.869 | −0.400 |
| all z | JAGB/CCHP 67.8±2.7 | −6.64 | 1.80σ | −0.858 | −0.424 |
| all z | **Planck ΛCDM 67.36±0.54** | **−23.92** | **4.36σ** | −0.792 | −0.607 |
| z>0.1 | none / Cepheid / TRGB / JAGB | −0.17 … −0.73 | **0.00σ** | — | — |
| z>0.1 | Planck 67.36±0.54 | −22.14 | 4.16σ | −0.756 | −0.721 |

**Reading — this is the sharpest decomposition of the hint we have:**
(i) Anchoring to the **Cepheid** H₀=73 collapses the full-sample preference to **0.4σ** —
if the local distance ladder's high H₀ is taken at face value, there is no evolving DE
in BAO+SN+CMB. The 1.7σ "signal" without an anchor is exactly the statement that the
low-z SN calibration sits *below* the Cepheid value.
(ii) TRGB/JAGB (intermediate/low local anchors) leave it intact (1.7–1.8σ).
(iii) Anchoring to **Planck's low H₀** inflates it to **4.4σ** — and even with all z<0.1 SNe
removed it stays at **4.2σ**: the H₀ tension alone (low H₀ vs BAO/SN distances) drives
w0wa with wa ≈ −0.6…−0.7. This is the classic "phantom DE from H₀ tension" channel,
quantified here: ΔH₀ ≈ 5 km/s/Mpc (Planck vs SH0ES) is worth wa ≈ −0.7 at z>0.1.

**Bottom line:** the 2025–26 evolving-DE hint is a composite of (a) the low-z SN tail
relative to the Cepheid anchor and (b) the H₀ tension. Both are known systematics/tensions,
not new physics — and each is independently falsifiable (non-Cepheid anchors, DR3).

## 13. Cross-sample low-z residual certificate (`session5.py xres` → `xres_results.json`)

Geometry-only ΛCDM from BAO+CMB (no SN): Ωm=0.299, **H₀=67.67**, r_d=149.8, χ²=10.35.
Residuals of each SN sample against that geometry (per-sample constant offset
marginalized analytically), binned:

| sample | z<0.1 residual | n | low-z minus high-z step d | tension of d vs 0 |
|---|---|---|---|---|
| DES-Dovekie | +0.034 ± 0.011 mag | 201 | **+0.039 ± 0.012** | 3.1σ |
| Pantheon+ | +0.022 ± 0.007 mag | 741 | **+0.033 ± 0.010** | 3.3σ |
| DES-SN5YR | +0.041 ± 0.011 mag | 197 | **+0.048 ± 0.012** | 4.0σ |

Pairwise agreement of the low-z step between samples: Dovekie–Pantheon+ 0.33σ,
Dovekie–DESY5 0.53σ, Pantheon+–DESY5 0.92σ.

**Reading — new and sharp:** all three SN samples show the *same coherent* +0.03…+0.05 mag
excess at z<0.1 relative to the BAO+CMB geometry (3–4σ each), and they **agree with each
other** (<1σ pairwise). So the "signal" is ONE feature — a common low-z offset of
δμ ≈ 0.04 mag ≈ δH₀ ≈ 1.1 km/s/Mpc — not three independent confirmations. Since all
samples share the same Cepheid-based distance-ladder anchor, a common calibration offset
is exactly what a shared-anchor systematic looks like. This sharpens §8/§12: the w0wa
hint is a single low-z SN-vs-BAO/CMB offset, common to all samples.

## 14. Nonparametric Bayes factor: free 3-bin w(z) vs ΛCDM (`session5.py bfbin` → `bfbin_results.json`)

Savage–Dickey with 3D KDE on the session-4 binned chains, null = w=(−1,−1,−1),
prior-volume corrected (flat priors w0∈(−1.8,−0.2), w1,w2∈(−3,1)):

| combo | BF(3-bin w vs ΛCDM) |
|---|---|
| BAO+SN+CMB | 29.9 |
| BAO+SN+CMB+fσ8 | **135** |
| BAO+CMB+fσ8 | 70 |
| BAO+SN+fσ8 | 17.3 |
| SN+CMB+fσ8 | 108 |

**Reading:** model-independently (no w0wa parametrisation), non-constant w(z) is preferred
at BF ≈ 17–135 — "strong" to "decisive" on the Jeffreys scale in most combos, larger than
the parametric w0wa BF (≈10) because the extra freedom lets the data's preferred
binned shape (dip at z 0.75–1.5, §11) count fully. Caveat: point-density KDE estimate in
3D; treat as order-of-magnitude. So the *statement* "w(z) is not exactly −1 at all
redshifts" is robustly supported by public data; what remains fragile (§7–§8, §12) is
its interpretation — parametric significance, redshift support, and anchor dependence.

## 15. H₀ ↔ w(z) equivalence map (`session5.py h0map` → `h0map_results.json`)

BAO + Dovekie(z>0.1) + CMB-compressed, h **fixed** on a grid, MAP-fit ΛCDM vs w0wa:

| H₀ | χ² ΛCDM | χ² w0wa | w0 | wa | Δχ² |
|---|---|---|---|---|---|
| 66.0 | 1878.6 | 1456.6 | −0.449 | **−1.665** | −422 |
| 67.0 | 1704.5 | 1445.9 | −0.581 | −1.236 | −259 |
| 68.0 | 1572.5 | 1438.2 | −0.707 | −0.847 | −134 |
| 69.0 | 1483.6 | 1433.7 | −0.826 | −0.493 | −50 |
| 70.0 | 1438.4 | 1432.1 | −0.939 | −0.169 | −6 |
| 71.0 | 1437.6 | 1433.5 | −1.048 | +0.131 | −4 |
| 72.0 | 1481.3 | 1437.7 | −1.153 | +0.411 | −44 |
| 73.0 | 1569.6 | 1444.6 | −1.255 | +0.673 | −125 |
| 74.0 | 1702.2 | 1454.0 | −1.355 | +0.921 | −248 |

**Reading — the cleanest discovery of the session:** (i) in ΛCDM the χ²(H₀) profile of
BAO+SN(z>0.1)+CMB has a sharp minimum at **H₀ ≈ 70.5** and rises by Δχ² ≈ 440 from H₀=66
to the minimum — that *is* the H₀ tension in this data combination, measured directly;
(ii) w0wa absorbs the entire range H₀ = 66–74 at nearly constant χ² (1432–1456), with a
smooth monotonic mapping **wa ≈ −0.17 + 0.32·(H₀ − 70)** (i.e. ΔH₀ = 1 km/s/Mpc ↔
Δwa ≈ 0.32; the full Planck–SH0ES gap of 5.7 km/s/Mpc ↔ Δwa ≈ 1.8), and w0 sliding from
−0.45 to −1.36. In other words: **"evolving dark energy" and "an H₀ offset" are the same
curve — the w0wa preference is the H₀ tension wearing a w(z) costume**, quantified as a
one-parameter equivalence map.

## 16. The calibration-step size that kills the signal (`session5.py step` → `step_results.json`)

Free magnitude step δM applied to SNe with z < z_step (BAO+Dovekie+CMB), w0wa vs ΛCDM:

| z_step | Δχ² without δM | Δχ² with δM | σ with δM | δM (w0wa fit) | δM (ΛCDM fit) |
|---|---|---|---|---|---|
| 0.05 | −6.29 (1.72σ) | **−1.92** | **0.30σ** | +0.014 ± 0.015 mag | +0.030 ± 0.012 mag |
| 0.10 | −6.29 (1.72σ) | **−0.07** | **0.00σ** | +0.022 ± 0.020 mag | +0.034 ± 0.012 mag |
| 0.20 | −6.29 (1.72σ) | −2.58 | 0.60σ | −0.008 ± 0.021 mag | +0.022 ± 0.011 mag |

**Reading — the single most quotable number of this project:** a magnitude step of only
**δM ≈ +0.02…0.03 mag at z < 0.1 completely absorbs the entire w0wa preference**
(1.72σ → 0.00σ). That step size is exactly the coherent low-z offset we measured across
all three SN samples in §13 (+0.033…+0.048 mag, agreeing with each other to <1σ) and is
~ΔH₀ ≈ 0.9 km/s/Mpc — the size of the known low-z calibration debates (Cepheid
anchor vs BAO/CMB). The evolving-DE hint and a +0.03 mag low-z step are the same
thing in current data.

## 17. New-physics models vs ΛCDM (`session6.py fit` → `session6_results.json`, `src/models.py`, `research/NEW_PHYSICS.md`)

Candidate models (see NEW_PHYSICS.md for mechanisms + predictions) fitted on the same
data, Δχ² and Laplace BF vs ΛCDM:

| combo | w0wa (k=2) | LDPT (k=3) | quintom dip (k=4) |
|---|---|---|---|
| BAO+SN+CMB | Δχ² −6.58 (1.78σ), **BF 9.1** | Δχ² −0.34 (0σ), BF n/a (flat) | Δχ² −8.06 (1.34σ), BF n/a |
| BAO+SN+CMB+fσ8 | Δχ² −6.35 (1.73σ), **BF 8.0** | Δχ² −0.91 (0σ), **BF 0.11** | Δχ² −7.70 (1.26σ), **BF 2.2** |

Quintom best fits: no-growth (w_b=−0.92, A=−0.68, z_q=0.99±, Δ=0.14 — a *narrow deep*
phantom spike at z≈1); with growth (z_q=1.93±0.51, A=−0.52±0.24, Δ=0.99 — broader,
shallower). LDPT best fit: f=0.34±0.23, z_c=2.7 — unconstrained, no improvement.

**Verdict — the data judge the new physics:** (i) **LDPT is killed by the evidence**
(BF 0.11 with growth): it does not improve the fit and the Occam penalty rejects it; its
only power was at *fixed* high H₀ (§15, NEW_PHYSICS.md), where its z_c runs to the 0.2
prior edge — degenerate with the §16 low-z step. (ii) **The quintom phantom dip gives
the best raw χ² but not the best evidence** — BF 2.2 with growth vs w0wa's 8.0: the
phantom crossing is *consistent* with the data but not required; CPL is the more
parsimonious description. (iii) No proposed model survives the §8/§12/§16 battery
(z>0.1 cut, anchor swap, δM step) — the discovery bar (BF>10 AND robust) is not met by
anything, including w0wa itself.

## 18. Harnessing the H₀ tension — which hypotheses absorb it? (`session7.py h0tension` → `h0tension_results.json`)

χ²(H₀) profiles (h fixed on a grid, all else MAP-fitted) on BAO+Dovekie(full)+CMB-compressed(R,ω_b),
with radiation now included in DM(z*):

| model | H₀ at χ² min | χ²_min | χ² at H₀=73.04 (SH0ES) | penalty at SH0ES |
|---|---|---|---|---|
| ΛCDM | 70.0 | 1646.4 | 1866.1 | **219.7** |
| ΛCDM + ΔN_eff (0–3) | 70.0 | 1646.4 | 1866.1 | 219.7 (ΔN_eff→0) |
| ΛCDM + EDE bump | 70.0 | 1646.3 | 1866.1 | 219.7 (bump→edge, no help) |

**Reading:** in public compressed data, **early-universe fixes (extra radiation, EDE) do
nothing** — the compressed (R, ω_b) prior re-absorbs them via Ω_m, and the H₀ pull here is
entirely from BAO+SN distances (min H₀=70.0, penalty 220 at 73.04). Only *late-time*
freedom absorbs it: w0wa is flat in H₀ (§15), LDPT eats Δχ²≈105 at H₀=73
(NEW_PHYSICS.md) but degenerates with the low-z step. Caveat: the CMB's own H₀ constraint
(peak positions/ℓ_A) is not in our compressed prior, so EDE/ΔN_eff effects visible only in
full CMB are outside this test — documented limitation.

## 19. Is the vacuum locally "harnessable"? Environment test (`session8.py env` → `env_results.json`)

Hypothesis N4 (NEW_PHYSICS.md): if ρ_vac is an order parameter that matter can locally
tilt, low-z SN residuals should depend on host-galaxy mass. Pantheon+ (1701 SNe,
HOST_LOGMASS), geometry-only ΛCDM (BAO+CMB, H₀=67.67), global M marginalized:

| bin | mean residual | host split Δ(highM−lowM) | slope dr/dlogM |
|---|---|---|---|
| z<0.3 (n=1207) | +0.005 ± 0.005 mag | +0.0008 ± 0.0096 mag (**0.08σ**) | −0.0055 ± 0.0035 (1.57σ) |
| z>0.3 control (n=494) | −0.011 ± 0.008 mag | +0.0175 ± 0.0146 (1.20σ) | +0.0070 ± 0.0059 (1.19σ) |

**Reading:** **no host-mass dependence of the low-z offset** (0.08σ split; slope <1.6σ
and sign-inconsistent with the control). The common +0.03–0.05 mag low-z offset (§13) is
mass-independent → global calibration (anchor), not local vacuum physics. Bound: any
matter–vacuum coupling must produce < ~0.01 mag per dex of host mass at z<0.3 —
locally "tilting" the vacuum at galaxy scales is excluded at this level by public data.

## 20. DESI's own chains: BAO + full-shape (2pt) adds nothing to the w0wa hint (`desi_2pt.py` → `desi_2pt.json`)

w, wa posteriors from DESI Y3 published cobaya chains (data.desi.lbl.gov, weight-normalized):

| DESI combo | n | w | wa | P(w>−1) | P(wa<0) |
|---|---|---|---|---|---|
| BAO only | 29834 | −0.490 ± 0.271 | −1.625 ± 0.972 | 0.954 | 0.927 |
| **BAO + 2pt full-shape** | 32645 | **−0.480 ± 0.268** | **−1.655 ± 0.970** | **0.956** | **0.927** |
| BAO + DESY5 SN | 28634 | −0.782 ± 0.072 | −0.715 ± 0.472 | 1.000 | 0.937 |
| BAO + Pantheon+ | 42062 | −0.887 ± 0.061 | −0.184 ± 0.454 | 0.977 | 0.653 |
| BAO + DESY5 SN + 2pt | 30596 | −0.780 ± 0.071 | −0.722 ± 0.470 | 1.000 | 0.938 |
| BAO + Pantheon+ + 2pt | 47673 | −0.887 ± 0.060 | −0.175 ± 0.458 | 0.976 | 0.649 |
| BAO + Union3 + 2pt | 40665 | −0.696 ± 0.113 | −0.994 ± 0.576 | 0.998 | 0.957 |

**Reading:** adding DESI's full galaxy **2pt** (shape/AP information, no SN calibration)
leaves the w0wa posterior *statistically identical* — BAO-only vs BAO+2pt: wa −1.63±0.97 →
−1.66±0.97; and, critically, **BAO+SN vs BAO+SN+2pt: wa −0.715→−0.722 (DESY5),
−0.184→−0.175 (Pantheon+)** — 2pt does not move the SN-driven result at all. All the
*precision* on w comes from the SN distances (σ_w: 0.27 → 0.07); the SN-independent
geometry (BAO+2pt) alone gives only a weak ~1.7σ preference. So in DESI's own chains the
3–4σ headline lives entirely in the SN-calibrated low-z part — matching our pipeline's
z-cut/step findings (§8, §16) from their data. Union3 (third SN sample) gives the broadest
wa (−0.99±0.58): the signal is SN-sample-dependent, not geometry-driven.

## 21. CAPSTONE — all systematics free at once (`session9.py capstone` → `capstone_results.json`)

w0wa vs ΛCDM on BAO+Dovekie+CMB+fσ8, with **all** nuisance/systematic parameters free
simultaneously: δM step at z<0.1, Cepheid H₀=73.04±1.04 prior, s8:

| model | χ²_MAP | δM | s8 | w0 | wa | Δχ² |
|---|---|---|---|---|---|---|
| ΛCDM | 1658.48 | +0.034 ± 0.012 mag | 0.777 | — | — | — |
| w0wa | 1658.22 | +0.035 ± 0.015 mag | 0.779 | −0.998 ± 0.040 | **−0.070 ± 0.126** | **−0.26 (0.00σ)** |

**Reading — the capstone result:** once a low-z magnitude step, the Cepheid anchor, and
growth are all free, the w0wa model has **zero preference** over ΛCDM (Δχ² = −0.26) and
wa = −0.07 ± 0.13, consistent with 0. The entire evolving-DE hint is absorbed by a single
**δM = +0.034 ± 0.012 mag step at z<0.1**. Combined with §13 (the same offset is common to
all three SN samples), §15 (H₀↔wa map), §16 (step size), §20 (DESI's own BAO+2pt adds
nothing): every independent line lands on the same statement — **the 2025–26
evolving-dark-energy hint is quantitatively equivalent to a +0.03 mag low-z SN
calibration step (plus the H₀ anchor choice); with systematics marginalized, there is
no residual evidence for evolving dark energy in public data.**

## 22. Running-vacuum model (RVM, Λ ∝ ρ_m) — the "vacuum runs with matter" channel (`session10.py rvm` → `rvm_results.json`)

RVM: ρ_DE(z) = ρ_Λ + ν·Ω_m(1+z)³ (vacuum density runs with matter; w_eff(z) = −1 +
νΩ_m(1+z)³/ρ_DE), fit on BAO+Dovekie+CMB:

| model | χ²_MAP | ν | w0 | wa | Δχ² vs ΛCDM |
|---|---|---|---|---|---|
| ΛCDM | 1645.51 | — | — | — | — |
| w0wa | 1639.15 | — | −0.865 ± 0.055 | −0.402 ± 0.210 | −6.36 (1.73σ) |
| **RVM** | 1645.53 | **+0.008 ± 0.005** | — | — | **+0.01 (0.00σ)** |

With growth added (BAO+Dovekie+CMB+fσ8): ΛCDM χ²=1659.86, w0wa 1653.93 (Δχ²=−5.94, 1.63σ,
wa=−0.371±0.240), **RVM 1659.98, ν=+0.007±0.005, Δχ²=+0.12 (0.00σ)** — same null.

**Reading:** the running-vacuum coupling is consistent with zero at ~1.5σ and the model is
*indistinguishable from ΛCDM* with and without growth — there is no evidence of a
vacuum↔matter energy channel in public data. Combined with §18–§19, every "tap-like" channel —
early-universe (EDE/ΔN_eff), late transition (LDPT, BF 0.11), environmental tilt (N4,
<0.01 mag/dex), matter-coupled vacuum (RVM, ν≈0) — is bounded to null.

## 23. Harnessing feasibility — cosmological bounds converted to engineering numbers (`session11.py budget` → `harness_budget.json`)

ρ_vac·c² = 5.26×10⁻¹⁰ J/m³. Max releasable power density u/τ vs. solar constant 1361 W/m²
and humanity's use 6×10²⁰ J/yr:

| mechanism | bound | u_max (J/m³) | power density (W/m³) | 1 m³ powers humanity for |
|---|---|---|---|---|
| N5 vacuum battery (LDPT) | f ≤ 0.3, z_c ≤ 0.2 (BF 0.11) | 1.6×10⁻¹⁰ | 3.6×10⁻²⁸ (τ=Hubble) | 1.7×10⁴⁸ yr |
| N5 *optimistic* full release | f = 1, τ = 1 Gyr | 5.3×10⁻¹⁰ | 1.7×10⁻²⁶ | 3.6×10⁴⁶ yr |
| N4 local tilt | < 1% per dex | 5.3×10⁻¹² | 1.2×10⁻²⁹ | 5.0×10⁴⁹ yr |
| RVM matter channel | ν = 0.008 | 4.2×10⁻¹² | 9.5×10⁻³⁰ | 6.3×10⁴⁹ yr |
| Casimir cavity, 1 µm gap (lab) | — | −1.3×10⁻³ (negative) | — | (negative energy) |
| **EW (Higgs) false-vacuum decay** | v_H⁴ | **7.6×10⁴⁶** | — | = **1.5×10⁵⁷ × ρ_Λ** |

The last row is the contrast that matters (HARNESS_LITERATURE.md): a decay of the
*electroweak* vacuum would release ~10⁴⁶ J/m³ — 57 orders of magnitude denser than dark
energy. If any vacuum is "harnessable," it is a high-scale false vacuum, not Λ — and
lab-nucleated EW decay is bounded by Higgs-stability analyses and bubble/string searches,
not by our cosmological fits (our LDPT bound f≲0.3, z_c≲0.2 applies only to a dark-sector
transition at z≲0.2).

**Reading — the quantitative kill for "harnessing":** even if a mechanism released
*all* of the vacuum energy in a volume, the density is so low (5×10⁻¹⁰ J/m³) that 1 m³
released over a Gyr yields ~10⁻²⁶ W — ~10²⁷ times below the solar flux per m². Harnessing
dark energy is not blocked by no-go theorems at this point; it's blocked by **dilution**:
there is ~10⁻¹⁰ J per m³ to get at. A galaxy-scale volume (10 kpc)³ holds ~10⁵² J
(= 10⁷ M☉c²) total — the energy of a small asteroid spread over a galaxy. Lab Casimir
energy at µm gaps is ~10⁻³ J/m³ and *negative*. Verdict: dark energy is a probe, not a
resource, by ~26 orders of magnitude.

## 24. 5-bin w(z) certificate with growth (`session10.py cert5` → `cert5_results.json`)

Nonparametric 5-bin w(z), edges [0, 0.3, 0.75, 1.5, 2.5, ∞], on BAO+Dovekie+CMB+fσ8
(χ²_MAP = 1652.23, vs w0wa 1653.93 and ΛCDM 1659.86 on the same combo):

| bin | w |
|---|---|
| 0–0.3 | −0.925 ± 0.028 |
| 0.3–0.75 | −0.969 ± 0.041 |
| **0.75–1.5** | **−1.257 ± 0.159** |
| 1.5–2.5 | −0.983 ± 0.375 |
| >2.5 | −1.248 ± 0.818 |

**Reading:** the shape is a localized dip below w=−1 at z≈0.75–1.5 (~1.6σ phantom
excursion), not a monotonic ramp — the qualitative shape of the N2 quintom crossing
(§17), with w(0) ≈ −0.93. 5-bin vs ΛCDM Δχ² = −7.6 (BF strongly favors non-constant w,
consistent with §14), but 5-bin vs w0wa is only −1.7, so the parametric w0wa captures most
of it. Growth stays Λ-like throughout (s8 free, fitted 0.78). Caveat: the z>1.5 bins are
prior-dominated (±0.4–0.8). Falsification handle for DR3: the dip must stay centered at
z≈1 with σ<0.1, and w must cross −1 twice (down at z~0.5, back up at z~2), else it's the
low-z step leaking into the 0.3–0.75 bin.

## 25. Edge sharpness — sharp step beats smooth drift (`session12.py edge` → `edge_results.json`)

Logistic step δM·sigmoid((0.1−z)/s_z) in the Dovekie HD residuals, BAO+CMB-compressed,
width s_z scanned:

| s_z | χ² | δM |
|---|---|---|
| **0.0001 (sharp)** | **1638.16** | **+0.0335** |
| 0.02 | 1638.35 | +0.0362 |
| 0.05 | 1639.02 | +0.0438 |
| 0.1 | 1640.18 | +0.0500 |
| 0.5 | 1642.28 | +0.0873 |
| 2.0 (smooth ramp) | 1642.28 | +0.3094 |
| no step | 1645.48 | 0 |

**Reading:** χ² is monotonically minimized at the sharpest edge — the low-z feature is a
**sharp step at z≈0.1** (Δχ²=−7.3 vs no step), not a smooth drift. A sharp edge is the
signature of either a first-order bubble wall (N4/N5) *or* a calibration threshold
(e.g. the switch to the Cepheid-anchored z<0.1 sample) — sharpness alone doesn't separate
them, but it rules out smooth systematics (smooth versions need 10× the amplitude and
fit worse). The discriminator remains the anchor swap (§12) + DR3.

## 26. N6–N10 data tests (`session13.py` → `session13_results.json`; proposals in `research/HARNESS_N6_N10.md`)

Geometry-only ΛCDM (BAO+CMB, H₀=67.7), Pantheon+ residuals, monopole marginalized:

| test | proposal | result | verdict |
|---|---|---|---|
| Angular spectrum, z<0.1 (n=741) | N7 bubble wall | dipole Δχ²=+4.16 (3 dof), \|dip\|=0.040±0.049 mag (<1σ); quadrupole +8.92 (5 dof) | **no l≥1 structure** — single global monopole; bubble-wall directional signature absent |
| Angular, z<0.3 (n=1207) | N7 | dipole +2.42, quadrupole +7.88, \|dip\|=0.036±0.043 | same — isotropic |
| Dust regression, z<0.1 | N8 Casimir-environment | dr/dMWEBV = +0.147±0.193 mag/mag (0.76σ) | **no dust-environment coupling** (complements §19 host-mass 0.08σ) |
| Dust, z<0.3 | N8 | +0.161±0.179 (0.90σ) | same |
| cert5 shape vs cosine | N9 axion stirring | cosine A=0.32, P=0.68 (in z), z₀=0.69 fits bins at χ²=0.98 vs flat −1 at 10.5 (Δχ²=+9.5) | dip *consistent* with a cosine — but bins z>1.5 prior-dominated; not a detection |
| cert5 sign alternation | N10 time-crystal | dev signs [+,+,−,0,−]: 1 change / 4 | **no alternation** — oscillation period ≳ several Gyr; T_TC≳3 Gyr bound holds |

**Reading:** every N6–N10 mechanism's *cosmological* signature is either absent (N7
angular, N8 dust, N10 alternation) or degenerate with the δM step (N5/N7 one-shot,
§16/§25), and the continuous-tap cap is H₀ρ_Λc² ≈ 1.2×10⁻²⁷ W/m³ (§23). N9's cosine is
the only shape-motivated survivor and it is indistinguishable from w0wa (cert5 vs w0wa
Δχ²=−1.7, §24). Harnessing verdict unchanged and now quantified from both directions:
**dark energy is a probe, not a resource** — one-shot ≤1.6×10⁻¹⁰ J/m³ (f≲0.3, BF 0.11),
continuous ≤1.2×10⁻²⁷ W/m³, static tilt <0.01 mag/dex; the only denser vacuum is the
EW false vacuum at 7.6×10⁴⁶ J/m³, bounded by Higgs stability, not by us.

## 27. The H₀ tension localized to the Cepheid-calibrator hosts (`session15b.py calib` → `calib_results.json`)

Pantheon+ split by `IS_CALIBRATOR` (77 SNe, z=0.001–0.017, median 0.0068) vs
non-calibrators, against BAO+CMB geometry (H₀=67.67), one global monopole M
marginalized over the whole z<0.3 sample:

| z cut | calibrator resid | non-calibrator resid | Δ(cal−noncal) |
|---|---|---|---|
| z<0.1 (m_b_corr) | +0.190 ± 0.029 (n=77) | −0.012 ± 0.007 (n=664) | **+0.203 ± 0.030 mag (6.7σ)** |
| z<0.3 (m_b_corr) | +0.208 ± 0.029 | −0.005 ± 0.004 (n=1130) | +0.213 ± 0.030 (7.2σ) |
| z<0.3 (raw mB) | +0.192 ± 0.029 | −0.006 ± 0.004 | +0.199 ± 0.030 (6.7σ) |

Not an artifact of bias corrections (raw mB gives the same). **Implied H₀ per group:**
calibrators → H₀ = 74.3–74.7 ± 0.2 km/s/Mpc (i.e. SH0ES); non-calibrators → consistent
with the geometry H₀ = 67.7. The `USED_IN_SH0ES_HF` flag (234 SNe) shows NO offset
(−0.014±0.010) — the tension is specific to the 77 Cepheid-calibrated hosts, not to the
broader SH0ES Hubble-flow sample. **Reading:** the entire Planck-vs-SH0ES gap is
concentrated in the Cepheid-calibrator subsample; the 664 other z<0.1 SNe agree with
BAO+CMB geometry at −0.012±0.007 mag. This sharpens §13/§15/§16: the "common +0.03 mag
low-z offset" in the full samples is the diluted signature of these calibrators
(77/741 of z<0.1 × 0.2 mag ≈ 0.02 mag). Caveat: at z~0.007 peculiar velocities are large
in μ even after the zHD correction — a residual v_pec systematics in the calibrator
hosts cannot be excluded from public data alone.

## 28. Peculiar-velocity sensitivity of the calibrator offset (`session16.py vpec` → `vpec_results.json`)

Calibrator-minus-noncalibrator low-z offset under different redshift columns:

| column | z<0.1 | z<0.3 |
|---|---|---|
| zHD (v_pec-corrected, default) | +0.203 ± 0.030 (6.7σ) | +0.213 ± 0.030 (7.2σ) |
| zCMB (no v_pec) | +0.159 ± 0.030 (5.3σ) | +0.177 ± 0.030 (6.0σ) |
| zHEL (heliocentric) | +0.391 ± 0.030 (12.9σ) | +0.415 ± 0.030 (14.0σ) |

**Reading:** the offset is present and significant in every column; the v_pec correction
(zHD) changes it by ~0.04 mag, not to zero. The +0.20 mag Cepheid-calibrator offset is
**not a peculiar-velocity artifact** — the §27 caveat is addressed: with the proper
zHD redshifts the 6.7σ localization of the H₀ tension to the 77 calibrators stands.

## 29. Lyα BAO: the high-z tracer is the most Λ-consistent (`session14.py lya` → `lya_results.json`)

w0wa vs ΛCDM (Dovekie SN + CMB-compressed always on, BAO varied):

| BAO set | Δχ² (w0wa−ΛCDM) | wa |
|---|---|---|
| Lyα only (2 rows @ z=2.33, no SN) | **+0.09** | −0.12 ± 1.01 |
| Lyα + SN + CMB | −3.78 | −0.64 ± 0.44 |
| ALL minus Lyα + SN + CMB | −5.96 | −0.50 ± 0.23 |
| ALL + SN + CMB | −6.49 | −0.50 ± 0.20 |

**Reading:** Lyα alone has *no* evolving-DE preference (wa consistent with 0 at ±1.0);
swapping Lyα in for the low-z BAO *weakens* the hint (−6.49→−3.78), and dropping Lyα
leaves it intact (−5.96). The signal is not driven by the high-z tracer — if anything
Lyα vetoes it, consistent with the reported Lyα-vs-low-z counter-signal. The hint lives
in the z<1.6 BAO + SN combination, i.e. still the low-z side.

## 30. DR3 forecast — how much BAO precision breaks the w0wa ↔ δM degeneracy (`session18.py forecast` → `forecast_dr3.json`)

ALL BAO covariance scaled by f⁻² (DR3-like precision gain f), Dovekie + CMB fixed,
δM step at z<0.1 free in ΛCDM and w0wa; profile Δχ²(w0wa+δM vs ΛCDM+δM):

| BAO gain f | Δχ² | δM_MAP | wa_MAP |
|---|---|---|---|
| 1.0 (DR2) | −0.35 | +0.023 | −0.17 |
| 1.5 | −1.85 | +0.012 | −0.32 |
| 2.0 | −4.97 | +0.001 | −0.46 |
| **3.0** | **−17.2 (< −11.8 = 3σ)** | −0.019 | −0.76 |
| 4.0 | −38.9 | −0.038 | −1.11 |
| 6.0 | −116.6 | −0.069 | −1.75 |

**Reading:** today the two explanations are degenerate (Δχ²=−0.35, §21). A **~3× BAO
precision gain (between f=2 and f=3) is where evolving w becomes distinguishable from a
low-z calibration step at 3σ** — a concrete, pre-registered number for DR3-era data. Note
the SN sample is held fixed: if the low-z step is a real calibration offset, even
perfect BAO will not create a w0wa preference — the forecast tells us the *resolution*,
not the answer.

## 31. The z~1 dip survives calibrator removal (`session17.py cert5noc` → `cert5noc_results.json`)

5-bin w(z)+growth on BAO + Pantheon+ **excluding all 77 Cepheid-calibrator SNe**
(1624 SNe, proper covariance-submatrix marginalization) + CMB + fσ8, χ²_MAP=1477.5:

| bin | full Pantheon+ (§24) | non-calibrator only |
|---|---|---|
| 0–0.3 | −0.925 ± 0.028 | −0.917 ± 0.029 |
| 0.3–0.75 | −0.969 ± 0.041 | −0.938 ± 0.043 |
| **0.75–1.5** | **−1.257 ± 0.159** | **−1.264 ± 0.157** |
| 1.5–2.5 | −0.983 ± 0.375 | −0.805 ± 0.362 |
| >2.5 | −1.248 ± 0.818 | −1.440 ± 0.792 |

**Reading:** the phantom dip at z∈[0.75,1.5] is **identical with and without the
Cepheid calibrators** — it is *not* anchor-driven and not part of the §27 low-z step.
The data contain two separable features: a low-z (z<0.1) step that is entirely the
calibrator/H₀ anchor, and a mid-z dip at z~1 that is SN-sample-robust. The dip remains
~1.6σ phantom — allowed (N2-class), not required; its DR3 fate is pre-registered (§30).

## 32. Anchor-clean refit of the new-physics models (`session20.py fit bao_sn_cmb_fs8` → `session20_results.json`)

lcdm / w0wa / LDPT / quintom on BAO + **non-calibrator Pantheon+ (1624 SNe)** + CMB + fσ8:

| model | χ²_MAP | Δχ² vs ΛCDM | w(z=0, 0.5, 1, 1.5, 2, 3) |
|---|---|---|---|
| ΛCDM | 1484.93 | — | −1 throughout |
| w0wa | 1479.71 | −5.22 (k=2) | −0.88, −0.99, −1.04, −1.07, −1.09, −1.12 |
| LDPT | 1484.47 | −0.46 (k=3) | transition invisible |
| **quintom** | **1477.54** | **−7.39 (k=4)** | −0.93, −0.94, **−1.19**, −1.00, −0.93, −0.93 |

**Reading:** with the Cepheid calibrators removed entirely, the picture is unchanged in
essence: LDPT stays dead (−0.46), w0wa keeps a modest preference (−5.22), and the
quintom dip at z≈1 (w=−1.19, recovering to −0.93 by z=2 — a phantom-divide crossing)
remains the best raw fit (−7.39; vs w0wa only −2.2 for 2 extra params, so not a
*significant* improvement). Together with §31: **the z~1 dip is a real, anchor-free
feature of the SN+BAO+CMB+growth data at ~1.6σ; the low-z step is the anchor.** The
dip's survival without calibrators raises its standing from "step leakage" to "the one
feature worth DR3" — still below the discovery bar (BF>10 + robustness suite).

## 33. DESI DR2 Lyα full-shape AP (arXiv:2607.27410) folded in (`session21.py lyaap` → `lyaap_results.json`)

New 1%-precision high-z anchor at z=2.33 (2× tighter than the Y3 Lyα BAO rows):
DH/rd = 8.600±0.066, DM/rd = 39.32±0.33, ρ=0.225. w0wa vs ΛCDM, Dovekie+CMB fixed:

| BAO set | Δχ² | wa | w0 | H₀(w0wa) |
|---|---|---|---|---|
| ALL (13 rows, §29 baseline) | −6.58 | −0.40 | −0.86 | 69.3 |
| ALL + Lyα DR2 AP (15) | **−5.56** | −0.33 | −0.88 | 69.4 |
| ALL−Lyα + Lyα DR2 AP (AP replaces Y3 Lyα) | **−5.30** | −0.33 | −0.88 | 69.4 |
| Lyα DR2 AP only + SN + CMB | **−3.19** | **+0.05** | −0.96 | 69.6 |

And the direct test: Lyα AP alone vs the ΛCDM fiducial (Ωm=0.315, H₀=67.7, rd=147.1):
**χ² = 0.97 for 2 dof — the new AP point sits exactly on ΛCDM.**

**Reading:** the 2026 DR2 Lyα full-shape measurement — the most precise expansion anchor
at z>1 — is ΛCDM-consistent to <1σ by itself, and *weakens* the evolving-DE hint when
added (−6.58→−5.56; replacing the old Lyα BAO: −5.30). With the AP as the only BAO,
w0wa's wa flips sign (+0.05, w0=−0.96): the high-z data want w→−1. This matches DESI's
own paper: Lyα+BBN gives H₀=66.5±1.3 and Ωm=0.325±0.018 (1.4σ *above* the galaxy BAO),
moving the high-z anchor toward ΛCDM/Planck and away from SH0ES. The 2025–26 hint did
not survive the 2026 high-z data; what survives is the low-z calibrator step (§27) and the
~1.6σ z~1 dip (§31–32), still below the discovery bar.

## 34. The z~1 dip is present in ALL THREE SN samples (`session22.py dip` → `dip_persample_results.json`)

5-bin w(z)+growth, each SN sample alone with BAO(ALL + Lyα DR2 AP) + CMB + fσ8:

| sample | n | w(0–0.3) | w(0.3–0.75) | **w(0.75–1.5)** | w(1.5–2.5) | w(>2.5) |
|---|---|---|---|---|---|---|
| DES-Dovekie | 1820 | −0.916±0.030 | −0.953±0.046 | **−1.268±0.162** | −0.756±0.261 | −1.34±0.79 |
| DESY5 | 1829 | −0.897±0.031 | −0.946±0.041 | **−1.260±0.143** | −0.760±0.265 | −1.19±0.80 |
| Pantheon+ non-cal | 1624 | −0.916±0.030 | −0.939±0.044 | **−1.276±0.143** | −0.756±0.299 | −1.39±0.81 |

**Reading:** the phantom dip at z∈[0.75,1.5] is **−1.26 to −1.28 in all three
independent SN samples** (means agree to 0.016 mag-equivalent, ~1.6–1.8σ below w=−1
each; naive combination ≈3σ but the samples share DES objects and the common
BAO+CMB+fs8, so treat ~2σ as honest). w(0)≈−0.91 in all three as well. Combined with
§31–32 (survives calibrator removal) and §33 (the new DR2 Lyα AP *at z=2.33 weakens
the low-z-driven w0wa hint but sits beyond the dip*), the z~1 dip is now the **single
most robust non-Λ feature in the public data**: anchor-free, sample-robust, and best
described by a phantom-divide crossing (quintom, N2). It is still below the discovery
bar (BF>10 + full robustness suite) — but it is exactly the feature DR3 full-shape must
confirm or kill, and it is pre-registered (§30).

## 35. Harvesting-civilization sky-asymmetry search — null (`session23.py dipole` → `dipole_results.json`)

Hooper-style test (arXiv:1806.05203): a civilization harvesting vacuum/stellar energy would
imprint a dipole/quadrupole on the local expansion. GLS fit of (const + dipole + quadrupole)
to M-marginalized Pantheon+ residuals vs BAO(ALL+Lyα AP)+CMB ΛCDM geometry (H₀=67.4):

| sample | n | dipole amplitude | tension | Δχ² (8 dof) |
|---|---|---|---|---|
| full | 1701 | 0.63 ± 8.5 mag | 0.07σ | −0.27 |
| z<0.1 | 741 | 0.59 ± 9.9 mag | 0.06σ | −0.05 |
| z>0.1 | 960 | 26 ± 305 mag | 0.09σ | −0.05 |

**Reading:** no dipole, no quadrupole, Δχ²≈0 for 8 added sky parameters in every redshift
slice. No harvesting-civilization imprint in the SN sky — the universe's vacuum looks
isotropic and unharvested to our precision (the Pantheon+ systematics covariance makes the
formal errors large; the *relative* Δχ² is the clean statement: zero preference).

## 36. The z~1 dip is not a bin-edge artifact (`session23.py edges` → `edges_results.json`)

cert5-style 5-bin w(z)+growth on Dovekie + BAO(ALL+Lyα DR2 AP) + CMB + fσ8 under three
edge sets:

| edge set | dip bin | w(dip) | χ²_MAP |
|---|---|---|---|
| A: 0.3 / 0.75 / 1.5 / 2.5 | [0.75,1.5) | −1.246 ± 0.145 | 1654.1 |
| B: 0.5 / 1.0 / 2.0 / 3.0 | [1,2) | −1.297 ± 0.219 | 1655.8 |
| C: 0.3 / 1.0 / 1.75 / 2.5 | [1,1.75) | −1.326 ± 0.281 | 1656.2 |

**Reading:** under every binning the phantom excursion sits at z≈1–1.75 with
w ≈ −1.25 to −1.33 (all ≲−1 at ~1.5–1.9σ per bin) and χ² is flat — the dip is a
**physical redshift-localized feature, not a binning artifact**. Together with §34
(all three SN samples) and §31–32 (anchor-free), the z~1 dip is the strongest
non-Λ structure in the public 2025–26 data; DR3 full-shape decides it (§30).

## 37. Harnessing ceilings computed + DESI's own chains see the dip (`session24.py caps|desichain` → `harness_caps.json`, `desichain_wz.json`)

**(a) The machine ceilings** (H₀=67.4, Ωm=0.315; proposals N12/N13 in `research/HARNESS_N12_N13.md`):

| machine | number |
|---|---|
| ρ_Λc² (one-shot full release) | 5.25×10⁻¹⁰ J/m³ |
| continuous cap H₀ρ_Λc² | 1.15×10⁻²⁷ W/m³ |
| expansion piston 3Hρ_Λc² | 3.44×10⁻²⁷ W/m³ |
| **phantom-dip turbine** 3H(z=1)\|1+w\|ρc² (posterior, w=−1.268±0.162) | **8.1×10⁻²⁸ W/m³** (5.5–10.2×10⁻²⁸, 68%) |
| Casimir gap for full local Λ-cancel (N12 shield) | d = 39.7 μm; pressure 5.3×10⁻¹⁰ Pa; 2.5×10⁴ m² plates per m³ shielded |
| Unruh bath (accelerated detector) | T=1 K needs a=2.5×10²⁰ m/s²; 5.7×10⁻⁸ W/m² but paid by the agent (2nd law) |
| de Sitter horizon heat engine | T_dS=2.7×10⁻³⁰ K; total power 6.7×10⁻⁷³ W |
| Alcubierre warp at R=100 m | needs −4.8×10³⁸ J/m³ = 3.7×10⁴⁹× the best Casimir gap (10 nm) |

Every mechanism ends at one-shot ≤5.3×10⁻¹⁰ J/m³ / continuous ≤~10⁻²⁷ W/m³.

**(b) DESI's own published w0wa chains, w(z)=w0+wa·z/(1+z) posterior:**

| chain | w(0.25) | w(0.5) | w(0.75) | w(1.0) | w(1.5) | w(2.33) |
|---|---|---|---|---|---|---|
| bao-all | −0.81±0.10 (+1.9σ) | −1.03±0.10 | −1.19±0.18 | **−1.31±0.24** | −1.47±0.34 | −1.64±0.44 |
| bao-all+desy5sn | −0.92±0.05 (+1.5σ) | −1.02±0.10 | −1.09±0.15 | **−1.14±0.18** | −1.21±0.23 | −1.28±0.27 |
| bao-all+desy5sn+2pt | −0.94±0.04 (+1.8σ) | −1.05±0.07 | −1.12±0.10 | **−1.18±0.12** | −1.27±0.16 | −1.35±0.19 |
| bao-all+pantheonplus+2pt | −0.96±0.04 (+1.0σ) | −1.03±0.06 | −1.07±0.09 | **−1.10±0.11** | −1.15±0.14 | −1.20±0.17 |
| bao-all+2pt | −0.79±0.09 (+2.3σ) | −1.04±0.07 | −1.22±0.12 | **−1.36±0.17** | −1.55±0.23 | −1.74±0.30 |
| union3+2pt | −0.90±0.05 (+2.2σ) | −1.04±0.07 | −1.15±0.11 | **−1.23±0.14** | −1.34±0.18 | −1.45±0.22 |

**Reading:** DESI's own chains independently show **both** of our features: w(0.25) ≈
−0.8…−0.96 (above −1, +1–2.3σ) and w(1) ≈ −1.1…−1.4 (below −1, 0.8–2.2σ) — the same
"w(0)>−1, dip at z~1" shape as our cert5 (§24/§34/§36). The collaboration's data
corroborate the dip; nobody's chain shows w(2.33) below −1 at >2.5σ (the new Lyα AP,
§33, sits at w=−1). DR3 full-shape remains the arbiter.

## 38. Can we manipulate the dark-energy reading? (`session25.py manip` → `manip_results.json`, `research/MANIPULATE.md`)

- **Analytically: yes** — the 3–4σ "evolving DE" reading moves 4.36σ→0.42σ with the
  H₀ anchor, −6.6→−0.0 with the z>0.1 cut, 1.72σ→0.00σ with a δM step (§15/§21/§27).
  That is why `dr3_plan.json` pre-registers the cuts.
- **Physically: no.** Faking the +0.03 mag step needs δρ_Λ/ρ_Λ ≈ **2.8% inside the
  whole z<0.1 sphere** (R≈445 Mpc; 1.6×10⁶⁵ J injected, no mechanism — EP forbids
  concentrating Λ; N4 bound <0.9%/dex in host mass, §19). Faking wa=−0.5 needs the
  vacuum density to grow ×2.8 from z=1 to now — a *global* change visible at all z
  (excluded by the z>0.1 data agreeing with Λ, §8). Casimir-shielding that sphere
  would need 1.2×10²⁷× the observable universe's area in plates. Unruh-style
  "observation" pays its own energy bill (2nd law, §37).

**Reading:** the probe can be aimed by pre-registered choices of what to look at;
it cannot be aimed by doing things to the vacuum. Λ reads the universe as it is.

## 39. Why Λ cannot be steered physically — the full argument (`research/WHY_LAMBDA_IS_UNSTEERABLE.md`)

- **Theorems:** Lovelock (Λ is a constant of nature, no adjustable material), NEC +
  passivity (any field has w≥−1; the vacuum ground state is passive — 2nd law), QI
  (sustained negative energy needs boundaries, which are zero-sum).
- **Knobs quantified:** faking the low-z step needs δρ_Λ/ρ_Λ≈2.8% across the whole
  445-Mpc z<0.1 sphere (N4 bound: <0.9%/dex host mass); local void/backreaction
  (δH/H~2.5% → δμ≈0.054 mag) steers ONLY z≈0 — the anchor, not the dip (at z~1,
  vpec gives δμ~0.002 mag); lensing σ_μ≈0.02 mag at z~1 vs dip signal ≈0.10 mag
  coherent in 3 samples — 5× above noise, cannot be faked; Casimir at cosmic scale
  needs 1.2×10²⁷× the observable universe in plates; global Λ change is a boundary
  condition of the whole universe, excluded by z>0.1 + Lyα AP (χ²=0.97/2).
- **Self-protection:** every tracer (BAO ratios, 3 SN samples, DESI's own chains,
  CMB, Lyα AP) reads the same Λ and they agree (§29/§33/§34/§37) — no slack to steer.
- **The one crack:** if the z~1 dip is real, w<−1 there (NEC violation, N2 phantom
  field) — a local coupling becomes conceivable but is capped at 8×10⁻²⁸ W/m³ and
  bounded by the <1%/dex tilt. DR3 decides if the vacuum is a moving field; even
  then we read it, we can't turn it.

## 40. Conceptual physical steering — the two rudders (`session26.py rudder` → `rudder_results.json`, `research/STEERING_N14_N15.md`)

- **N14 Vacuum Rudder** (steer Λ_eff locally via the phantom sector): requires the z~1
  dip to be real (w<−1, NEC violation — DR3 gate). Posterior phantom fraction
  |1+w| = 0.27 (0.11–0.48); local energy-creation 8×10⁻²⁸ W/m³; max local tilt
  0.6 mag/e-fold if fully coupled — but any astrophysical tilt is bounded <0.01
  mag/dex (§19). Best case = a **lab-scale Λ_eff modulator**; it can never steer the
  cosmological reading (needs coherence over 10⁷⁶ m³).
- **N15 Optical Rudder** (steer the reading via an engineered line-of-sight medium —
  no new physics): faking δμ=0.03 mag at z=0.05 needs δt~3×10¹⁴ s → DM~2×10²⁸
  pc/cm³ (10²⁶× the whole IGM) or a refractive index −1.4% over 200 Mpc (~10¹⁷× IGM
  density); both are frequency-dependent (Δt∝1/ν²) → caught by multi-band SN + radio
  BAO + Lyα AP + CMB all reading the same geometry.
- **Door 3 (global Λ, unimodular boundary condition):** changing it changes the
  universe at all z — excluded by z>0.1 agreement and Lyα AP (χ²=0.97/2, §33).

**Reading:** best-case physical steering is the N14 phantom rudder, gated on DR3. If
the dip dies, all doors close in principle (Lovelock+NEC+passivity+QI). If it lives,
we get a local lever — but the universe's reading stays unfakeable due to
multi-tracer averaging.

## 41. Bar raised to 5σ + strongest honest dip from public data (`session28.py union` → `union_dip_results.json`)

Per request the discovery bar moved **3σ → 5σ** (Δχ²(k=2) ≥ 25.0; dip |w+1|/σ ≥ 5; 3–5σ = evidence tier only) and the passive DR3 watcher was removed (`dr3_plan.json`, `dr3_dial.py`).

**(a) Union deduplicated cert5** — DESY5 (1829, HD form) + Pantheon+ non-calibrators **deduplicated by CID** (375 duplicate CIDs and all 77 calibrators removed → 1249 added), with BAO(ALL)+Lyα DR2 AP+CMB+fσ8, 5 bins edges A:

w bins: −0.905±0.026 | −0.942±0.047 | **−1.331±0.166** | −0.690±0.271 | −1.122±0.927
→ **dip z∈[0.75,1.5]: w = −1.331±0.166, |w+1|/σ = 2.00** (χ²_map=2852, acc 0.19)

**(b) Inter-sample correlation scan** combining the three s22 per-sample dips over ρ∈[0,1]:

| ρ | combined w | σ | significance |
|---|---|---|---|
| 0.00 (independent) | −1.268 | 0.086 | **3.10σ** |
| 0.25 | −1.268 | 0.106 | 2.54σ |
| 0.50 | −1.268 | 0.122 | 2.20σ |
| 0.75 | −1.268 | 0.136 | 1.97σ |
| 1.00 (fully shared systematics) | −1.268 | 0.149 | 1.80σ |

**Reading against the raised bar:** even the most optimistic combination (ρ=0, assuming the three SN samples are statistically independent — they are not) reaches only **3.10σ**, and the honest union fit gives **2.00σ**. **Nothing reaches 5σ: the N14 phantom-rudder door stays CLOSED on current public data.** The dip remains the one live feature (anchor-free, all samples, DESI's own chains, bin-robust), but it is evidence-tier at best, not discovery.

## Caveats
- Diagonal (R, ω_b) CMB prior instead of full Planck / PR4 compressed triple.
- r_s from EH98 z_drag + numerical integral (~2–3% vs CAMB) — only enters via optional ℓ_A (off by default).
- Single adaptive-Metropolis chain per model, 24k steps (acc ≈ 0.2–0.3); MAP Δχ² quoted, not Bayesian evidence.
- Dovekie HD assumed absolute-magnitude-marginalized (χ²/dof ≈ 0.9 confirms consistency).
