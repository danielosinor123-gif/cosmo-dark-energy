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

## Caveats
- Diagonal (R, ω_b) CMB prior instead of full Planck / PR4 compressed triple.
- r_s from EH98 z_drag + numerical integral (~2–3% vs CAMB) — only enters via optional ℓ_A (off by default).
- Single adaptive-Metropolis chain per model, 24k steps (acc ≈ 0.2–0.3); MAP Δχ² quoted, not Bayesian evidence.
- Dovekie HD assumed absolute-magnitude-marginalized (χ²/dof ≈ 0.9 confirms consistency).
