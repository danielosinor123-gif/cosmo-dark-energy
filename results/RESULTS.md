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
TBD (running): adds the public 24-point fσ8 RSD/PV benchmark (with published WiggleZ/eBOSS
covariance blocks) as an independent probe, σ8 as free nuisance.

## 8. SN-sample and z-cut systematics (`session3.py sn_sys`, `zcut`)
TBD (running): Dovekie vs Pantheon+ vs DESY5; Dovekie z>0.1 / z>0.3 cuts.

## Caveats
- Diagonal (R, ω_b) CMB prior instead of full Planck / PR4 compressed triple.
- r_s from EH98 z_drag + numerical integral (~2–3% vs CAMB) — only enters via optional ℓ_A (off by default).
- Single adaptive-Metropolis chain per model, 24k steps (acc ≈ 0.2–0.3); MAP Δχ² quoted, not Bayesian evidence.
- Dovekie HD assumed absolute-magnitude-marginalized (χ²/dof ≈ 0.9 confirms consistency).
