# Results — independent re-analysis of the DESI DR2 evolving-dark-energy hint

Pipeline: `src/cosmo.py` (own numpy/scipy code; flat w0waCDM / wCDM / LCDM, binned-w option).
Data: DESI DR2 BAO (CobayaSampler/bao_data, `desi_bao_dr2`), DES-Dovekie SN Ia
(CobayaSampler/sn_data; 1820 SNe, precomputed inverse total covariance),
Planck 2018 TT,TE,EE+lowE compressed prior as diagonal Gaussian on (R, ω_b)
— ℓ_A omitted because our EH98-z_d numerical r_s is ~2–3% off CAMB (systematic noted).
Reference: DESI's own published DR2 chains (`data.desi.lbl.gov/public/papers/y3/bao-cosmo-params`).

## Main result (fit_models.py)
TBD — see main_results.json

## Per-tracer leave-one-out (fit_models.py tracers)
TBD — see tracer_results.json

## Probe-consistency certificates (probe_consistency.py)
TBD — see probe_results.json

## Reference values from DESI's own chains (w0waCDM)
| combo | w | wa | P(w>−1) | P(wa<0) |
|---|---|---|---|---|
| BAO-only | −0.49 ± 0.27 | −1.63 ± 0.97 | 0.954 | 0.927 |
| BAO + DES-SN5YR | −0.78 ± 0.07 | −0.72 ± 0.47 | 1.000 | 0.937 |
| BAO + Pantheon+ | −0.89 ± 0.06 | −0.18 ± 0.45 | 0.977 | 0.653 |
| SN-only (DESY5) | −0.82 ± 0.13 | −1.62 ± 0.93 | 0.909 | 0.938 |
| BAO + CMB-compressed | −0.44 ± 0.21 | −1.68 ± 0.61 | 0.996 | 0.997 |
(base ΛCDM, BAO+CMB-compressed: Ωm=0.3008±0.0043, H0=68.34±0.34, r_d=147.75±0.20)
