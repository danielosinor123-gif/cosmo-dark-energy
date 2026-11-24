# Cosmology research landscape — dark sector (as of Mon Nov 23, 2026)

Live sources checked today: desi.lbl.gov (2026 posts), data.desi.lbl.gov/doc/papers (DR2 list),
arXiv API (Nov 2026 submissions). Purpose: establish the exact 2026 state so our multi-session
work lands on a real, open, falsifiable question that needs NO new instruments — only public
data + analysis.

## 1. The two problems

- **Dark matter (~27% of the universe).** No particle detected despite decades of direct
  (LZ/XENONnT/PandaX/DarkSide), indirect, and collider searches. Nov 2026 status: LZ WS2024 run
  is null (arXiv 2603.26877); DAMA/LIBRA annual modulation is excluded against LZ/PICO/CDMS in
  an extended EFT analysis (Heikinheimo et al., arXiv 2611.04890, Nov 3 2026); WIMP parameter
  space is being pushed into the neutrino fog (Arcadi et al. 2606.23457); SuperCDMS SNOLAB
  starts 2026 in the 0.5–10 GeV range; non-WIMP candidates (UFO/relativistic decoupling,
  arXiv 2605.03014; axions) are the growth area. On the astrophysical side, DESI DR2 weighed the
  Milky Way at M200c = 1.2e12 Msun, 11% uncertainty, with a key lesson: method-vs-method
  systematics dominate over statistical errors unless validated on common mocks
  (arXiv 2611.05518, Nov 4 2026).
- **Dark energy (~68%) + the cosmological constant problem.** Observed: late-time acceleration.
  The vacuum-energy prediction from QFT is off by ~10^120 — still no resolution; it is a
  theory-side problem (needs new physics ideas, not more laptop data). The *observational*
  question — is w(z) = −1 or evolving? — is live and data-driven (below).

## 2. The live 2026 question: is dark energy evolving?

- DESI DR2 ("Y3": first 3 years) BAO + CMB + Type Ia SN (Pantheon+, Union3, DES-Dovekie)
  favor w0waCDM (CPL) over ΛCDM at **~3.1–3.2σ** (Δχ² ≈ 12.8), and this survives
  covariance-rescaling robustness tests (Kolonia, Li, Zhang, Afshordi, arXiv 2611.13829,
  Nov 10 2026: 3.15σ → 3.07σ with free covariance rescaling).
- **Counter-signal:** DESI DR2 Lyα-forest *full-shape* Alcock–Paczynski analysis
  (arXiv 2607.27410, Jul 30 2026) shifts back toward ΛCDM — "the current hints of evolving
  dark energy may fade away, or a more complex model may be needed."
- Nov 2026 theory responses fitting the public DESI DR2 chains: Fading Dark Sector
  (Obied–Vafa–Wu, arXiv 2611.17422, 4.1σ over ΛCDM with full-shape), VCDM one-parameter
  histories F42/Ag/Qui (De Felice et al. 2611.16407), HoloEoS Hubble-running (2611.17374),
  thermodynamic TPDE (2611.23340), unified BDE-BDM dark sector (2611.13567), coupled-DE
  potential-independent lower bound β_min ≈ 0.047–0.058 (Sapone & Arratia 2611.12792),
  fast-scrambling cosmology (Afshordi et al. 2611.13100), PINN reconstructions finding phantom
  crossing at z ≈ 0.41–0.53 (Paliathanasis 2611.19117), and joint background+gravity
  reconstructions with an Ω_X peak at z ≈ 0.6–0.7 (Gallavotti & Silvestri 2611.18557/58).

**Public data available to us (no instrument needed):**
- DESI DR2 BAO cosmology chains + likelihood inputs:
  https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/ (DR2 spectra not yet public,
  but BAO chains/posteriors are).
- DESI DR1 full public data release: https://data.desi.lbl.gov/doc/releases/dr1/
- Planck 2018 chains + compressed distance priors (public).
- Pantheon+ (github.com/dscolnic/PantheonPlus), Union3, DES-SN5YR/Dovekie — public SN
  likelihoods with covariances.
- Cosmic chronometers and fσ8/RSD compilations — public.

## 3. The gap we can actually fill from a laptop

Everyone is fitting models to the same public chains. What is NOT yet done (checked against
the Nov 2026 arXiv wave above):

1. **Probe-by-probe consistency "certificates" for w(z):** reconstruct w(z) (or distance
   history) independently from BAO-only, SN-only, Lyα-only, cosmic-chronometer-only, RSD-only
   and quantify *where they disagree* — disagreement localizes systematics, which is exactly
   what decides whether the 3σ hint is real. (Analog of ARIADNE's verification inheritance:
   each probe certifies the others.)
2. **Independent re-derivation of the 3.1σ with our own pipeline** (our own MCMC on the
   public BAO covariances + SN + Planck compressed priors), plus leave-one-probe-out influence
   and a proper Bayesian evidence comparison — reproducibility + new robustness angles.
3. **Novel falsifiable theory tests on public chains** (like Sapone's β_min, but for another
   mechanism class), and/or a small differentiable emulator of the distance integrals to make
   nested sampling cheap.

## 4. Multi-session plan (proposal)

- S1 (today): landscape + workspace + environment (numpy 2.5.1, scipy 1.18.0; pip available).
- S2: download public data (DESI DR2 BAO chains/covariances, Planck compressed priors,
  Pantheon+/Union3/DES-Dovekie, CC, RSD); build numpy/scipy likelihood pipeline; our own MCMC.
- S3: reproduce ΛCDM vs w0waCDM on public data; quantify the preference; leave-one-out.
- S4: probe-by-probe w(z)/distance reconstruction agreement certificates (the novel test).
- S5: evidence-based model comparison incl. 1-parameter physics models (FDS-like, VCDM-F42).
- S6: write results.md + paper-style draft; push to GitHub.

Honesty note: we cannot "solve" dark energy from a laptop — the 10^120 problem needs new
physics and new instruments. What we CAN do is produce a rigorous, novel, falsifiable analysis
of the live 3σ evolving-DE hint using only public data — and kill it or confirm it ruthlessly.

---

## Session log

### S1 (Nov 23, 2026) — landscape + workspace
- Workspace `cosmo/` (research/, src/, data/, experiments/, results/).
- `src/cosmo.py`: flat w0waCDM background (CPL), cached-grid distances, numerical sound
  horizon, loaders for DESI DR2 BAO / DES-Dovekie SN / Pantheon+, Planck-2018 compressed
  CMB prior, adaptive Metropolis sampler. numpy/scipy only.
- Data cloned: `CobayaSampler/bao_data` (DESI DR2 BAO, `desi_bao_dr2/`, 13-dim ALL vector:
  BGS DV/rd + 6×(DM/rd, DH/rd)), `CobayaSampler/sn_data` (Pantheon+, Union3, DESY5,
  DES-Dovekie — Dovekie ships precomputed *inverse* total covariance as packed triangle).
- CMB prior decision: Planck 2018 TT,TE,EE+lowE compressed (R=1.7502±0.0046, ωb=0.02237±0.00015)
  used as diagonal Gaussian on (R, ωb) only — ℓA needs rs and our EH98-z_d+numerical-integral
  rs is ~2-3% off CAMB (would be ~100σ in ℓA); documented as a systematic, robustness run with
  ℓA optional. DESI's own "CMB-compressed" = PR4 (θ*, ωb, ωc) — data file is NERSC-only, so we
  use PR3 (R, ωb) instead.
- Validation reference downloaded: DESI DR2 published chains
  (`data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/cobaya/`) for base and base_w_wa with
  BAO+CMB-compressed: w = −0.44±0.21, wa = −1.68±0.61, P(w>−1)=0.996, P(wa<0)=0.997
  (base: Ωm=0.3008±0.0043, H0=68.34±0.34, rdrag=147.75±0.20). See
  `experiments/desi_reference.py`.
- Bugs fixed this session: missing factor of c in rs integral; H prefactor (100, not H0, when
  using ωh² densities); Dovekie packed-triangle covariance unpack; 1×1/diagonal per-tracer BAO
  covariances.

### S2 (Nov 23, 2026) — own pipeline running
- `experiments/fit_models.py`: ΛCDM / wCDM / w0waCDM × combos {BAO+SN+CMB, BAO+SN, BAO+CMB,
  BAO-only, SN+CMB} + per-tracer leave-one-out; our own MCMC; outputs `results/*.json` + chains.
- Physics sanity at Planck-fiducial: rs=150.5 Mpc (CAMB 147.1), BAO χ²=28.8/13, SN χ²=1936/1820,
  CMB(R,ωb) χ²=2.9; ~5 ms per likelihood eval.
- Bug fixed: Ω_b prior range (0.010,0.040) capped the ω_b prior → pile-up at boundary and a
  spurious wa>0 solution; widened to (0.005,0.10). After fix: BAO+SN+CMB w0wa = w0 −0.87±0.07,
  wa −0.35±0.30 (1.73σ vs ΛCDM), matching DESI's BAO+DES-SN5YR chain within ~1σ.

### S2 results (see results/RESULTS.md, repo danielosinor123-gif/cosmo-dark-energy)
- w0wa preferred over ΛCDM in every combo (0.7–1.7σ), always w0>−1, wa<0; per-tracer LOO:
  1.0–1.5σ each, w0≈−0.87±0.11 — no single tracer drives it.
- Binned-w probe certificates (probe_consistency.py): joint w(z) dips to −1.31±0.21 in
  z∈[0.75,1.5] (BAO-driven); max pairwise tension vs joint < 1σ — internally consistent.
- Gap to DESI's 3.1σ is mostly CMB-prior strength (we use diagonal R,ω_b; DESI uses PR4
  θ*,ω_b,ω_c compressed triple). SN sample choice also matters (Pantheon+ weakens wa).

### S3 (Nov 23, 2026) — stress tests + evidence + growth probe
- `experiments/session3.py`: SN-sample systematics (Dovekie/Pantheon+/DESY5 — fixed the
  packed-covariance loaders: first entry = dimension, then n×n row-major), z-cut systematics,
  multi-seed stability, BAO covariance-inflation stress (profile over chains), and Bayes
  factors via Savage-Dickey (KDE at w0=-1,wa=0) + Laplace at refined MAP.
- Findings so far (results/RESULTS.md §5-6): BAO-only hint dies under 2× cov inflation, but
  BAO+SN+CMB survives 5× at ~1.2σ; **BF(w0wa/ΛCDM) ≈ 9–12 (strong, Jeffreys) for BAO+SN+CMB**,
  while wCDM is not preferred over ΛCDM → the data want *evolution* (wa≠0), not constant w.
- Growth probe DONE: public 24-pt fσ8 benchmark (AparajeetS/fsigma8-growth-consistency =
  Sagredo+2018 Table I + extensions, WiggleZ/eBOSS covariance blocks) in `data/fs8/`; exact
  growing-mode integral (Heath) for D(a), f(z) in `src/cosmo.py` (validated vs ODE +
  Ωm^0.55; ~20 ms/call). Result: fσ8 kills the *evolution* (wa→0, w0wa=wcdm) and BAO+fσ8
  shows no w0wa preference.
- SN systematics DONE: Dovekie 1.73σ / Pantheon+ 3.60σ / DESY5 2.91σ (w0wa vs ΛCDM,
  BAO+CMB) — same sign everywhere, significance swings with sample. Loader fixes:
  Pantheon+ = m_b_corr + zHD (mB gives χ²/dof 2.6); DESY5 total cov = covsys +
  diag(MUERR_FINAL²) (covsys alone not PD).
- z-cut DONE: signal entirely in z<0.1 SNe (z>0.1 cut → Δχ² −0.01, 0σ).
- Seeds DONE: Δχ² −6.34..−6.48 across 4 seeds — converged.
- **Session-3 verdict (results/RESULTS.md §9): hint is real in data (BF≈10, probe-consistent)
  but not robust — lives in low-z SNe, erased by growth, sample-dependent → most
  parsimonious explanation is a low-z SN calibration systematic.** Falsifiable prediction:
  DR3/full-shape + non-Cepheid low-z anchor should collapse the hint.
- S4 started: (a) DESI's own iminuit minima parsed (`experiments/desi_minima.py`,
  RESULTS.md §10) — **their** BAO+DESY5SN+CMB zmin=0.1 cut drops w0wa from 2.09σ to 1.06σ,
  confirming our z<0.1 smoking gun in DESI's own pipeline; (b) anchor test running
  (H0 priors Cepheid 73.04 / TRGB 70.4 / JAGB 67.8 / Planck 67.36 × z-cut 0/0.1,
  `session4.py anchor`); (c) growth-inclusive binned-w certificates DONE
  (`session4.py cert`, RESULTS.md §11): all probe combos incl. fσ8 recover the same binned
  w(z) shape (w≈−0.93 today, dip ≈−1.2 at z 0.75–1.5) with max tension 0.51σ — the
  non-constant w(z) shape is probe-robust; only its parametric significance is fragile.
- S4 anchor test DONE (RESULTS.md §12): Cepheid H0=73 anchor collapses the full-sample
  w0wa hint to 0.42 sigma; TRGB/JAGB leave it at 1.7-1.8 sigma; Planck H0=67.36 inflates it
  to 4.36 sigma (4.16 sigma even with z>0.1 SNe cut) — the hint is a composite of the
  low-z-SN-vs-Cepheid channel and the H0-tension channel (dH0~5 km/s/Mpc ~ wa~-0.7).
- S5 results so far: (b) xres DONE (RESULTS.md §13) — all 3 SN samples show the SAME
  coherent low-z excess vs BAO+CMB geometry (d=+0.033..0.048 mag, 3-4 sigma each) and
  agree with each other (<1 sigma pairwise) => the hint is ONE common low-z offset
  (shared Cepheid anchor), delta_mu~0.04 mag ~ delta_H0~1.1 km/s/Mpc; (d) bfbin DONE
  (§14) — nonparametric BF(3-bin w vs LCDM) = 17-135, decisive in most combos.
  (c) h0map DONE (§15): at fixed H0, LCDM chi2(H0) has a sharp minimum at H0~70.5 (H0
  tension measured directly, dchi2~440 to H0=66) while w0wa is flat in H0 (1432-1456)
  with wa ~ -0.17 + 0.32*(H0-70) — evolving DE and an H0 offset are the same curve.
- S5 (a) step DONE (§16): a dM=+0.02..0.03 mag step at z<0.1 absorbs the WHOLE w0wa
  preference (1.72 sigma -> 0.00 sigma at z_step=0.1) — same size as the §13 common
  low-z offset. THE quotable number of the project.
- S6 NEW PHYSICS (research/NEW_PHYSICS.md, §17): N1 LDPT (late dark-sector phase
  transition) — absorbs dchi2~105 at fixed H0=73 (h0map_s6_ldpt) but z_c->0.2 edge,
  degenerate with the low-z step; killed by evidence on the full sample (BF 0.11).
  N2 quintom phantom-dip — best raw chi2 (-8.06) but BF 2.2 < w0wa 8.0: allowed, not
  required. No model meets the discovery bar (BF>10 AND robust to z>0.1/anchor/dM).
- S7 DONE (§18): h0tension — in public compressed data, dNeff/EDE absorb NOTHING of the
  H0 tension (penalty at H0=73.04 = 219.7 for lcdm/neff/ede alike); only late-time freedom
  (w0wa flat in H0; LDPT dchi2~105) does. CMB peak-side constraint not in our prior — caveat.
- S8 DONE (§19): N4 environment test — Pantheon+ host-mass split of low-z residuals
  0.08 sigma, slope <1.6 sigma -> vacuum NOT locally tiltable at galaxy scales
  (<0.01 mag/dex); the low-z offset is global calibration, not environment.
- N5 "vacuum battery" (NEW_PHYSICS.md): energy budget computed; capped by same f<~0.3,
  z_c<~0.2 bounds as N1 (BF 0.11) — releasing vacuum energy is constrained by the same
  step that mimics calibration.
- S9 DONE: DESI own chains — BAO+2pt full-shape identical to BAO-only (wa +-0.97);
  SN distances carry all w precision (§20). CAPSTONE (§21): w0wa+{dM step + Cepheid
  anchor + fs8} vs LCDM+same: dchi2 = -0.26 (0.00 sigma), wa = -0.07+-0.13, dM=+0.034
  mag — the hint is fully absorbed by the low-z step.
- S10 running: RVM (running vacuum Lambda~rho_m) fits; cert5 (5-bin w(z)+growth) queued.
  Subagents: 2025-26 literature sweep for vacuum-energy extraction; N6-N10 harnessing
  proposals (NEW_PHYSICS.md to be extended).

### S4+ candidates (open)
- Lyα BAO vs Lyα full-shape tension (shape-vs-scale) — needs the DR2 Lyα full-shape likelihood.
- Emulator of distance integrals for cheap dense sampling; nested-sampling-grade evidence.
- Combine everything: full joint with growth + per-tracer certificates + SN systematics.
