# cosmo-dark-energy — independent re-analysis of the 2025–26 evolving-dark-energy hint + harnessing bounds

Independent pipeline (numpy/scipy only) on public data: DESI DR2 BAO (13-pt, full cov),
DES-Dovekie / Pantheon+ / DESY5 supernovae, Planck-2018 compressed CMB prior (R, ω_b),
public 24-pt fσ8 growth compilation, plus DESI's own published Y3 chains.

## Headline results (`results/RESULTS.md` §1–§32, `research/VERDICT.md`)

1. The 2025–26 evolving-DE hint (w0wa Δχ²=−6.3, BF 9–12) lives **entirely in the z<0.1
   SN tail** and is **equivalent to a +0.03 mag sharp step there** (capstone with all
   systematics free: Δχ²=−0.26, wa=−0.07±0.13) — §8, §16, §21, §25.
2. The step **is the Cepheid calibrators**: 77 calibrator SNe sit +0.20±0.03 mag (6.7σ)
   off BAO+CMB geometry → implied H₀=74.3±0.2; the other 664 low-z SNe agree at
   −0.012±0.007 — §27, v_pec-robust §28. H₀↔wₐ map: wₐ≈−0.17+0.32(H₀−70) — §15.
3. SN-independent geometry does not carry the signal: DESI's own BAO+2pt = BAO-only;
   BAO+SN+2pt = BAO+SN; Lyα alone Δχ²=+0.09 — §20, §29.
4. One anchor-free feature survives: a **phantom dip w=−1.26±0.16 at z∈[0.75,1.5]**
   (5-bin w(z)+growth), identical with/without calibrators (−1.264 vs −1.257), best
   fit by a quintom crossing (Δχ²=−7.4) — §24, §31, §32. ~1.6σ; DR3 will decide.
5. New physics judged by data: LDPT killed (BF 0.11), quintom allowed-not-required,
   RVM null (ν=0.008±0.005) — §17, §22, §32.
6. **Harnessing caps** (§23, §26, `research/HARNESS_N6_N10.md`, `HARNESS_LITERATURE.md`):
   one-shot ≤1.6×10⁻¹⁰ J/m³, continuous ≤H₀ρ_Λc²≈1.2×10⁻²⁷ W/m³ (1 m³ powers humanity
   for ~10⁴⁶ yr — dilution, ~26 orders below solar flux), static tilt <0.01 mag/dex;
   no angular/dust/alternation signatures. Only the EW false vacuum is denser
   (7.6×10⁴⁶ J/m³), bounded by Higgs stability. **Dark energy is a probe, not a resource.**
7. DR3 pre-registered plan + forecast: w0wa↔δM degeneracy breaks at 3σ only at ~3× BAO
   precision (`results/dr3_plan.json`, §30).

## Layout

- `src/cosmo.py` — background, BAO/SN/CMB/fσ8 likelihoods, loaders, MCMC
- `src/models.py` — LDPT, quintom, generic BackgroundFn
- `experiments/session1..20.py`, `desi_2pt.py`, `desi_minima.py` — all analyses
- `results/` — RESULTS.md (§1–§32 + caveats), all JSON/npy outputs
- `research/` — LANDSCAPE.md, NEW_PHYSICS.md (N1–N5), HARNESS_N6_N10.md (N6–N10),
  HARNESS_LITERATURE.md, PAPER_DRAFT.md, VERDICT.md
- `data/` — public datasets (DESI DR2 BAO, SN samples, DESI Y3 chains, fσ8)

## Reproduce

```
python experiments/session3.py evidence     # BF w0wa vs LCDM
python experiments/session9.py capstone     # all-systematics joint fit
python experiments/session17.py cert5noc    # dip without calibrators
python experiments/session14.py lya         # Ly-alpha-only test
python experiments/desi_2pt.py              # DESI's own chains, 2pt invariance
```
