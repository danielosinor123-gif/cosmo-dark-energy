# N14–N15 — conceptual PHYSICAL STEERING of the dark-energy reading (session 26)

Goal changed from *harnessing* to *steering the physical itself* — best-case scenario,
built on the three cracks identified in `WHY_LAMBDA_IS_UNSTEERABLE.md`:
(1) NEC/passivity (phantom sector), (2) QI/boundaries (Casimir), (3) global boundary
condition (unimodular Λ). Two rudders, each with its minimal assumption-to-break and
its falsifiable signature. Numbers: `results/rudder_results.json`.

---

## N14 — the VACUUM RUDDER (steer Λ_eff locally by coupling to the dark-sector field)

**Minimal assumption to break:** NEC — requires w < −1 somewhere, i.e. **the z~1 dip
must be real** (our cert5: w = −1.27±0.16 in z∈[0.75,1.5], §34/§36; DESI's own chains
w(1)=−1.1…−1.4, §37). DR3 is the gate.

**Mechanism.** If the dip is real, the dark sector is a two-field phantom/quintom
system (N2) with ρ+p < 0 in that epoch. A device with a portal coupling
(g φ F_μνF^μν, or nonminimal ξRφ²) can locally *excite* the phantom field. Because the
phantom has negative kinetic energy, exciting it **adds** energy density instead of
costing it — the field is a "phantom battery": the local effective Λ_eff = Λ + δV(φ)
can be pushed *up or down* by driving φ.

**Steering capability (posterior, §37/rudder_results.json):**
- phantom fraction |1+w| = **0.27 (0.11–0.48, 68%)** in the dip bin
- local energy-creation rate 3H|1+w|ρ_Λc² = **8×10⁻²⁸ W/m³** (5.5–10×10⁻²⁸)
- max local tilt imprint: δμ ≈ 0.6 mag per e-fold *if fully coupled* — but the data
  bound any astrophysical local tilt at **<0.01 mag/dex (<0.9%/dex, §19)**.

**Why it still can't steer the probe:** the cosmological reading is the average over
~10⁷⁶ m³; a rudder tilts only its own device volume. To steer the *cosmological* w(z)
you would need the phantom excitation coherent across a 445-Mpc sphere — the same
2.8%-of-Λ requirement as before, with the same EP/Lovelock obstruction at every point.
**Best case: a lab-scale Λ_eff modulator** — real local physics, invisible to the
cosmological probe (which is exactly why the probe is safe).

**Signature if it existed at astrophysical scale:** host-mass / position-dependent μ
tilt >0.01 mag/dex — **already excluded** (§19). DR3 + a host-mass re-test with the
dip assumed real is the falsification plan.

---

## N15 — the OPTICAL RUDDER (steer the READING, not the vacuum: engineered photon path)

**Minimal assumption to break:** *none* — photon propagation in a medium is standard
physics. You don't touch Λ; you put a designed medium along the line of sight so the
apparent distance changes. This is the honest "steer the reading physically" option.

**Quantified impossibility (`rudder_results.json`):**
- Faking δμ = 0.03 mag at z=0.05 (DL≈200 Mpc) needs a light-travel difference of
  δt ≈ 3×10¹⁴ s (~10 Myr). As a dispersive delay (Δt ∝ DM/ν²) that needs
  **DM ~ 2×10²⁸ pc/cm³ ≈ 10²⁶× the entire IGM column**; as a refractive shift
  (n−1 ≈ −1.4% over 200 Mpc) it needs electron density ~10¹⁷× the IGM, spread with
  exactly the right z-dependence over 445 Mpc.
- **Signature:** both mechanisms are frequency-dependent (Δt ∝ 1/ν², n(ν)). Radio
  BAO, optical SN, X-ray/GRB timing would all disagree — multi-band + multi-tracer
  consistency (BAO ratios, Lyα AP at z=2.33, CMB) reads the same geometry and catches
  it. Distance-duality (DL = (1+z)² DA) tests would also break.

**Verdict N15:** conceptually valid, physically absurd by 16–26 orders, and
self-defeating (dispersion signature) — the reading cannot be steered through the
light path without every other tracer noticing.

---

## The best-case scenario, honestly stated

- **Door 1 (N14, phantom rudder):** exists *only if* the z~1 dip is real phantom
  physics. Even then it steers only a local volume (≤8×10⁻²⁸ W/m³, tilt bounded
  <0.9%/dex by our own data). It can never steer the cosmological reading, because
  that requires coherence over 10⁷⁶ m³.
- **Door 2 (N15, optical rudder):** no new physics needed but needs 10¹⁶–10²⁶× the
  IGM along 445 Mpc with the right z-shape, and betrays itself by frequency
  dispersion across bands/tracers.
- **Door 3 (global Λ = boundary condition, unimodular):** changing it is changing
  the universe — visible at *all* redshifts at once; excluded by z>0.1 agreement +
  Lyα AP (χ²=0.97/2, §33).

**Conclusion:** the best-case physical steering device is the N14 phantom rudder,
gated on DR3 confirming the z~1 dip. If the dip dies in DR3, all three doors close
and Λ is confirmed unsteerable in principle (Lovelock + NEC + passivity + QI), not
just in practice. If the dip lives, we get a lab-scale Λ_eff modulator — a genuine
new physical lever — while the cosmological probe stays protected by multi-tracer
averaging. Either way the probe's reading of the *universe* cannot be faked; only a
local bubble can be steered, and only if phantom physics is real.
