# The redistribution channel, derived to its end — a theorem, not a toy fit
(session 29, `experiments/session29.py redist` → `results/redistribution_theorem.json`)

We did **not** re-fit someone else's paper here. We took the redistribution idea
literally — "move Λ out of a volume into somewhere else" — and derived what a
redistributing boundary/curvature must look like, then computed whether it can exist.

## 1. The universal micro-scale of Λ-redistribution

Any local stand-in for Λ made from boundaries or curvature has an energy density set
by the smallest length scale ℓ of the device:

- parallel plates (Casimir): ρ_C = −K/d⁴, K = π²ħc/240 = 1.300×10⁻²⁷ J·m
- generic curvature (trace anomaly / vacuum polarization): ρ_curv ~ −ħc/ℓ⁴

Matching the cosmological density ρ_Λc² = 5.253×10⁻¹⁰ J/m³ gives a **universal length**:

- plates: d* = (K/ρ_Λc²)¼ = **39.7 μm**
- curvature: ℓ* = (ħc/ρ_Λc²)¼ = **88.1 μm**

**Discovery of the derivation:** Λ has exactly one natural material scale, **40–90 μm**.
Above that scale any boundary/curvature vacuum term is diluted below ρ_Λ; below it, it
overshoots. So redistribution can only ever be a *microscopic* phenomenon — the
cosmological volume averages it away.

## 2. Killer A — passivity (already known, restated as a law)

A Casimir cavity that cancels Λ in its gap stores the removed energy in the plates'
stress. Re-opening the gap returns it. Per cycle:

W_net = ∮ p dV = 0 .

Only expansion does net work on the cavity: dW = ρ_Λc² dV → ceiling 3Hρ_Λc²
= 3.4×10⁻²⁷ W/m³ (§37). Redistribution without expansion is zero-sum.

## 3. Killer B — the material-cost bound (new, quantitative)

Suppose plates of area A, gap d, thickness t try to cancel Λ over volume V. The
canceled energy is

E_C = K A / d³ = ρ_Λc² V  →  A = ρ_Λc² V d³ / K .

But the plates themselves carry positive rest energy E_m = ρ_m c² (2At). The ratio is

R = E_m / E_C = 2 ρ_m c² t d³ / K .

Scan (`redistribution_theorem.json`, aluminum, one-atom layer t = 0.3 nm):

| gap d | Casimir pressure | fill factor to cancel Λ in 1 m³ | plate area per m³ | R = plate mass-energy / canceled Λ |
|---|---|---|---|---|
| 0.3 nm | 1.6×10¹¹ Pa | 3.3×10⁻²¹ | 1.1×10⁻¹¹ m² | **3.0×10⁹** (absolute min, atomic contact) |
| 1 nm | 1.3×10⁹ Pa | 4.0×10⁻¹⁹ | 4.0×10⁻¹⁰ m² | 1.1×10¹¹ |
| 10 nm | 1.3×10⁵ Pa (~1 atm) | 4.0×10⁻¹⁵ | 4.0×10⁻⁷ m² | **1.1×10¹⁴** |
| 1 μm | 1.3×10⁻³ Pa | 4.0×10⁻⁷ | 0.40 m² | 1.1×10²⁰ |
| **39.7 μm** | **5.25×10⁻¹⁰ Pa = ρ_Λ** | 1 | 2.5×10⁴ m² | 7.0×10²⁴ |

So even in the most absurd limit — one-atom-thick plates at atomic contact, where the
Casimir formula is already being pushed — the plates' rest energy is **≥3×10⁹ times**
the Λ energy they cancel. For any mechanically stable gap (≳10 nm, pressure ≲1 atm)
the penalty is **≥10¹⁴**. **No material can pay for the redistribution it performs.**
This is independent of Killer A: even a massless plate would still be zero-sum per
cycle; ordinary matter makes it worse by ≥9 orders.

## 4. The curvature / black-hole channel

Set curvature vacuum polarization ~ħc/R_s⁴ equal to ρ_Λc². A horizon with
R_s = ℓ* = 88 μm has mass M = R_s c²/2G = **5.93×10²² kg ≈ 0.81 Moon**. Its Hawking
power is P = ħc⁶/(15360πG²M²) ≈ **1.0×10⁻¹³ W**, lifetime ~5×10⁴⁴ yr. A black hole
is the only object whose *local* curvature can carry Λ-scale vacuum energy — and it
radiates 10⁻¹³ W. It redistributes nothing.

## 5. Dynamical Casimir (moving boundaries make real photons — still not a source)

A moving mirror can turn virtual pairs into real photons (dynamical Casimir effect).
Energy conservation is exact: the actuator supplies at least the photon energy; the
vacuum only provides the virtual intermediate state. Creating one 1 eV photon costs
≥1 eV of drive work. **Net extractable from the vacuum: 0.**

## 6. Theorem and falsifiable prediction

**Redistribution theorem.** Let Λ_eff(x) = Λ + δΛ(x) be produced by any local,
passive boundary/curvature device of feature scale ℓ and material density ρ_m. Then
either (i) ℓ ≳ 40–90 μm, in which case |δΛ| ≤ ρ_Λ only in a microscopic volume and
the device mass-energy exceeds the shifted energy by R ≥ 3×10⁹ (≥10¹⁴ for stable
gaps); or (ii) the only net-positive term is cosmic expansion, capped at
3Hρ_Λc² ≈ 3.4×10⁻²⁷ W/m³. Hence **macroscopic redistribution of Λ is impossible, and
the maximum any redistribution device can yield is the expansion term.**

**Falsifiable prediction:** a cavity whose gap distribution has volume-averaged
Casimir density −ρ_Λ over a macroscopic volume V would show a Casimir pressure of
5.25×10⁻¹⁰ Pa at the matching gap ~40 μm — about 10⁹ below present force sensitivity —
and would require plate mass energy ≥10⁹ (realistically ≥10¹⁴) times the canceled
energy. If a future experiment measures a stable macroscopic region with
Λ_eff measurably different from Λ, this theorem is wrong and DR3-grade data should
show a sky dipole in μ (§26 already null at Δχ²≈−0.05…−0.27).

**Bottom line:** the redistribution idea is not a toy — it is a real channel with a
universal 40–90 μm scale, and it dies twice: zero-sum per cycle (passivity) and
material cost ≥10⁹–10¹⁴× (rest energy). The only surviving output is expansion work
itself, capped at 10⁻²⁷ W/m³. Dark energy can be *moved* only inside microscopic
cavities at enormous material cost; it cannot be moved across the sky.
