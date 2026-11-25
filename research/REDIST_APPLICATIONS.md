# What the 40–90 μm redistribution channel can actually do
(session 30, `experiments/session30.py apps` → `results/redist_apps.json`)

Given the theorem (`REDISTRIBUTION_THEOREM.md`), the honest applied envelope is small
but real. Numbers below are computed, not asserted.

## 1. It is a Λ-pressure gauge (the one clean application)

At the matching gap d = 39.7 μm the Casimir pressure equals the cosmological pressure:

P(d*) = ρ_Λc² = **5.253×10⁻¹⁰ Pa**.

Force = P·A:

| plate area | force |
|---|---|
| 1 mm² | 0.53 fN |
| 1 cm² | 53 fN |
| 1 m² | 0.53 nN |
| 1 km² | 0.53 mN |

To resolve a 1 fN threshold you need only **1.9 mm²** of plates; 1 pN needs **19 cm²**;
1 nN needs **1.9 m²**. So the Λ-scale pressure is *in principle* within reach of
precision force metrology — the catch is that at 40 μm the absolute Casimir signal is
buried under surface forces, so what you are really doing is calibrating the matching
point P(d)=K/d⁴, i.e. using a Casimir cavity as a direct gauge of ρ_Λ.

## 2. It can store and return negative energy — like a spring with negative efficiency

Collapsing a gap from d_eq down to d yields work per plate area
W/A = (K/3)(1/d³ − 1/d_eq³):

- to 1 μm: 4.3×10⁻¹⁰ J/m²
- to 100 nm: 4.3×10⁻⁷ J/m²
- to **10 nm: 4.33×10⁻⁴ J/m²**

Getting **1 J** this way needs **2307 m²** of 10 nm gaps — and the plate material cost
ratio there is R = 1.12×10¹⁴, so the net is **−1.12×10¹⁴ J**. It is a vacuum spring:
real, reversible, and strictly loss-making once matter is counted (and zero-sum even
before that, by passivity).

## 3. The only net-positive output is expansion — and it is absurdly small

Power density 3Hρ_Λc² = **3.44×10⁻²⁷ W/m³**. To get **1 W** you need a volume of
**2.9×10²⁶ m³** — a cube of side **6.6×10⁸ m** (~660,000 km), about **1.3×10⁷ Moon
volumes**. Energy delivered per m³ over a Hubble time ≈ ρ_Λc² ≈ **5×10⁻¹⁰ J**.

## 4. No propulsion, no local gravity lever, no macro shielding

Canceling Λ over a radius r = d*/2 = 19.8 μm changes the Λ-acceleration by
δg = (8πGρ_Λ/3c²) r ≈ **6.5×10⁻⁴¹ m/s²** — nothing. A *homogeneous* Λ-free region
would have local expansion H_local/H₀ = √Ω_m = **0.56** (44% slower), but making that
homogeneous requires fill factor 1 at 40 μm (2.5×10⁴ m² of plates per m³) with material
cost R = 7×10²⁴; sparse nano-gap averaging is not homogeneous and does not change
local expansion. So no anti-expansion bubble, no levitation, no warp.

## Bottom line — what redistribution can do

1. **Measure Λ directly** as a Casimir pressure at the 40 μm matching point (cm²–m²
   plates at fN–pN sensitivity; calibration, not a new force).
2. **Act as a reversible nanoscale negative-energy spring** (vacuum battery) — zero-sum
   by passivity, and ≥10¹⁴× loss-making with real plates.
3. **Harvest expansion work** only: 3.4×10⁻²⁷ W/m³ (1 W ≈ 13 million Moon volumes).

It cannot power anything, push anything, shield anything macroscopic, or change the
sky. Redistribution exists only as a microscopic dial and a sensor — which is exactly
why the cosmological probe stays readable.
