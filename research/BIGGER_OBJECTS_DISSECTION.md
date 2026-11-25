# "Move bigger objects" — dissected: every scaling road, and where it breaks
(session 32, `experiments/session32.py scale` → `results/bigger_objects.json`)

Novelty: yes as a design challenge. Physics: no — bigger mass makes every term worse,
not better. Here is the dissection.

## The master scaling law
For Casimir/Λ pressure at gap d, P = K/d⁴, so to hold mass m:

**A = m·g·d⁴ / K** grows as d⁴, while the gap tolerance for ±10% force is **δd = d/40**
(shrinks with d). Smaller gap = more force but far tighter flatness. For **1 kg**:

| gap | pressure | area for 1 kg | gap tolerance | plate mass (1 μm, 2 plates) | collapse energy | 1-mm-cell actuators |
|---|---|---|---|---|---|---|
| 10 nm | 1.3×10⁵ Pa | 0.75 cm² | **0.25 nm** | 0.4 mg | 3×10⁻⁸ J | 75 |
| 100 nm | 13 Pa | 0.75 m² | **2.5 nm** | 4 g | 3×10⁻⁷ J | 7.5×10⁵ |
| 1 μm | 1.3×10⁻³ Pa | 7.5×10³ m² | 25 nm | 40 kg | 3×10⁻⁶ J | 7.5×10⁹ |
| 10 μm | 1.3×10⁻⁷ Pa | 7.5×10⁷ m² | 250 nm | 4×10⁵ kg | 3×10⁻⁵ J | 7.5×10¹³ |
| **40 μm = d\*** | **5.25×10⁻¹⁰ Pa = ρ_Λ** | **1.9×10¹⁰ m² (~19,000 km²)** | 1 μm | 1×10⁸ kg | 1.3×10⁻⁴ J | 1.9×10¹⁶ |
| 100 μm | 1.3×10⁻¹¹ Pa | 7.5×10¹¹ m² | 2.5 μm | 4×10⁹ kg | 3×10⁻⁴ J | 7.5×10¹⁷ |

## Where each road breaks
1. **Area scaling**: 1 kg needs continent-scale plates already at 40 μm, and 7500 m²
   even at 1 μm. Plate mass and flatness become the object itself.
2. **Gap tolerance**: at 100 nm you must hold **2.5 nm** flatness over the whole area —
   below real surface roughness; at 10 nm it’s 0.25 nm, i.e. atomic. Impossible without
   active leveling.
3. **Actuator explosion**: controlling A with 1 mm cells needs N = A/10⁻⁶ — 750,000 cells
   for 1 kg at 100 nm, 7.5×10⁹ at 1 μm, 10¹⁶+ at d\*. Wiring, power, heat, and control
   complexity kill it.
4. **Instability/pull-in**: attractive Casimir is unstable (gap opens → force drops).
   Repulsive Lifshitz helps but doesn’t change the d⁻⁴ scaling or tolerance problem.
5. **Not an energy source**: collapse energy for 1 kg is 10⁻⁸–10⁻⁴ J — you get back
   almost nothing; the plates always cost more (R ≥ 10⁹–10¹⁴ at d\*).
6. **Dark energy cannot lift**: Λ acceleration reaches g only at
   r = 3.0×10³⁶ m ≈ **2×10¹⁰ × the observable radius** — never at any accessible place.
7. **Vacuum balloon**: canceling Λ inside a bubble gives density contrast ρ_Λ = 5.8×10⁻²⁷
   kg/m³; buoyancy lifts 1 kg only with V = 1.7×10²⁵ m³ (~0.8 million Moon volumes),
   and the shell costs R ≥ 10⁹–10¹⁴.
8. **Horizon piston**: enclosing the observable universe would theoretically give
   3.7×10⁵² W of expansion work, but it requires a piston shell at the horizon — not
   physical, not causal. Per m³ it is still 3.4×10⁻²⁷ W.
9. **Fractal/nested plates & metamaterials**: raise area per volume but are bounded by
   the same gap-tolerance, actuator-count, QI, and material-strength limits.

## What actually moves bigger objects
Ordinary forces: **magnetic/electromagnetic levitation** (maglev trains), diamagnetic
levitation (frogs), electrostatics, acoustics, aerodynamics. Those move kg–ton objects
today — but they are not dark energy, and adding Λ to them changes nothing measurable.

## Conclusion
Bigger-object movement via dark energy/redistribution is **impossible by scaling**: force
∝ d⁻⁴ while required flatness ∝ d, area ∝ m·d⁴, actuators ∝ area, material cost
≥10⁹–10¹⁴, and Λ’s own gravity is ~10³⁶ too weak. The only surviving use remains
**G1: measure Λ at 40 μm**; S2 stays a **milligram-scale micro-actuator/levitator**.
Practical large objects are moved by electromagnetism etc., not by dark energy.
