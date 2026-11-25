# Universal non-contact "invisible hand" — what in every solid can be grabbed, and the computed payload frontier
(session 33, `experiments/session33.py universal` → `results/universal_levitation.json`)

Assistive-only concept: a mouth-retainer/tongue or non-invasive BCI sets the **goal**; a
field system moves **objects** for disabled users. No implants, no harm, no
non-consensual use; see `ASSISTIVE_SAFETY.md`. This is **not** dark-energy levitation.

## The answer to "magnets need iron": grab charge, not magnetism
Every solid contains **electrons and nuclei**, so every solid has **permittivity** and can
be polarized or inductively charged. Electric fields therefore act on **all** materials:

- conductors → free electrons → **induction**: always attracted to a charged electrode
  (works on metal, water, tissue).
- dielectrics → bound electrons → **DEP**: Clausius–Mossotti CM=(ε_p−ε_m)/(ε_p+2ε_m).
  In **air** (ε_m=1) every solid has ε_p>1 → CM>0 → **attraction** to the electrode.
  Repulsion (CM<0) is obtained in a **higher-ε medium** (e.g. ε_p=2 object in ε_m=3
  fluid: CM=−0.125) or by a push field (acoustic/ion) from the other side.

**Dark matter / dark energy are not handles.** Dark matter has no known non-gravitational
coupling and cannot be localized; dark energy is uniform and repulsive (Λ pressure
5.25×10⁻¹⁰ Pa → load **5×10⁻¹¹ kg/m²**; Λ acceleration reaches g only at
r=3×10³⁶ m ≈ 2×10¹⁰× the observable radius). Using them is closed by
`REDISTRIBUTION_THEOREM.md` and `BIGGER_OBJECTS_DISSECTION.md`.

## Computed force frontier (`universal_levitation.json`)

| field | limit used | pressure | load | what it means |
|---|---|---|---|---|
| gravity (overhead mass at r) | a=g | — | need M=g r²/G: **3.7×10⁸ kg at 5 cm, 1.5×10⁹ kg at 10 cm** | impossible portable |
| electrostatic in **air** | E_bd=3 MV/m | **39.8 Pa** | **4.06 kg/m²** | 10 g needs 24 cm², 100 g needs 0.025 m², **1 kg needs 0.246 m²**, 70 kg needs 17.2 m²; hand pad 0.05 m² holds **0.20 kg**, panel 0.25 m² holds **1.0 kg** (at breakdown; design with margin) |
| electrostatic in **oil** | ε_r=3, E_bd=20 MV/m | 5.3 kPa | 541 kg/m² | 0.05 m² pad → 27 kg, but liquid-immersed, not wearable in air |
| DEP in air (ε_p=3) | CM=0.4 at 3 MV/m | 15.9 Pa | 1.6 kg/m² | attraction for all ordinary solids in air |
| acoustic | safe SPL | 120 dB → 0.29 g/m²; 140 dB → 28.8 g/m² | — | holding 1 kg over 0.05 m² needs **168 dB**, 70 kg over 0.1 m² needs **184 dB** — unsafe; acoustic = mm-scale/lab traps, or a push jet |
| magnetic diamagnetic (universal) | χ≈−10⁻⁵ | 1 T→4 Pa (0.41 kg/m²); 10 T→398 Pa (40.6 kg/m²); 16 T→1.0 kPa (104 kg/m²) | 0.05 m² pad: 20 g @1 T, 2 kg @10 T, 5.2 kg @16 T | 10–16 T = lab superconducting magnet, not portable |
| ion wind (pushes via air, any object) | v_e=10 m/s | thrust=mg | power ideal = ½ mg v_e | 0.1 kg ≈ 5 W ideal / 25 W real; 1 kg ≈ 49/245 W; **70 kg ≈ 3.4/17 kW** — not wearable |
| Λ / dark handles | — | 5.25×10⁻¹⁰ Pa | 5×10⁻¹¹ kg/m² | dead |

## Genuine assistive design: "UG-1 invisible hand"
**Hardware**
- Electrode array (hand pad or room panel), insulated, high-voltage low-current supply
  with interlocks; grip = induction on conductors / positive DEP on dielectrics in air.
- Ultrasonic phased array for fine positioning and **push/release** (opposite-side
  node) and for mm-scale levitation in a chamber.
- Optional ion micro-jets for gentle contactless nudging of light objects.
- Cameras/LiDAR/IMU + force sensing; closed loop <100 ms.

**Control**
- Mouth-retainer/tongue-drive or non-invasive EEG BCI = **setpoint only** ("cup", "left",
  speed). Autonomy does path planning and collision avoidance.
- Impedance control renders a virtual hand (stiffness/damping) so it feels like a hand.
- Hard limits: payload, altitude/speed, geofence/no-fly zones, kill switch, auto-land on
  signal loss, cybersecurity, medical supervision, consent (`ASSISTIVE_SAFETY.md`).

**Payload budget in air (safe, with margin)**
- Hand pad 500 cm² at 3 kV/1 mm: ~**0.2 kg absolute**, design **50–100 g**.
- Panel 0.25 m² at 30 kV/1 cm: ~**1 kg absolute**, design **0.3–0.5 kg**.
- Heavier than that → soft robotic gripper or maglev track; heavier-than-air human
  levitation is not possible with any of these fields safely (acoustic 184 dB,
  electrostatic 17 m²/300 kV, ion wind 17 kW, gravity 10⁹ kg overhead, Λ dead).

**Pre-registered tests**: F(V,d) for plastic/wood/metal/glass/fruit; payload vs area;
breakdown margin; DEP sign in air (attract) vs oil (can repel); acoustic SPL safety;
BCI latency <100 ms; fail-safe on signal loss; payload bar: ≥100 g handheld, ≥0.3 kg
panel, zero safety incidents.
