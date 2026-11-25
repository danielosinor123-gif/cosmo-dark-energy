# Assistive-safety rules for any mouth-BCI "invisible hand" (non-negotiable)

Scope: assistive use only (e.g. helping a disabled user move objects / navigate a
mobility base). No implants. Not for harm, surveillance, or non-consensual control.

1. **Consent and medical oversight.** Any mouth-worn neural interface (tongue-drive,
   intraoral EEG, etc.) is a medical device: clinician approval, informed consent,
   data-privacy and cybersecurity review before any trial.
2. **BCI is a setpoint, never the actuator.** The safety layer owns the field drivers:
   hard payload/force limits, max altitude and speed, geofenced/no-fly zones, kill
   switch, automatic fail-safe landing/grip-release on signal loss or fault.
3. **High-voltage interlocks.** Electrostatic gripper: insulated electrodes, current
   limiting, arc detection, enclosure, emergency discharge; operate below breakdown
   with margin (air ~3 MV/m; design well under it).
4. **Acoustic safety.** Keep SPL within safe hearing limits outside the work volume;
   beamform away from ears; no macro-payload ultrasound (needs >160 dB).
5. **No weapons, no targeting of people.** Payloads are objects, not humans; no
   human levitation; geofencing and remote disable.
6. **Fail-safe by default.** Loss of BCI, vision, power, or comms → release object
   gently and stop/land.
7. **Auditability.** Logs of intent, setpoints, field commands, and faults; independent
   safety review; pre-registered tests (`UNIVERSAL_LEVITATION.md`).
