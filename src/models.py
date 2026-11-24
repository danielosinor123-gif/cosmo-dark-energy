"""models.py — candidate NEW dark-energy physics (session 6).

All models are smooth, non-interacting DE backgrounds specified by
lnF(z) = ln(rho_DE(z)/rho_DE(0)); the generic BackgroundFn and the Heath-integral
growth (valid for any smooth E(a)) apply unchanged.

N1  LDPT — Late Dark-sector Phase Transition:
    rho_DE(z) = rho_DE0 * g(z), g(0) = 1,
    g(z) = [1 - f*sig((z-zc)/d)] / [1 - f*sig(-zc/d)],  sig = logistic.
    Physics: a first-order phase transition in the dark sector at z ~ zc converts a
    fraction f of the vacuum energy into dark matter / radiation ("late vacuum
    metamorphosis"). Predicts a step-like feature in rho_DE(z), hence in w_eff(z),
    and extra matter at z > zc -> modified growth (tested against f sigma8).
    w_eff(z) = -1 - (1/3) d ln rho_DE / d ln a.

N2  QUINTOM-DIP — localized phantom-divide crossing:
    w(z) = wb + A * exp(-((z-zc)/d)^2)
    A < 0 gives a phantom dip (w < -1) around z ~ zc, with w(0) > -1 and w -> wb at
    high z. Single-field scalar models cannot cross w = -1 (ghost/gradient
    instability), so the physical claim is a two-field (quintom) dark sector with a
    divide crossing localized at z ~ 1.
"""
import numpy as np

from cosmo import _ZG, C_KMS


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -500, 500)))


# --- N1: LDPT ---------------------------------------------------------------
def lnF_ldpt(z, f, zc, d):
    z = np.asarray(z, dtype=float)
    g = (1.0 - f * sigmoid((z - zc) / d)) / (1.0 - f * sigmoid(-zc / d))
    return np.log(np.maximum(g, 1e-300))


def w_eff_ldpt(z, f, zc, d):
    """w_eff(z) = -1 - (1/3) d ln rho_DE / d ln a = -1 + (1+z)/3 * d ln g / dz."""
    z = np.asarray(z, dtype=float)
    h = 1e-4 * np.maximum(1.0, z)
    lp = lnF_ldpt(z + h, f, zc, d)
    lm = lnF_ldpt(z - h, f, zc, d)
    dlnf_dz = (lp - lm) / (2.0 * h)
    return -1.0 + (1.0 + z) / 3.0 * dlnf_dz


# --- N2: quintom dip --------------------------------------------------------
def w_quintom(z, wb, A, zc, d):
    z = np.asarray(z, dtype=float)
    return wb + A * np.exp(-((z - zc) / d) ** 2)


def make_lnF_quintom(wb, A, zc, d):
    """lnF(z) = 3 ∫_0^z (1+w(z'))/(1+z') dz', precomputed on the global grid."""
    y = 3.0 * (1.0 + w_quintom(_ZG, wb, A, zc, d)) / (1.0 + _ZG)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(_ZG))])

    def lnF(z):
        return np.interp(np.asarray(z, dtype=float), _ZG, cum, left=0.0)
    return lnF


# --- generic background from any lnF ----------------------------------------
class BackgroundFn:
    """Distances for arbitrary lnF(z) (flat universe)."""

    def __init__(self, Om, lnF, Ok=0.0):
        self.Om, self.Ok, self.lnF = Om, Ok, lnF
        Ode = 1.0 - Om - Ok
        E = np.sqrt(Om * (1.0 + _ZG) ** 3 + Ok * (1.0 + _ZG) ** 2
                    + Ode * np.exp(lnF(_ZG)))
        self.zg = _ZG
        self.Ig = np.concatenate([[0.0], np.cumsum(0.5 * (1 / E[1:] + 1 / E[:-1]) * np.diff(_ZG))])

    def E(self, z):
        z = np.asarray(z, dtype=float)
        Ode = 1.0 - self.Om - self.Ok
        return np.sqrt(self.Om * (1 + z) ** 3 + self.Ok * (1 + z) ** 2
                       + Ode * np.exp(self.lnF(z)))

    def DM(self, z, H0):
        z = np.asarray(z, dtype=float)
        return (C_KMS / H0) * np.interp(z, self.zg, self.Ig, left=0.0)

    def DH(self, z, H0):
        return C_KMS / (H0 * self.E(z))

    def DV(self, z, H0):
        z = np.asarray(z, dtype=float)
        return (z * self.DM(z, H0) ** 2 * self.DH(z, H0)) ** (1.0 / 3.0)

    def mu(self, z, H0):
        return 5.0 * np.log10((1 + np.asarray(z, float)) * self.DM(z, H0)) + 25.0
