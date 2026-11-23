"""cosmo.py — flat w0waCDM background distances + public-likelihood plumbing (numpy/scipy only).

Cosmology:
  E(z)^2 = Om (1+z)^3 + Ok (1+z)^2 + Ode * f(z),
  f(z) = (1+z)^{3(1+w0+wa)} exp(-3 wa z/(1+z))          (CPL / w0wa)
  DM(z) = (c/H0) ∫0^z dz'/E(z')   (transverse comoving, flat)
  DH(z) = c/(H0 E(z)),  DL = (1+z) DM,  mu = 5 log10(DL/Mpc) + 25
  D_V(z) = (z DM^2 DH)^{1/3}
  Sound horizon: exact numerical integral of c_s/H from z_drag (Eisenstein & Hu 1998 fit
  for z_drag) to infinity, including radiation + neutrinos.

Data (all public, cloned into ../data):
  - DESI DR2 BAO, per-tracer and ALL: data/bao_data/desi_bao_dr2/*_{mean,cov}.txt
    mean rows: "z value quantity" with quantity in {DM_over_rs, DH_over_rs, DV_over_rs}
  - DES-Dovekie SN Hubble diagram: data/sn_data/DES-Dovekie/{DES-Dovekie_HD.csv, covtot_inv_000.npz}
  - Pantheon+ (cross-check): data/sn_data/PantheonPlus/{Pantheon+SH0ES.dat, Pantheon+SH0ES_STAT+SYS.cov}
  - Planck 2018 (PR3) TT,TE,EE+lowE compressed prior on (R, lA, omega_b) — used as a diagonal
    Gaussian (documented approximation): R=1.7502+-0.0046, lA=301.471+-0.090, wb=0.02237+-0.00015
    (Planck 2018 results VI, Table 2 values; correlations neglected -> conservative).
"""
import os
import numpy as np

C_KMS = 299792.458          # km/s
TCMB = 2.7255               # K
OGH2 = 2.4728e-5            # photon density parameter (T=2.7255 K)
NEFF = 3.046
ORH2 = OGH2 * (1.0 + 0.2271 * NEFF)   # photons + neutrinos (massless approx)


# ----------------------------------------------------------------------------
# Background distances with a cached integration grid
# ----------------------------------------------------------------------------

_ZG = np.concatenate([np.linspace(0.0, 3.0, 3001), np.geomspace(3.0, 1200.0, 3001)[1:]])


def Ez(z, Om, w0=-1.0, wa=0.0, Ok=0.0, Ode=None):
    z = np.asarray(z, dtype=float)
    if Ode is None:
        Ode = 1.0 - Om - Ok
    f = (1 + z) ** (3 * (1 + w0 + wa)) * np.exp(-3 * wa * z / (1 + z))
    return np.sqrt(Om * (1 + z) ** 3 + Ok * (1 + z) ** 2 + Ode * f)


def _cumint(Om, w0, wa, Ok):
    """I(z) = ∫0^z dz'/E(z') on the fixed grid (vectorized, cached by caller)."""
    E = Ez(_ZG, Om, w0, wa, Ok)
    return _ZG, np.concatenate([[0.0], np.cumsum(0.5 * (1 / E[1:] + 1 / E[:-1]) * np.diff(_ZG))])


class Background:
    """Distances for one (Om, w0, wa, Ok); H0 only rescales DM."""

    def __init__(self, Om, w0=-1.0, wa=0.0, Ok=0.0):
        self.Om, self.w0, self.wa, self.Ok = Om, w0, wa, Ok
        self.zg, self.Ig = _cumint(Om, w0, wa, Ok)

    def DM(self, z, H0):
        z = np.asarray(z, dtype=float)
        return (C_KMS / H0) * np.interp(z, self.zg, self.Ig, left=0.0)

    def DH(self, z, H0):
        return C_KMS / (H0 * Ez(z, self.Om, self.w0, self.wa, self.Ok))

    def DV(self, z, H0):
        z = np.asarray(z, dtype=float)
        return (z * self.DM(z, H0) ** 2 * self.DH(z, H0)) ** (1.0 / 3.0)

    def mu(self, z, H0):
        return 5.0 * np.log10((1 + np.asarray(z, float)) * self.DM(z, H0)) + 25.0


def lnF_binned(z, edges, w_bins):
    """ln of the DE density evolution factor for piecewise-constant w(z).
    edges = [z0, z1, ..., zn] (z0 = 0, zn = inf allowed); last bin extends to infinity."""
    z = np.asarray(z, dtype=float)
    out = np.zeros_like(z)
    for i, wi in enumerate(w_bins):
        a = edges[i]
        b = edges[i + 1] if i + 1 < len(edges) else np.inf
        zhi = np.minimum(z, b)
        out = out + np.where(z > a, 3.0 * (1.0 + wi) * np.log((1.0 + zhi) / (1.0 + a)), 0.0)
    return out


class BackgroundBinned:
    """Distances for piecewise-constant w(z) in bins `edges` (flat universe)."""

    def __init__(self, Om, w_bins, edges, Ok=0.0):
        self.Om, self.Ok = Om, Ok
        self.edges, self.w_bins = list(edges), list(w_bins)
        Ode = 1.0 - Om - Ok
        E = np.sqrt(Om * (1 + _ZG) ** 3 + Ok * (1 + _ZG) ** 2 + Ode * np.exp(lnF_binned(_ZG, edges, w_bins)))
        self.zg = _ZG
        self.Ig = np.concatenate([[0.0], np.cumsum(0.5 * (1 / E[1:] + 1 / E[:-1]) * np.diff(_ZG))])

    def DM(self, z, H0):
        return (C_KMS / H0) * np.interp(z, self.zg, self.Ig, left=0.0)

    def DH(self, z, H0):
        z = np.asarray(z, dtype=float)
        Ode = 1.0 - self.Om - self.Ok
        E = np.sqrt(self.Om * (1 + z) ** 3 + self.Ok * (1 + z) ** 2
                    + Ode * np.exp(lnF_binned(z, self.edges, self.w_bins)))
        return C_KMS / (H0 * E)

    def DV(self, z, H0):
        z = np.asarray(z, dtype=float)
        return (z * self.DM(z, H0) ** 2 * self.DH(z, H0)) ** (1.0 / 3.0)

    def mu(self, z, H0):
        return 5.0 * np.log10((1 + np.asarray(z, float)) * self.DM(z, H0)) + 25.0


# ----------------------------------------------------------------------------
# Linear growth: D(a) and f = dlnD/dlna (radiation neglected in E^2)
# ----------------------------------------------------------------------------

def growth(Om, z_eval, w0=-1.0, wa=0.0, Ok=0.0, lnF=None):
    """Return (D, f) at z_eval, D normalized to 1 at z=0, f = dlnD/dlna.
    Uses the exact growing-mode integral (Heath 1977, radiation neglected):
    D(a) = (5/2) Om E(a) ∫0^a da' / (a'^3 E(a')^3)."""
    ag = np.geomspace(1e-4, 1.0, 4001)
    zg = 1.0 / ag - 1.0
    if lnF is None:
        lnF_z = 3.0 * (1 + w0 + wa) * np.log1p(zg) - 3.0 * wa * zg / (1 + zg)
    else:
        lnF_z = lnF(zg)
    Ode = 1.0 - Om - Ok
    E = np.sqrt(Om * ag ** -3 + Ok * ag ** -2 + Ode * np.exp(lnF_z))
    integ = 1.0 / (ag ** 3 * E ** 3)
    I = np.concatenate([[0.0], np.cumsum(0.5 * (integ[1:] + integ[:-1]) * np.diff(ag))])
    D = 2.5 * Om * E * I
    D = D / D[-1]
    x = np.log(ag)
    f = np.gradient(np.log(np.maximum(D, 1e-300)), x)
    xe = np.log(1.0 / (1.0 + np.asarray(z_eval, dtype=float)))
    return np.interp(xe, x, D), np.interp(xe, x, f)


def fs8_model(Om, s8, z, w0=-1.0, wa=0.0, Ok=0.0, lnF=None):
    """f sigma8(z) = s8 * f(z) * D(z) with D(0)=1."""
    D, f = growth(Om, z, w0, wa, Ok, lnF)
    return s8 * f * D


def load_fs8():
    """Public f sigma8 benchmark (AparajeetS/fsigma8-growth-consistency, Sagredo+2018
    Table I + extensions): z, fs8, inverse total covariance. Published covariance blocks
    (WiggleZ 3x3, eBOSS DR14 QSO 4x4) replace the diagonal entries for those points."""
    d = os.path.join(DATA, "fs8")
    import csv
    rows = list(csv.DictReader(open(os.path.join(d, "catalogue.csv"))))
    rows = [r for r in rows if r["included_in_primary_vector"] == "True"]
    z = np.array([float(r["z"]) for r in rows])
    fs8 = np.array([float(r["fsigma8"]) for r in rows])
    sig = np.array([float(r["sigma"]) for r in rows])
    n = len(z)
    cov = np.diag(sig ** 2)
    blocks = {"wigglez_fs8": "wigglez_fsigma8_covariance.csv",
              "eboss_dr14_qso_fs8": "eboss_dr14_qso_fsigma8_covariance.csv"}
    for bid, fname in blocks.items():
        idx = [i for i, r in enumerate(rows) if r["covariance_block_id"] == bid]
        if not idx:
            continue
        idx = sorted(idx, key=lambda i: z[i])   # block rows are in ascending z
        B = np.loadtxt(os.path.join(d, "covariance_blocks", fname), delimiter=",")
        assert B.shape == (len(idx), len(idx)), (bid, B.shape, len(idx))
        for a, i in enumerate(idx):
            for b, j in enumerate(idx):
                cov[i, j] = B[a, b]
    invC = np.linalg.inv(cov)
    return z, fs8, invC


def fs8_chi2(Om, s8, z, fs8v, invC, w0=-1.0, wa=0.0, Ok=0.0, lnF=None):
    r = fs8v - fs8_model(Om, s8, z, w0, wa, Ok, lnF)
    return float(r @ invC @ r)


# ----------------------------------------------------------------------------
# Sound horizon (drag epoch), numerical
# ----------------------------------------------------------------------------

def z_drag_EH98(Om, Ob, h, Tcmb=TCMB):
    omh2, obh2 = Om * h * h, Ob * h * h
    theta = Tcmb / 2.7
    b1 = 0.313 * omh2 ** -0.419 * (1 + 0.607 * omh2 ** 0.674)
    b2 = 0.238 * omh2 ** 0.223
    return 1291 * omh2 ** 0.251 / (1 + 0.659 * omh2 ** 0.828) * (1 + b1 * obh2 ** b2)


def sound_horizon(Om, Ob, h, Tcmb=TCMB):
    """r_s(z_drag) in Mpc via numerical integration of c_s/H including radiation."""
    zd = z_drag_EH98(Om, Ob, h, Tcmb)
    H0 = 100.0 * h
    omh2, obh2 = Om * h * h, Ob * h * h
    # integrate in x = ln(1+z) from ln(1+zd) to large x
    x = np.linspace(np.log1p(zd), np.log1p(1e7), 40001)
    z = np.expm1(x)
    # H(z) = 100 * sqrt(omh2 (1+z)^3 + orh2 (1+z)^4 + ...): densities are *_h^2, so the
    # prefactor is 100 km/s/Mpc, NOT H0 (which would double-count h).
    Hz = 100.0 * np.sqrt(omh2 * (1 + z) ** 3 + ORH2 * (1 + z) ** 4)
    # NOTE: late-time DE negligible at z >> zd; radiation dominates the integral tail.
    Rb = 0.75 * (obh2 / OGH2) / (1 + z)          # 3 rho_b / 4 rho_gamma
    cs = 1.0 / np.sqrt(3.0 * (1.0 + Rb))          # in units of c
    integrand = cs / Hz * (1 + z)                  # dz = (1+z) dx; rs = c * ∫ (cs/c) dz / H
    rs = C_KMS * np.trapezoid(integrand, x)
    return rs, zd


# ----------------------------------------------------------------------------
# DESI DR2 BAO likelihood
# ----------------------------------------------------------------------------

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "data")


def load_desi_bao(dr2_dir=None, which="ALL"):
    """Returns list of dicts: {name, z[], vals[], inv_cov}. ALL = combined 13-dim vector."""
    dr2_dir = dr2_dir or os.path.join(DATA, "bao_data", "desi_bao_dr2")
    out = []
    if which == "ALL":
        mean = np.genfromtxt(os.path.join(dr2_dir, "desi_gaussian_bao_ALL_GCcomb_mean.txt"),
                             dtype=None, encoding="utf-8", comments="#")
        cov = np.loadtxt(os.path.join(dr2_dir, "desi_gaussian_bao_ALL_GCcomb_cov.txt"))
        z = np.array([float(r[0]) for r in mean])
        vals = np.array([float(r[1]) for r in mean])
        qty = [str(r[2]) for r in mean]
        out.append(dict(name="ALL", z=z, vals=vals, qty=qty, inv_cov=np.linalg.inv(cov)))
    else:
        for f in sorted(os.listdir(dr2_dir)):
            if not f.endswith("_mean.txt") or f.startswith("desi_gaussian_bao_ALL"):
                continue
            name = f.replace("desi_gaussian_bao_", "").replace("_mean.txt", "")
            mean = np.genfromtxt(os.path.join(dr2_dir, f), dtype=None, encoding="utf-8", comments="#")
            if mean.ndim == 0:
                mean = mean.reshape(1)
            nrows = len(mean)
            covf = f.replace("_mean.txt", "_cov.txt")
            cov = np.loadtxt(os.path.join(dr2_dir, covf))
            if cov.ndim == 0:
                cov = cov.reshape(1, 1)
            elif cov.ndim == 1:
                cov = np.diag(cov) if cov.size == nrows else cov.reshape(1, 1)
            out.append(dict(name=name, z=np.array([float(r[0]) for r in mean]),
                            vals=np.array([float(r[1]) for r in mean]),
                            qty=[str(r[2]) for r in mean], inv_cov=np.linalg.inv(cov)))
    return out


def bao_model_vec(bg, H0, rd, rows):
    m = np.empty(len(rows["z"]))
    for i, q in enumerate(rows["qty"]):
        z = rows["z"][i]
        if q == "DM_over_rs":
            m[i] = bg.DM(z, H0) / rd
        elif q == "DH_over_rs":
            m[i] = bg.DH(z, H0) / rd
        elif q == "DV_over_rs":
            m[i] = bg.DV(z, H0) / rd
        else:
            raise ValueError(q)
    return m


def bao_chi2(bg, H0, rd, rows):
    r = bao_model_vec(bg, H0, rd, rows) - rows["vals"]
    return float(r @ rows["inv_cov"] @ r)


# ----------------------------------------------------------------------------
# SN likelihoods
# ----------------------------------------------------------------------------

def load_dovekie():
    """DES-Dovekie Hubble diagram: z, MU, inverse total covariance (npz)."""
    d = os.path.join(DATA, "sn_data", "DES-Dovekie")
    raw = np.genfromtxt(os.path.join(d, "DES-Dovekie_HD.csv"), delimiter=",", names=True)
    z = np.asarray(raw["zHD"], float)
    mu = np.asarray(raw["MU"], float)
    npz = np.load(os.path.join(d, "covtot_inv_000.npz"))
    packed = np.asarray(npz["cov"], float)
    n = len(z)
    assert packed.size == n * (n + 1) // 2, packed.shape
    invC = np.zeros((n, n))
    iu = np.triu_indices(n)          # packed upper triangle WITH diagonal
    invC[iu] = packed
    invC = invC + invC.T - np.diag(np.diag(invC))
    return z, mu, invC


def load_desy5():
    """DES-SN5YR Hubble diagram: zHD, MU, inverse total covariance (cached).
    covsys_000.txt alone is NOT PD (832 negative eigenvalues) — it is the systematics
    block; total cov = covsys + diag(MUERR_FINAL^2) (validated: PD, chi2/dof~1.1)."""
    d = os.path.join(DATA, "sn_data", "DESY5")
    cache = os.path.join(d, "covtot_inv.npy")
    raw = np.genfromtxt(os.path.join(d, "DES-SN5YR_HD.csv"), delimiter=",", names=True)
    z = np.asarray(raw["zHD"], float)
    mu = np.asarray(raw["MU"], float)
    merr = np.asarray(raw["MUERR_FINAL"], float)
    if os.path.exists(cache):
        invC = np.load(cache)
    else:
        raw_cov = np.loadtxt(os.path.join(d, "covsys_000.txt"))
        if raw_cov.ndim == 1:
            n = int(raw_cov[0])
            cov = raw_cov[1:1 + n * n].reshape(n, n)
        else:
            cov = raw_cov
        assert cov.shape == (len(z), len(z)), cov.shape
        cov = cov + np.diag(merr ** 2)
        invC = np.linalg.inv(cov)
        np.save(cache, invC)
    return z, mu, invC


def load_pantheonplus():
    """Pantheon+ SH0ES: zHD, m_b_corr, STAT+SYS covariance (inverse cached to .npy).
    Note: the mB column does NOT give a sane chi2 with this covariance; m_b_corr + zHD
    gives chi2/dof ~ 1.04 at a Planck fiducial (validated in session 3)."""
    d = os.path.join(DATA, "sn_data", "PantheonPlus")
    cache = os.path.join(d, "cov_statsys.npy")
    if os.path.exists(cache):
        invC = np.load(cache)
    else:
        raw_cov = np.loadtxt(os.path.join(d, "Pantheon+SH0ES_STAT+SYS.cov"))
        if raw_cov.ndim == 1:
            # format: first entry = dimension n, then n*n matrix row-major
            n = int(raw_cov[0])
            cov = raw_cov[1:1 + n * n].reshape(n, n)
        else:
            cov = raw_cov
        invC = np.linalg.inv(cov)
        np.save(cache, invC)
    raw = np.genfromtxt(os.path.join(d, "Pantheon+SH0ES.dat"), names=True, dtype=None, encoding="utf-8")
    z = np.asarray(raw["zHD"], float)
    mB = np.asarray(raw["m_b_corr"], float)
    return z, mB, invC


def sn_hd_chi2(bg, H0, z, mu, invC):
    """Hubble-diagram form (M already marginalized): chi2 = r^T C^-1 r."""
    r = mu - bg.mu(z, H0)
    return float(r @ invC @ r)


def sn_mB_chi2(bg, H0, z, mB, invC):
    """mB form with analytic marginalization over M: chi2 = rCr - (1Cr)^2/(1C1)."""
    r = mB - bg.mu(z, H0)
    one = np.ones_like(r)
    C1 = invC @ one
    return float(r @ invC @ r - (r @ C1) ** 2 / (one @ C1))


# ----------------------------------------------------------------------------
# Planck 2018 compressed prior (diagonal Gaussian approximation, documented)
# ----------------------------------------------------------------------------

PLANCK18_COMPRESSED = dict(
    R=(1.7502, 0.0046), lA=(301.471, 0.090), wb=(0.02237, 0.00015),  # TT,TE,EE+lowE
    zstar=1089.9,
)


def cmb_compressed_chi2(Om, Ob, h, bg=None, use_lA=False):
    """Planck 2018 compressed prior as diagonal Gaussian.
    Default uses (R, omega_b) only: lA = pi DM(z*)/rs needs rs, and our EH98-z_d +
    numerical-integral rs is ~2-3% off CAMB, which would be ~100 sigma in lA.
    R = sqrt(Om) H0 DM(z*)/c and omega_b = Ob h^2 are rs-independent."""
    p = PLANCK18_COMPRESSED
    H0 = 100.0 * h
    if bg is None:
        bg = Background(Om)
    DMstar = bg.DM(p["zstar"], H0)
    R = np.sqrt(Om) * H0 * DMstar / C_KMS
    wb = Ob * h * h
    chi2 = ((R - p["R"][0]) / p["R"][1]) ** 2
    chi2 += ((wb - p["wb"][0]) / p["wb"][1]) ** 2
    if use_lA:
        rs, _ = sound_horizon(Om, Ob, h)
        lA = np.pi * DMstar / rs
        chi2 += ((lA - p["lA"][0]) / p["lA"][1]) ** 2
    return float(chi2)


# ----------------------------------------------------------------------------
# Tiny adaptive Metropolis sampler (numpy only)
# ----------------------------------------------------------------------------

def metropolis(logpost, x0, nsteps=30000, burn=6000, step=None, seed=0, adapt_every=250):
    rng = np.random.default_rng(seed)
    x = np.array(x0, dtype=float)
    lp = logpost(x)
    dim = len(x)
    if step is None:
        step = 0.02 * np.ones(dim)
    cov = np.diag(step ** 2)
    chol = np.linalg.cholesky(cov + 1e-12 * np.eye(dim))
    samples = []
    acc = 0
    recent = []
    for i in range(nsteps):
        prop = x + chol @ rng.standard_normal(dim)
        lpp = logpost(prop)
        if np.log(rng.random()) < lpp - lp:
            x, lp = prop, lpp
            acc += 1
        if i >= burn:
            samples.append(x.copy())
        recent.append(x.copy())
        if i > burn and i % adapt_every == 0:
            arr = np.array(recent[-adapt_every:])
            if len(arr) > dim + 2:
                cov = np.cov(arr.T) + 1e-10 * np.eye(dim)
                chol = np.linalg.cholesky(cov)
    return np.array(samples), acc / nsteps
