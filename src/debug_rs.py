import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import cosmo
from cosmo import C_KMS, OGH2, ORH2

Om, Ob, h = 0.315, 0.0224 / 0.674 ** 2, 0.674
omh2, obh2 = Om * h * h, Ob * h * h
H0 = 100 * h
zd = cosmo.z_drag_EH98(Om, Ob, h)
print("zd =", zd, " omh2 =", omh2, " obh2 =", obh2)

# brute-force integral in z (not x)
z = np.geomspace(zd, 1e7, 400000)
Hz = H0 * np.sqrt(omh2 * (1 + z) ** 3 + ORH2 * (1 + z) ** 4)
Rb = 0.75 * (obh2 / OGH2) / (1 + z)
cs = 1.0 / np.sqrt(3.0 * (1.0 + Rb))
rs_z = C_KMS * np.trapezoid(cs / Hz, z)
print("rs (dz integral from zd) =", rs_z)

# EH98 analytic rs for comparison (the formula from the original cosmo.py)
omh2_, obh2_ = omh2, obh2
theta = cosmo.TCMB / 2.7
zeq = 2.50e4 * omh2_ * theta ** -4
keq = 7.46e-2 * omh2_ * theta ** -2
b1 = 0.313 * omh2_ ** -0.419 * (1 + 0.607 * omh2_ ** 0.674)
b2 = 0.238 * omh2_ ** 0.223
zd2 = 1291 * omh2_ ** 0.251 / (1 + 0.659 * omh2_ ** 0.828) * (1 + b1 * obh2_ ** b2)
RE = 31.5 * obh2_ * theta ** -4 * (1e3 / zeq)
RD = 31.5 * obh2_ * theta ** -4 * (1e3 / zd2)
rs_eh = (2.0 / (3.0 * keq)) * np.sqrt(6.0 / RE) * np.log((np.sqrt(1 + RD) + np.sqrt(RD + RE)) / (1 + np.sqrt(RE)))
print("EH98: zd =", zd2, " rs =", rs_eh)

# integral from EH98 zd2 instead of my zd
z = np.geomspace(zd2, 1e7, 400000)
Hz = H0 * np.sqrt(omh2 * (1 + z) ** 3 + ORH2 * (1 + z) ** 4)
Rb = 0.75 * (obh2 / OGH2) / (1 + z)
cs = 1.0 / np.sqrt(3.0 * (1.0 + Rb))
print("rs (dz integral from EH98 zd2) =", C_KMS * np.trapezoid(cs / Hz, z))

# also check: maybe include DE/curvature term properly (should be negligible)
# and check neutrino contribution to H (massless approx ok)

# --- dovekie cov unpack test ---
d = os.path.join(cosmo.DATA, "sn_data", "DES-Dovekie")
raw = np.genfromtxt(os.path.join(d, "DES-Dovekie_HD.csv"), delimiter=",", names=True)
n = len(raw)
npz = np.load(os.path.join(d, "covtot_inv_000.npz"))
packed = np.asarray(npz["cov"], float)
print("n =", n, "packed size =", packed.size, "n(n+1)/2 =", n * (n + 1) // 2)
M = np.zeros((n, n))
iu = np.triu_indices(n)
M[iu] = packed
M = M + M.T - np.diag(np.diag(M))
print("symmetry err:", np.abs(M - M.T).max(), "diag min:", np.diag(M).min())
w = np.linalg.eigvalsh(M)
print("eig min/max:", w.min(), w.max())
