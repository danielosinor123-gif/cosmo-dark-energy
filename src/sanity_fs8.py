import os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cosmo
from cosmo import growth, fs8_model, load_fs8, fs8_chi2

z, fs8v, invC = load_fs8()
print(f"fs8: {len(z)} points, z range {z.min()}..{z.max()}, invC shape {invC.shape}")

t0 = time.time()
D, f = growth(0.31, np.array([0.0, 0.5, 1.0, 2.0]))
print(f"growth LCDM Om=0.31: D={np.round(D,4)} f={np.round(f,4)}  ({time.time()-t0:.3f}s)")
# LCDM checks: D(0)=1, f(0)=Om^0.55 ~ 0.31^0.55 = 0.525
print("f(0) vs Om^0.55:", f[0], 0.31 ** 0.55)

t0 = time.time()
m = fs8_model(0.31, 0.81, z)
chi2 = fs8_chi2(0.31, 0.81, z, fs8v, invC)
print(f"fs8 model vs data at Planck-ish fiducial: chi2={chi2:.2f}/{len(z)} ({time.time()-t0:.3f}s)")
print("z:      ", np.round(z, 3))
print("data:   ", np.round(fs8v, 3))
print("model:  ", np.round(m, 3))

# w0wa sanity: evolving DE changes growth
m2 = fs8_model(0.31, 0.81, z, w0=-0.8, wa=-0.8)
print("model w0wa(-0.8,-0.8):", np.round(m2, 3))
D2, f2 = growth(0.31, np.array([0.5, 1.0]), w0=-0.8, wa=-0.8)
print("w0wa D, f at 0.5,1.0:", np.round(D2, 4), np.round(f2, 4))
