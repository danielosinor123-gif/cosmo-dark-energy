import os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import cosmo
from cosmo import Background, bao_chi2, sn_hd_chi2, cmb_compressed_chi2, load_desi_bao, load_dovekie

# fiducial Planck-ish LCDM
Om, Ob, h, rd = 0.315, 0.0224 / 0.674 ** 2, 0.674, 147.1
H0 = 100 * h
t0 = time.time()
rs, zd = cosmo.sound_horizon(Om, Ob, h)
print(f"rs = {rs:.3f} Mpc (expect ~147), zd = {zd:.1f} (expect ~1059), took {time.time()-t0:.2f}s")

bg = Background(Om)
bao = load_desi_bao(which="ALL")[0]
print("BAO chi2 (LCDM fid):", bao_chi2(bg, H0, rd, bao))
z, mu, invC = load_dovekie()
t0 = time.time()
print("SN Dovekie chi2:", sn_hd_chi2(bg, H0, z, mu, invC), f"({time.time()-t0:.2f}s, n={len(z)})")
t0 = time.time()
print("CMB compressed chi2:", cmb_compressed_chi2(Om, Ob, h, bg), f"({time.time()-t0:.2f}s)")
# w0wa point
bg2 = Background(Om, w0=-0.8, wa=-0.8)
print("BAO chi2 (w0=-0.8,wa=-0.8):", bao_chi2(bg2, H0, rd, bao))
print("SN chi2 (w0wa):", sn_hd_chi2(bg2, H0, z, mu, invC))
print("CMB chi2 (w0wa):", cmb_compressed_chi2(Om, Ob, h, bg2))
# timing of one full logpost-like eval
t0 = time.time()
for _ in range(20):
    b = Background(0.31, -0.9, 0.4)
    bao_chi2(b, H0, rd, bao); sn_hd_chi2(b, H0, z, mu, invC); cmb_compressed_chi2(Om, Ob, h, b)
print(f"full eval: {(time.time()-t0)/20*1000:.1f} ms")
