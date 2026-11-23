import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from cosmo import Background, BackgroundBinned

# binned w that mimics CPL w0=-0.9, wa=0.5 at low z: w(0)=-0.9, w(1)=-0.4, w(2)~-0.067
edges = [0.0, 0.8, 1.6, np.inf]
wb = [-0.9, -0.4, -0.1]
Om, H0 = 0.31, 67.5
bb = BackgroundBinned(Om, wb, edges)
bc = Background(Om, w0=-0.9, wa=0.5)
z = np.array([0.3, 0.8, 1.2, 2.0, 2.33])
print("z      DM_binned  DM_cpl   DH_binned  DH_cpl   mu_binned mu_cpl")
for zi in z:
    print(f"{zi:5.2f} {bb.DM(zi,H0):9.2f} {bc.DM(zi,H0):8.2f} {bb.DH(zi,H0):9.2f} {bc.DH(zi,H0):8.2f} "
          f"{bb.mu(zi,H0):9.3f} {bc.mu(zi,H0):8.3f}")
# consistency: constant w=-1 bins must reproduce LCDM exactly
bl = BackgroundBinned(Om, [-1.0, -1.0, -1.0], edges)
b0 = Background(Om)
print("max |DM diff| vs LCDM:", np.abs(bl.DM(z, H0) - b0.DM(z, H0)).max())
