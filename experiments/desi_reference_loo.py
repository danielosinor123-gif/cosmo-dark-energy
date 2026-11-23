"""w0/wa posterior stats from DESI's own DR2 w0wa chains for different data combos."""
import numpy as np
import os

D = r"C:\Users\OSINOR TRILLIONS\cosmo\data\desi_y3_chains"
files = {
    "BAO-only (desi-bao-all)": "w0wa__desi-bao-all.txt",
    "BAO+DES-SN5YR": "w0wa__desi-bao-all_desy5sn.txt",
    "SN-only (desy5sn)": "w0wa__desy5sn.txt",
    "BAO+Pantheon+": "w0wa__desi-bao-all_pantheonplus.txt",
    "BAO+CMB-compressed": "desi_w0wa_compressed_chain.1.txt",
}
for label, fn in files.items():
    path = os.path.join(D, fn)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"{label}: missing"); continue
    with open(path) as f:
        cols = f.readline().strip().lstrip("#").split()
    data = np.loadtxt(path)
    idx = {c: i for i, c in enumerate(cols)}
    w, wa = data[:, idx["w"]], data[:, idx["wa"]]
    print(f"{label:28s} n={data.shape[0]:6d}  w={w.mean():+.3f}+-{w.std():.3f}  "
          f"wa={wa.mean():+.3f}+-{wa.std():.3f}  P(w>-1)={(w>-1).mean():.3f} P(wa<0)={(wa<0).mean():.3f}")
