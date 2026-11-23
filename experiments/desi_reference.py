"""Extract w0/wa (and base params) posterior stats from DESI's own published DR2 chains."""
import numpy as np
import os

D = r"C:\Users\OSINOR TRILLIONS\cosmo\data\desi_y3_chains"

for fn, label in [("desi_base_compressed_chain.1.txt", "base LCDM (BAO+CMB-compressed)"),
                  ("desi_w0wa_compressed_chain.1.txt", "base_w_wa (BAO+CMB-compressed)")]:
    path = os.path.join(D, fn)
    with open(path) as f:
        header = f.readline().strip()
    cols = header.lstrip("#").split()
    data = np.loadtxt(path)
    print(f"\n=== {label}: {data.shape[0]} samples, {len(cols)} cols")
    idx = {c: i for i, c in enumerate(cols)}
    for p in ["w", "wa", "omegam", "H0", "rdrag", "omegamh2"]:
        if p in idx:
            x = data[:, idx[p]]
            print(f"  {p}: mean={x.mean():.4f} std={x.std():.4f} "
                  f"16/50/84 pct={np.percentile(x, [16, 50, 84]).round(4)}")
    if "wa" in idx and "w" in idx:
        w, wa = data[:, idx["w"]], data[:, idx["wa"]]
        # LambdaCDM point w=-1, wa=0: posterior mass in the quadrant w>-1 & wa<0 etc.
        frac_w_gt_m1 = (w > -1).mean()
        frac_wa_lt_0 = (wa < 0).mean()
        print(f"  P(w>-1)={frac_w_gt_m1:.3f}  P(wa<0)={frac_wa_lt_0:.3f}")
