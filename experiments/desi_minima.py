"""Parse DESI's own iminuit bestfit.minimum.txt files -> dchi2_MAP per combo/model."""
import glob
import os
import numpy as np

D = r"C:\Users\OSINOR TRILLIONS\cosmo\data\desi_y3_minima"
rows = {}
for f in sorted(glob.glob(os.path.join(D, "*.txt"))):
    base = os.path.basename(f)[:-4]
    model, combo = base.split("__", 1)
    with open(f) as fh:
        lines = [l for l in fh if not l.startswith("#")]
    parts = lines[0].split()
    hdr = [t for t in open(f).readline().split() if not t.startswith("#")]
    i_chi2 = hdr.index("chi2")
    i_mlp = hdr.index("minuslogpost")
    i_w = hdr.index("w") if "w" in hdr else None
    i_wa = hdr.index("wa") if "wa" in hdr else None
    rows[(combo, model)] = dict(
        chi2=float(parts[i_chi2]), minuslogpost=float(parts[i_mlp]),
        w=float(parts[i_w]) if i_w is not None else None,
        wa=float(parts[i_wa]) if i_wa is not None else None)

from scipy.stats import chi2 as chi2d, norm
combos = sorted({c for c, m in rows})
print(f"{'combo':60s} {'chi2(base)':>10s} {'chi2(w0wa)':>10s} {'dchi2':>8s} {'sigma':>7s}   w0wa best (w, wa)")
for c in combos:
    if (c, "base") in rows and (c, "base_w_wa") in rows:
        d = rows[(c, "base_w_wa")]["chi2"] - rows[(c, "base")]["chi2"]
        p = chi2d.sf(abs(d), 2)
        sig = abs(norm.isf(p)) if p < 0.5 else 0.0
        w, wa = rows[(c, "base_w_wa")]["w"], rows[(c, "base_w_wa")]["wa"]
        print(f"{c:60s} {rows[(c,'base')]['chi2']:10.3f} {rows[(c,'base_w_wa')]['chi2']:10.3f} "
              f"{d:+8.2f} {sig:6.2f}s  ({w:+.3f}, {wa:+.3f})")
    if (c, "base") in rows and (c, "base_w") in rows:
        d = rows[(c, "base_w")]["chi2"] - rows[(c, "base")]["chi2"]
        p = chi2d.sf(abs(d), 1)
        sig = abs(norm.isf(p)) if p < 0.5 else 0.0
        w = rows[(c, "base_w")]["w"]
        print(f"{c:60s} {rows[(c,'base')]['chi2']:10.3f} {rows[(c,'base_w')]['chi2']:10.3f} "
              f"{d:+8.2f} {sig:6.2f}s  (w={w:+.3f})  [wCDM]")
