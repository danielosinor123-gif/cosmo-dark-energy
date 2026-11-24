"""desi_2pt.py — DESI's own chains: does BAO + full-shape (2pt) keep the w0wa signal?

Parses DESI Y3 cobaya chains (data.desi.lbl.gov) and compares w, wa posteriors:
  w0wa__desi-bao-all.txt            BAO only
  w0wa__bao_all_desy3x2pt.txt       BAO + DESI Y3 galaxy 2pt (full shape, no SN)
  w0wa__desi-bao-all_desy5sn.txt    BAO + DESY5 SN
  w0wa__desi-bao-all_pantheonplus.txt  BAO + Pantheon+
plus chain.margestats for the same combos when available.
"""
import os
import numpy as np

D = r"C:\Users\OSINOR TRILLIONS\cosmo\data\desi_y3_chains"

FILES = {
    "BAO only": "w0wa__desi-bao-all.txt",
    "BAO + 2pt (full shape)": "w0wa__bao_all_desy3x2pt.txt",
    "BAO + DESY5 SN": "w0wa__desi-bao-all_desy5sn.txt",
    "BAO + Pantheon+": "w0wa__desi-bao-all_pantheonplus.txt",
}

print(f"{'combo':28s} {'n':>8s} {'w mean+-sd':>18s} {'wa mean+-sd':>18s} "
      f"{'P(w>-1)':>8s} {'P(wa<0)':>8s} {'corr(w,wa)':>10s}")
out = {}
for tag, fn in FILES.items():
    p = os.path.join(D, fn)
    if not os.path.exists(p):
        print(f"{tag}: MISSING {fn}")
        continue
    with open(p) as f:
        hdr = ""
        for _ in range(5):
            line = f.readline()
            if "weight" in line:
                hdr = line
                break
    names = hdr.lstrip("#").split()
    iw, iwa = names.index("w"), names.index("wa")
    iwt = names.index("weight") if "weight" in names else 0
    ch = np.loadtxt(p, comments="#")
    w = ch[:, iw] * ch[:, iwt] if iwt else ch[:, iw]
    wa = ch[:, iwa] * ch[:, iwt] if iwt else ch[:, iwa]
    # normalize by weights
    if iwt:
        w = ch[:, iw]; wa = ch[:, iwa]; wt = ch[:, iwt]
        wm = np.average(w, weights=wt); wam = np.average(wa, weights=wt)
        ws = np.sqrt(np.average((w - wm) ** 2, weights=wt))
        was = np.sqrt(np.average((wa - wam) ** 2, weights=wt))
        pw = float(np.sum(wt[w > -1]) / np.sum(wt))
        pwa = float(np.sum(wt[wa < 0]) / np.sum(wt))
        corr = float(np.average((w - wm) * (wa - wam), weights=wt) / (ws * was))
    else:
        wm, ws = w.mean(), w.std(); wam, was = wa.mean(), wa.std()
        pw = float((w > -1).mean()); pwa = float((wa < 0).mean())
        corr = float(np.corrcoef(w, wa)[0, 1])
    out[tag] = dict(w_mean=wm, w_std=ws, wa_mean=wam, wa_std=was, P_w_gt_m1=pw,
                    P_wa_lt_0=pwa, corr=corr, n=len(w))
    print(f"{tag:28s} {len(w):8d} {wm:+.3f}+-{ws:.3f}   {wam:+.3f}+-{was:.3f}   "
          f"{pw:8.3f} {pwa:8.3f} {corr:+10.2f}")

import json
with open(os.path.join(D, os.pardir, os.pardir, "results", "desi_2pt.json"), "w") as f:
    json.dump(out, f, indent=2)
print("2PT DONE")
