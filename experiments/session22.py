"""session22.py — the z~1 dip, per SN sample: 5-bin w(z)+growth with BAO(ALL+Lya DR2 AP)
+ CMB + fsigma8, run separately on Dovekie, DESY5, and non-calibrator Pantheon+.

If the phantom dip at z in [0.75,1.5] (sections 24/31/32) appears in ALL THREE
independent SN samples, it is the strongest surviving evidence for non-Lambda DE.
If it lives in only one sample, it is sample/calibration-like.

Usage: python experiments/session22.py dip
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import BackgroundBinned, lnF_binned, bao_chi2, sn_hd_chi2, sn_mB_chi2, \
    cmb_compressed_chi2, fs8_chi2, metropolis

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

# Ly-a DR2 AP rows (session21)
Z_AP = 2.33
cov_ap = np.array([[0.066 ** 2, 0.225 * 0.066 * 0.33],
                   [0.225 * 0.066 * 0.33, 0.33 ** 2]])
LYA_AP = dict(name="Lya_DR2_AP", z=np.array([Z_AP, Z_AP]),
              vals=np.array([8.600, 39.32]), qty=["DH_over_rs", "DM_over_rs"],
              inv_cov=np.linalg.inv(cov_ap))


def concat_rows(a, b):
    n = len(a["z"])
    inv = np.zeros((n + len(b["z"]),) * 2)
    inv[:n, :n] = a["inv_cov"]
    inv[n:, n:] = b["inv_cov"]
    return dict(name=a["name"] + "+lyaAP", z=np.concatenate([a["z"], b["z"]]),
                vals=np.concatenate([a["vals"], b["vals"]]),
                qty=list(a["qty"]) + list(b["qty"]), inv_cov=inv)


print("loading data ...", flush=True)
bao = concat_rows(cosmo.load_desi_bao(which="ALL")[0], LYA_AP)
z_d, mu_d, invC_d = cosmo.load_dovekie()
z_y, mu_y, invC_y = cosmo.load_desy5()
z_p, mB_p, invC_p = cosmo.load_pantheonplus()
raw = np.genfromtxt(os.path.join(cosmo.DATA, "sn_data", "PantheonPlus", "Pantheon+SH0ES.dat"),
                    names=True, dtype=None, encoding="utf-8")
keep = np.asarray(raw["IS_CALIBRATOR"], float) <= 0.5
idx = np.where(keep)[0]
z_pn = z_p[keep]
mB_pn = mB_p[keep]
invC_pn = np.linalg.inv(np.linalg.inv(invC_p)[np.ix_(idx, idx)])
z_f, fs8v, invC_f = cosmo.load_fs8()

SAMPLES = {
    "dovekie": ("hd", z_d, mu_d, invC_d),
    "desy5": ("hd", z_y, mu_y, invC_y),
    "pantheon_nocal": ("mB", z_pn, mB_pn, invC_pn),
}

EDGES5 = [0.0, 0.3, 0.75, 1.5, 2.5, np.inf]
PRIORS = {"Om": (0.05, 0.9), "Ob": (0.005, 0.10), "h": (0.55, 0.85), "rd": (120.0, 165.0),
          "s8": (0.2, 1.4)}
for i in range(5):
    PRIORS[f"w{i}"] = (-3.0, 1.0)
names = ["Om", "Ob", "h", "rd", "s8"] + [f"w{i}" for i in range(5)]
LBL = ["0-0.3", "0.3-0.75", "0.75-1.5", "1.5-2.5", ">2.5"]


def make_logpost(kind, z_s, m_s, invC_s):
    def logpost(x):
        d = dict(zip(names, x))
        for n in names:
            lo, hi = PRIORS[n]
            if not (lo <= d[n] <= hi):
                return -np.inf
        Om, Ob, h, rd, s8 = d["Om"], d["Ob"], d["h"], d["rd"], d["s8"]
        wb = [d[f"w{i}"] for i in range(5)]
        bg = BackgroundBinned(Om, wb, EDGES5)
        lnF = lambda z: lnF_binned(z, EDGES5, wb)
        H0 = 100.0 * h
        chi2 = bao_chi2(bg, H0, rd, bao) + cmb_compressed_chi2(Om, Ob, h, bg) \
            + fs8_chi2(Om, s8, z_f, fs8v, invC_f, lnF=lnF)
        if kind == "hd":
            chi2 += sn_hd_chi2(bg, H0, z_s, m_s, invC_s)
        else:
            chi2 += sn_mB_chi2(bg, H0, z_s, m_s, invC_s)
        return -0.5 * chi2
    return logpost


out = {"edges": EDGES5, "bao": "ALL + Lya DR2 AP", "samples": {}}
x0 = np.array([0.31, 0.049, 0.68, 147.0, 0.81] + [-0.93, -0.97, -1.2, -1.0, -1.1])
step = np.array([0.01, 0.005, 0.01, 0.5, 0.02] + [0.06] * 5)
for sname, (kind, z_s, m_s, invC_s) in SAMPLES.items():
    logpost = make_logpost(kind, z_s, m_s, invC_s)
    t0 = time.time()
    samples, acc = metropolis(logpost, x0, nsteps=22000, burn=6000, step=step, seed=202)
    stride = max(1, len(samples) // 2500)
    sub = samples[::stride]
    lps = np.array([logpost(s) for s in sub])
    wi = sub[:, 5:10]
    entry = {"n_sn": int(len(z_s)), "chi2_map": float(-2 * lps.max()), "acceptance": float(acc),
             "w_bins_mean": wi.mean(axis=0).tolist(), "w_bins_std": wi.std(axis=0).tolist(),
             "seconds": round(time.time() - t0, 1)}
    out["samples"][sname] = entry
    np.save(os.path.join(RESULTS, f"chain_s22_cert5_{sname}.npy"), np.column_stack([sub, lps]))
    print(f"  [{sname}] n={len(z_s)} chi2_map={entry['chi2_map']:.2f} acc={acc:.2f}: "
          + " | ".join(f"{l}: {m:+.3f}+-{s:.3f}" for l, m, s in
                      zip(LBL, wi.mean(axis=0), wi.std(axis=0)))
          + f" ({entry['seconds']}s)", flush=True)

with open(os.path.join(RESULTS, "dip_persample_results.json"), "w") as f:
    json.dump(out, f, indent=2)
print("DIP DONE", flush=True)
