"""session28.py — move the bar up: strongest honest dip measurement from public data.

(a) UNION cert5: DESY5 (HD form) + Pantheon+ non-calibrator (mB form) DEDUPLICATED by
    CID (Pantheon+ objects already in DESY5 removed — no double counting), with
    BAO(ALL + Ly-a DR2 AP) + CMB + fsigma8. 5-bin w(z), edges A.
(b) Correlation scan: combine the three per-sample dip estimates
    (dovekie/desy5/pantheon_nocal from s22) over inter-sample correlation rho in
    [0,1] -> combined |w+1|/sigma range. rho=1 = fully shared systematics (most
    conservative), rho=0 = fully independent.

Usage: python experiments/session28.py union
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

Z_AP = 2.33
cov_ap = np.array([[0.066 ** 2, 0.225 * 0.066 * 0.33],
                   [0.225 * 0.066 * 0.33, 0.33 ** 2]])
LYA_AP = dict(z=np.array([Z_AP, Z_AP]), vals=np.array([8.600, 39.32]),
              qty=["DH_over_rs", "DM_over_rs"], inv_cov=np.linalg.inv(cov_ap))


def concat(a, b):
    n = len(a["z"])
    inv = np.zeros((n + len(b["z"]),) * 2)
    inv[:n, :n] = a["inv_cov"]
    inv[n:, n:] = b["inv_cov"]
    return dict(z=np.concatenate([a["z"], b["z"]]), vals=np.concatenate([a["vals"], b["vals"]]),
                qty=list(a["qty"]) + list(b["qty"]), inv_cov=inv)


def union_mode():
    print("loading + deduping ...", flush=True)
    bao = concat(cosmo.load_desi_bao(which="ALL")[0], LYA_AP)
    z_y, mu_y, invC_y = cosmo.load_desy5()
    z_p, mB_p, invC_p = cosmo.load_pantheonplus()
    raw = np.genfromtxt(os.path.join(cosmo.DATA, "sn_data", "PantheonPlus", "Pantheon+SH0ES.dat"),
                        names=True, dtype=None, encoding="utf-8")
    # CID column has mixed numeric/string ids -> read as str, column 0, skip header
    desy5_cids = set(np.loadtxt(os.path.join(cosmo.DATA, "sn_data", "DESY5", "DES-SN5YR_HD.csv"),
                                delimiter=",", usecols=0, dtype=str, skiprows=1).tolist())
    cal = np.asarray(raw["IS_CALIBRATOR"], float) <= 0.5
    pan_cid = raw["CID"].astype(str)
    keep = np.array([c not in desy5_cids for c in pan_cid]) & cal
    idx = np.where(keep)[0]
    z_pn, mB_pn = z_p[keep], mB_p[keep]
    invC_pn = np.linalg.inv(np.linalg.inv(invC_p)[np.ix_(idx, idx)])
    z_f, fs8v, invC_f = cosmo.load_fs8()
    print(f"  DESY5: {len(z_y)} | Pantheon+ nocal not in DESY5: {keep.sum()} "
          f"(dropped {int((~keep & cal).sum())} dup CIDs + {int((~cal).sum())} calibrators)",
          flush=True)

    EDGES = [0.0, 0.3, 0.75, 1.5, 2.5, np.inf]
    PRIORS = {"Om": (0.05, 0.9), "Ob": (0.005, 0.10), "h": (0.55, 0.85), "rd": (120.0, 165.0),
              "s8": (0.2, 1.4)}
    for i in range(5):
        PRIORS[f"w{i}"] = (-3.0, 1.0)
    names = ["Om", "Ob", "h", "rd", "s8"] + [f"w{i}" for i in range(5)]

    def logpost(x):
        d = dict(zip(names, x))
        for n in names:
            lo, hi = PRIORS[n]
            if not (lo <= d[n] <= hi):
                return -np.inf
        Om, Ob, h, rd, s8 = d["Om"], d["Ob"], d["h"], d["rd"], d["s8"]
        wb = [d[f"w{i}"] for i in range(5)]
        bg = BackgroundBinned(Om, wb, EDGES)
        lnF = lambda z: lnF_binned(z, EDGES, wb)
        H0 = 100.0 * h
        return -0.5 * (bao_chi2(bg, H0, rd, bao) + sn_hd_chi2(bg, H0, z_y, mu_y, invC_y)
                       + sn_mB_chi2(bg, H0, z_pn, mB_pn, invC_pn)
                       + cmb_compressed_chi2(Om, Ob, h, bg)
                       + fs8_chi2(Om, s8, z_f, fs8v, invC_f, lnF=lnF))

    x0 = np.array([0.31, 0.049, 0.68, 147.0, 0.81] + [-0.93, -0.97, -1.2, -1.0, -1.1])
    step = np.array([0.01, 0.005, 0.01, 0.5, 0.02] + [0.06] * 5)
    t0 = time.time()
    samples, acc = metropolis(logpost, x0, nsteps=24000, burn=7000, step=step, seed=288)
    stride = max(1, len(samples) // 2500)
    sub = samples[::stride]
    lps = np.array([logpost(s) for s in sub])
    wi = sub[:, 5:10]
    LBL = ["0-0.3", "0.3-0.75", "0.75-1.5", "1.5-2.5", ">2.5"]
    dip_m, dip_s = wi[:, 2].mean(), wi[:, 2].std()
    print(f"  UNION cert5: chi2_map={-2*lps.max():.2f} acc={acc:.2f} "
          + " | ".join(f"{l}: {m:+.3f}+-{s:.3f}" for l, m, s in zip(LBL, wi.mean(0), wi.std(0)))
          + f" ({time.time()-t0:.0f}s)", flush=True)
    print(f"  UNION DIP z in [0.75,1.5]: w = {dip_m:+.3f} +- {dip_s:.3f} "
          f"-> |w+1|/sigma = {abs(dip_m+1)/dip_s:.2f}", flush=True)
    np.save(os.path.join(RESULTS, "chain_s28_union_cert5.npy"), np.column_stack([sub, lps]))

    # (b) correlation scan combining the three s22 per-sample dip estimates
    d = json.load(open(os.path.join(RESULTS, "dip_persample_results.json")))
    ests = [(e["w_bins_mean"][2], e["w_bins_std"][2]) for e in d["samples"].values()]
    rows = []
    for rho in [0.0, 0.25, 0.5, 0.75, 1.0]:
        ms = np.array([e[0] for e in ests])
        ss = np.array([e[1] for e in ests])
        n = len(ms)
        var = (ss ** 2).sum() + 2 * rho * sum(ss[i] * ss[j] for i in range(n) for j in range(i + 1, n))
        m = ms.mean()
        sig = abs(m + 1) / (np.sqrt(var) / n)
        rows.append({"rho": rho, "combined_w": float(m), "combined_sigma": float(np.sqrt(var) / n),
                     "significance": float(sig)})
        print(f"  rho={rho:.2f}: w={m:+.3f} +- {np.sqrt(var)/n:.3f} -> {sig:.2f}s", flush=True)

    out = {"union_cert5": {"chi2_map": float(-2 * lps.max()), "acceptance": float(acc),
                           "w_bins_mean": wi.mean(0).tolist(), "w_bins_std": wi.std(0).tolist(),
                           "dip_w": float(dip_m), "dip_sigma": float(dip_s),
                           "dip_significance": float(abs(dip_m + 1) / dip_s),
                           "n_desy5": int(len(z_y)), "n_pantheon_nocal_dedup": int(keep.sum())},
           "correlation_scan": rows}
    with open(os.path.join(RESULTS, "union_dip_results.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("UNION DONE", flush=True)


if __name__ == "__main__":
    union_mode()
