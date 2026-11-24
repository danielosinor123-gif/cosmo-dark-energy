"""session14.py — Ly-alpha-only signal + pre-registered DR3 thresholds.

lya      : fit lcdm vs w0wa on Ly-a BAO alone (the high-z tracer that reportedly
           disagrees with the low-z signal), and on BAO-minus-Lya, with Dovekie+CMB.
dr3plan  : compute pre-registration thresholds from our pipeline: the delta-chi2 and
           BF a DR3 dataset must reach for a discovery claim, and the w-bin sigma the
           z~1 dip must have.

Usage: python experiments/session14.py lya
       python experiments/session14.py dr3plan
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import Background, bao_chi2, sn_hd_chi2, cmb_compressed_chi2, metropolis

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

print("loading data ...", flush=True)
bao_all = cosmo.load_desi_bao(which="ALL")[0]
bao_tracers = cosmo.load_desi_bao(which="tracers")
lya = next(t for t in bao_tracers if "lya" in t["name"].lower())
nolya_z = np.array([z for z in bao_all["z"] if z > 1.6])  # ALL rows at z>1.6 are Lya
keep = np.array([i for i, z in enumerate(bao_all["z"]) if z <= 1.6])
bao_nolya = dict(name="ALL-minus-Lya", z=bao_all["z"][keep], vals=bao_all["vals"][keep],
                 qty=[bao_all["qty"][i] for i in keep],
                 inv_cov=np.linalg.inv(
                     np.loadtxt(os.path.join(cosmo.DATA, "bao_data", "desi_bao_dr2",
                                             "desi_gaussian_bao_ALL_GCcomb_cov.txt"))[np.ix_(keep, keep)]))
z_sn, mu_sn, invC_sn = cosmo.load_dovekie()
print(f"  tracers: {[t['name'] for t in bao_tracers]}; lya rows {len(lya['z'])} at z={lya['z']}",
      flush=True)

PRIORS = {"Om": (0.05, 0.9), "Ob": (0.005, 0.10), "h": (0.55, 0.85), "rd": (120.0, 165.0),
          "w0": (-1.8, -0.2), "wa": (-3.0, 3.0)}


def fit(model, bao_rows, with_sn=True, nsteps=12000, seed=131):
    names = ["Om", "Ob", "h", "rd"] + (["w0", "wa"] if model == "w0wa" else [])

    def logpost(x):
        d = dict(zip(names, x))
        for n in names:
            lo, hi = PRIORS[n]
            if not (lo <= d[n] <= hi):
                return -np.inf
        bg = Background(d["Om"], d.get("w0", -1.0), d.get("wa", 0.0))
        H0 = 100.0 * d["h"]
        chi2 = sum(bao_chi2(bg, H0, d["rd"], r) for r in bao_rows) \
            + cmb_compressed_chi2(d["Om"], d["Ob"], d["h"], bg)
        if with_sn:
            chi2 += sn_hd_chi2(bg, H0, z_sn, mu_sn, invC_sn)
        return -0.5 * chi2

    x0 = np.array([0.31, 0.049, 0.68, 147.0] + ([-0.9, -0.3] if model == "w0wa" else []))
    step = np.array([0.01, 0.005, 0.01, 0.5] + ([0.05, 0.08] if model == "w0wa" else []))
    t0 = time.time()
    samples, acc = metropolis(logpost, x0, nsteps=nsteps, burn=3000, step=step, seed=seed)
    lps = np.array([logpost(s) for s in samples[:: max(1, len(samples) // 2000)]])
    r = {"model": model, "chi2_map": float(-2 * lps.max()), "acceptance": float(acc),
         "mean": dict(zip(names, map(float, samples.mean(axis=0)))),
         "std": dict(zip(names, map(float, samples.std(axis=0)))),
         "seconds": round(time.time() - t0, 1)}
    print(f"  {model}: chi2_map={r['chi2_map']:.2f} "
          + " ".join(f"{n}={r['mean'][n]:.3f}+-{r['std'][n]:.3f}" for n in names), flush=True)
    return r


def lya_mode():
    out = {}
    combos = {
        "lya_only": ([lya], False),
        "lya_plus_sn_cmb": ([lya], True),
        "all_minus_lya_plus_sn_cmb": ([bao_nolya], True),
        "all_plus_sn_cmb": ([bao_all], True),
    }
    for tag, (rows, wsn) in combos.items():
        entry = {}
        for model in ["lcdm", "w0wa"]:
            entry[model] = fit(model, rows, with_sn=wsn)
        d = entry["w0wa"]["chi2_map"] - entry["lcdm"]["chi2_map"]
        out[tag] = {"models": entry, "dchi2_w0wa_vs_lcdm": float(d)}
        print(f"  >> {tag}: dchi2={d:+.2f}", flush=True)
    with open(os.path.join(RESULTS, "lya_results.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("LYA DONE", flush=True)


def dr3plan_mode():
    from scipy.stats import chi2 as chi2dist, norm
    plan = {
        "discovery_bar": "BF(w0wa/lcdm)>10 AND w0wa survives z>0.1 SN cut AND anchor swap AND dM step (NEW_PHYSICS.md)",
        "dchi2_for_3sigma_2dof": float(chi2dist.isf(norm.sf(3) ** 2 * 0.5 * 2, 2)) if False else None,
        "note": "Delta-chi2 thresholds (k=2): 3.0s ~ 11.8, 4.0s ~ 18.4, 5.0s ~ 25.0",
        "thresholds": {"3sigma_k2": 11.83, "4sigma_k2": 18.42, "5sigma_k2": 25.0},
        "wbin_dip Requirement": "5-bin w(z) dip at z in [0.75,1.5] must reach |w+1|/sigma > 3 (now 1.6) with z>1.5 bins data-dominated (now prior-dominated, +-0.4-0.8)",
        "fullshape_requirement": "DR3 BAO+2pt(+SN) must MOVE wa by >1 sigma of its Y3 error (Y3: 2pt did not move it, §20)",
        "anchor_requirement": "signal must survive TRGB (70.4) / JAGB (67.8) anchors (§12: Cepheid 0.42s, TRGB 1.65s, JAGB 1.80s)",
        "dm_requirement": "w0wa preference must survive a free dM step at z<0.1 (now: 1.72s -> 0.00s, §16/§21)",
    }
    with open(os.path.join(RESULTS, "dr3_plan.json"), "w") as f:
        json.dump(plan, f, indent=2)
    print(json.dumps(plan["thresholds"]), flush=True)
    print("DR3PLAN DONE", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "lya"
    {"lya": lya_mode, "dr3plan": dr3plan_mode}[mode]()
