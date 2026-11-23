"""session4.py — (A) low-z anchor test, (B) growth-inclusive w(z) certificates.

anchor : does the w0wa hint survive different local H0 anchors?
         H0 priors: none / Cepheid 73.04+-1.04 / CCHP TRGB 70.4+-1.9 / CCHP JAGB
         67.8+-2.7 / Planck LCDM 67.36+-0.54, each with Dovekie all-z and z>0.1.
cert   : binned w(z) [0,0.75,1.5,inf] certificates with f sigma8 in the joint
         (geometry vs growth tension per redshift bin).

Usage: python experiments/session4.py anchor [zcut]   (zcut in {0.0, 0.1})
       python experiments/session4.py cert
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import (Background, BackgroundBinned, lnF_binned, bao_chi2, sn_hd_chi2,
                   cmb_compressed_chi2, fs8_chi2, metropolis)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")
os.makedirs(RESULTS, exist_ok=True)

print("loading data ...", flush=True)
bao_all = cosmo.load_desi_bao(which="ALL")[0]
z_sn, mu_sn, invC_sn = cosmo.load_dovekie()
z_f, fs8v, invC_f = cosmo.load_fs8()
print(f"  BAO {len(bao_all['z'])} rows, SN {len(z_sn)}, fs8 {len(z_f)}", flush=True)

PRIORS = {
    "Om": (0.05, 0.9), "Ob": (0.005, 0.10), "h": (0.55, 0.85), "rd": (120.0, 165.0),
    "s8": (0.2, 1.4),
    "w": (-1.6, -0.4), "w0": (-1.8, -0.2), "wa": (-3.0, 3.0),
    "w1": (-3.0, 1.0), "w2": (-3.0, 1.0),
}
MODELS = {
    "lcdm": ["Om", "Ob", "h", "rd"],
    "w0wa": ["Om", "Ob", "h", "rd", "w0", "wa"],
    "binned": ["Om", "Ob", "h", "rd", "s8", "w0", "w1", "w2"],
}
X0 = {"Om": 0.31, "Ob": 0.049, "h": 0.67, "rd": 147.0, "s8": 0.81,
      "w": -1.0, "w0": -0.9, "wa": -0.3, "w1": -0.8, "w2": -0.6}

H0_PRIORS = {   # (H0, sigma) — public 2025-26 values
    "none": None,
    "cepheid_shoes": (73.04, 1.04),     # Riess+2022
    "trgb_cchp": (70.4, 1.9),           # Freedman+2025 (2408.06153): 70.39 +- 1.22 stat +- 1.33 sys +- 0.70 SN
    "jagb_cchp": (67.8, 2.7),           # Lee+2024/Freedman+2025: 67.80 +- 2.17 stat +- 1.64 sys
    "planck_lcdm": (67.36, 0.54),       # Planck 2018
}

EDGES = [0.0, 0.75, 1.5, np.inf]


def sigma_from_dchi2(dchi2, k):
    from scipy.stats import chi2 as chi2dist, norm
    p = chi2dist.sf(abs(dchi2), k)
    z = float(norm.isf(p)) if p < 0.5 else 0.0
    return z * float(np.sign(dchi2)) if dchi2 != 0 else 0.0, float(p)


def make_logpost(model, zmask, h0prior, use_fs8):
    names = MODELS[model]
    if zmask is not None:
        m = z_sn > zmask
        z_s, mu_s, invC_s = z_sn[m], mu_sn[m], invC_sn[np.ix_(m, m)]
    else:
        z_s, mu_s, invC_s = z_sn, mu_sn, invC_sn

    def logpost(x):
        d = dict(zip(names, x))
        for n in names:
            lo, hi = PRIORS[n]
            if not (lo <= d[n] <= hi):
                return -np.inf
        Om, Ob, h, rd = d["Om"], d["Ob"], d["h"], d["rd"]
        if model == "binned":
            wb = [d["w0"], d["w1"], d["w2"]]
            bg = BackgroundBinned(Om, wb, EDGES)
            lnF = lambda zz: lnF_binned(zz, EDGES, wb)
            w0, wa = None, None
        else:
            w0 = d.get("w0", -1.0)
            wa = d.get("wa", 0.0)
            if "w" in d:
                w0, wa = d["w"], 0.0
            bg = Background(Om, w0, wa)
            lnF = None
        H0 = 100.0 * h
        chi2 = bao_chi2(bg, H0, rd, bao_all) + sn_hd_chi2(bg, H0, z_s, mu_s, invC_s) \
            + cmb_compressed_chi2(Om, Ob, h, bg)
        if use_fs8:
            chi2 += fs8_chi2(Om, d["s8"], z_f, fs8v, invC_f, lnF=lnF)
        if h0prior is not None:
            chi2 += ((H0 - h0prior[0]) / h0prior[1]) ** 2
        return -0.5 * chi2

    return logpost, names


def run(model, tag, zmask, h0prior, use_fs8, nsteps, seed):
    logpost, names = make_logpost(model, zmask, h0prior, use_fs8)
    x0 = np.array([X0[n] for n in names])
    step = np.array([0.01 if n not in ("rd", "s8") else 0.5 if n == "rd" else 0.02 for n in names])
    t0 = time.time()
    samples, acc = metropolis(logpost, x0, nsteps=nsteps, burn=max(4000, nsteps // 4),
                              step=step, seed=seed)
    stride = max(1, len(samples) // 2500)
    sub = samples[::stride]
    lps = np.array([logpost(s) for s in sub])
    chi2_map = float(-2.0 * lps.max())
    summ = {"model": model, "tag": tag, "acceptance": float(acc), "chi2_map": chi2_map,
            "mean": dict(zip(names, map(float, samples.mean(axis=0)))),
            "std": dict(zip(names, map(float, samples.std(axis=0)))),
            "seconds": round(time.time() - t0, 1)}
    np.save(os.path.join(RESULTS, f"chain_s4_{model}_{tag}_s{seed}.npy"),
            np.column_stack([sub, lps]))
    print(f"  [{tag}] {model}: chi2_map={chi2_map:.2f} acc={acc:.2f} ({summ['seconds']}s) "
          + " ".join(f"{n}={m:.3f}+-{s:.3f}" for n, m, s in
                     zip(names, samples.mean(axis=0), samples.std(axis=0))), flush=True)
    return summ


def anchor_mode():
    zcut = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
    zmask = zcut if zcut > 0 else None
    try:
        with open(os.path.join(RESULTS, "anchor_results.json")) as f:
            out = [e for e in json.load(f) if e["zcut"] != zcut]
    except Exception:
        out = []
    for pname, hp in H0_PRIORS.items():
        tag = f"zcut{zcut}_{pname}"
        res = {}
        for model in ["lcdm", "w0wa"]:
            res[model] = run(model, tag, zmask, hp, False, nsteps=16000, seed=21)
        d = res["w0wa"]["chi2_map"] - res["lcdm"]["chi2_map"]
        sig, p = sigma_from_dchi2(d, 2)
        res["w0wa_vs_lcdm"] = dict(dchi2=float(d), sigma=sig, p=float(p))
        print(f"  >> {tag}: dchi2={d:+.2f} -> {abs(sig):.2f} sigma "
              f"({'prefers w0wa' if d < 0 else 'prefers lcdm'})", flush=True)
        out.append({"zcut": zcut, "h0_prior": pname, "runs": res})
        with open(os.path.join(RESULTS, "anchor_results.json"), "w") as f:
            json.dump(out, f, indent=2)
    print("ANCHOR DONE", flush=True)


def cert_mode():
    combos = [
        ("joint_all", False, False),        # BAO+SN+CMB
        ("joint_fs8", False, True),         # BAO+SN+CMB+f8
        ("bao_cmb_fs8", True, True),        # BAO+CMB+f8 (no SN)
        ("bao_sn_fs8", False, True),        # placeholder replaced below
    ]
    # (tag, bao_only?, use_fs8) — implement via explicit flags
    combos = [
        ("joint_all", dict(fs8=False)),
        ("joint_fs8", dict(fs8=True)),
        ("bao_cmb_fs8", dict(fs8=True, drop_sn=True)),
        ("bao_sn_fs8", dict(fs8=True, drop_cmb=True)),
        ("fs8_sn_cmb", dict(fs8=True, drop_bao=True)),
    ]
    try:
        with open(os.path.join(RESULTS, "cert_fs8_results.json")) as f:
            done = {e["tag"] for e in json.load(f)}
    except Exception:
        done = set()
    out = []
    try:
        with open(os.path.join(RESULTS, "cert_fs8_results.json")) as f:
            out = json.load(f)
    except Exception:
        pass
    for tag, cfg in combos:
        if tag in done:
            continue
        use_fs8 = cfg.get("fs8", False)

        def logpost_factory():
            names = MODELS["binned"]

            def logpost(x):
                d = dict(zip(names, x))
                for n in names:
                    lo, hi = PRIORS[n]
                    if not (lo <= d[n] <= hi):
                        return -np.inf
                Om, Ob, h, rd, s8 = d["Om"], d["Ob"], d["h"], d["rd"], d["s8"]
                wb = [d["w0"], d["w1"], d["w2"]]
                bg = BackgroundBinned(Om, wb, EDGES)
                lnF = lambda zz: lnF_binned(zz, EDGES, wb)
                H0 = 100.0 * h
                chi2 = 0.0
                if not cfg.get("drop_bao", False):
                    chi2 += bao_chi2(bg, H0, rd, bao_all)
                if not cfg.get("drop_sn", False):
                    chi2 += sn_hd_chi2(bg, H0, z_sn, mu_sn, invC_sn)
                if not cfg.get("drop_cmb", False):
                    chi2 += cmb_compressed_chi2(Om, Ob, h, bg)
                if use_fs8:
                    chi2 += fs8_chi2(Om, s8, z_f, fs8v, invC_f, lnF=lnF)
                return -0.5 * chi2
            return logpost, names

        logpost, names = logpost_factory()
        x0 = np.array([X0[n] for n in names])
        step = np.array([0.01 if n not in ("rd", "s8") else 0.5 if n == "rd" else 0.02
                         for n in names])
        t0 = time.time()
        samples, acc = metropolis(logpost, x0, nsteps=20000, burn=5000, step=step, seed=31)
        stride = max(1, len(samples) // 3000)
        sub = samples[::stride]
        lps = np.array([logpost(s) for s in sub])
        chi2_map = float(-2.0 * lps.max())
        wi = samples[:, 5:8]
        entry = {"tag": tag, "chi2_map": chi2_map, "acceptance": float(acc),
                 "w_bins_mean": wi.mean(axis=0).tolist(), "w_bins_std": wi.std(axis=0).tolist(),
                 "mean": dict(zip(names, map(float, samples.mean(axis=0)))),
                 "std": dict(zip(names, map(float, samples.std(axis=0)))),
                 "seconds": round(time.time() - t0, 1)}
        np.save(os.path.join(RESULTS, f"chain_s4_binned_{tag}.npy"),
                np.column_stack([sub, lps]))
        print(f"  [{tag}] chi2_map={chi2_map:.2f} acc={acc:.2f} w_bins="
              + " ".join(f"{m:.2f}+-{s:.2f}" for m, s in zip(wi.mean(axis=0), wi.std(axis=0)))
              + f" ({entry['seconds']}s)", flush=True)
        out = [e for e in out if e["tag"] != tag] + [entry]
        with open(os.path.join(RESULTS, "cert_fs8_results.json"), "w") as f:
            json.dump(out, f, indent=2)
    # tension vs joint_fs8
    jf = next(e for e in out if e["tag"] == "joint_fs8")
    for e in out:
        if e["tag"] == "joint_fs8":
            continue
        e["tension_vs_joint_fs8"] = [(a - b) / np.hypot(sa, sb) for a, b, sa, sb in zip(
            e["w_bins_mean"], jf["w_bins_mean"], e["w_bins_std"], jf["w_bins_std"])]
        print(f"  tension {e['tag']} vs joint_fs8: "
              + " ".join(f"{t:+.2f}" for t in e["tension_vs_joint_fs8"]), flush=True)
    with open(os.path.join(RESULTS, "cert_fs8_results.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("CERT DONE", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "anchor"
    {"anchor": anchor_mode, "cert": cert_mode}[mode]()
