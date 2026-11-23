"""session3.py — stress tests + Bayesian evidence for the evolving-DE hint.

Modes:
  sn_sys     : SN-sample systematics — Dovekie vs Pantheon+ vs DESY5 x {bao_sn_cmb, bao_sn}
  zcut       : Dovekie redshift-cut systematics (z>0.1, z>0.3) on bao_sn_cmb
  seeds      : multi-seed stability of the headline combo (bao_sn_cmb), seeds 11-13
  cov_stress : BAO covariance-inflation robustness, dchi2(alpha) from saved chains (fast)
  evidence   : Savage-Dickey + Laplace Bayes factors from saved chains (fast)

All on public data via src/cosmo.py. Usage: python experiments/session3.py <mode>
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import (Background, bao_chi2, sn_hd_chi2, sn_mB_chi2, cmb_compressed_chi2,
                   metropolis)

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")
os.makedirs(RESULTS, exist_ok=True)

print("loading data ...", flush=True)
bao_all = cosmo.load_desi_bao(which="ALL")[0]
z_dov, mu_dov, invC_dov = cosmo.load_dovekie()
z_pp, mB_pp, invC_pp = cosmo.load_pantheonplus()
z_y5, mu_y5, invC_y5 = cosmo.load_desy5()
print(f"  BAO ALL: {len(bao_all['z'])} rows; SN: dovekie={len(z_dov)} pantheon+={len(z_pp)} desy5={len(z_y5)}", flush=True)

PRIORS = {
    "Om": (0.05, 0.9), "Ob": (0.005, 0.10), "h": (0.55, 0.85), "rd": (120.0, 165.0),
    "w": (-1.6, -0.4), "w0": (-1.8, -0.2), "wa": (-3.0, 3.0),
}
MODELS = {
    "lcdm": ["Om", "Ob", "h", "rd"],
    "wcdm": ["Om", "Ob", "h", "rd", "w"],
    "w0wa": ["Om", "Ob", "h", "rd", "w0", "wa"],
}
X0 = {"Om": 0.31, "Ob": 0.049, "h": 0.67, "rd": 147.0, "w": -1.0, "w0": -0.9, "wa": -0.3}

# SN sample registry: name -> (chi2_fn, label)
SN = {
    "dovekie": (lambda bg, H0: sn_hd_chi2(bg, H0, z_dov, mu_dov, invC_dov), "DES-Dovekie HD"),
    "pantheonplus": (lambda bg, H0: sn_mB_chi2(bg, H0, z_pp, mB_pp, invC_pp), "Pantheon+ mB"),
    "desy5": (lambda bg, H0: sn_hd_chi2(bg, H0, z_y5, mu_y5, invC_y5), "DES-SN5YR HD"),
}


def make_logpost(model, use_bao, sn_name, use_cmb, zmask=None):
    names = MODELS[model]
    sn_fn = SN[sn_name][0] if sn_name is not None else None

    def logpost(x):
        d = dict(zip(names, x))
        for n in names:
            lo, hi = PRIORS[n]
            if not (lo <= d[n] <= hi):
                return -np.inf
        Om, Ob, h, rd = d["Om"], d["Ob"], d["h"], d["rd"]
        w0 = d.get("w0", -1.0)
        wa = d.get("wa", 0.0)
        if "w" in d:
            w0, wa = d["w"], 0.0
        H0 = 100.0 * h
        bg = Background(Om, w0, wa)
        chi2 = 0.0
        if use_bao:
            chi2 += bao_chi2(bg, H0, rd, bao_all)
        if sn_name is not None:
            chi2 += sn_fn(bg, H0)
        if use_cmb:
            chi2 += cmb_compressed_chi2(Om, Ob, h, bg)
        return -0.5 * chi2

    return logpost, names


def sigma_from_dchi2(dchi2, k):
    from scipy.stats import chi2 as chi2dist, norm
    p = chi2dist.sf(abs(dchi2), k)
    z = float(norm.isf(p)) if p < 0.5 else 0.0   # tiny |dchi2| -> no significance
    return z * float(np.sign(dchi2)) if dchi2 != 0 else 0.0, float(p)


def run(model, use_bao, sn_name, use_cmb, nsteps=18000, burn=5000, seed=5, tag="", zmask=None):
    logpost, names = make_logpost(model, use_bao, sn_name, use_cmb, zmask)
    x0 = np.array([X0[n] for n in names])
    step = np.array([0.01 if n != "rd" else 0.5 for n in names])
    t0 = time.time()
    samples, acc = metropolis(logpost, x0, nsteps=nsteps, burn=burn, step=step, seed=seed)
    stride = max(1, len(samples) // 2500)
    sub = samples[::stride]
    lps = np.array([logpost(s) for s in sub])
    chi2_map = float(-2.0 * lps.max())
    summ = {"model": model, "sn": sn_name, "bao": use_bao, "cmb": use_cmb, "tag": tag,
            "acceptance": float(acc), "chi2_map": chi2_map,
            "mean": dict(zip(names, map(float, samples.mean(axis=0)))),
            "std": dict(zip(names, map(float, samples.std(axis=0)))),
            "seconds": round(time.time() - t0, 1)}
    fname = f"chain_s3_{model}_{sn_name or 'nosn'}{'_bao' if use_bao else ''}{'_cmb' if use_cmb else ''}{('_' + tag) if tag else ''}_s{seed}.npy"
    np.save(os.path.join(RESULTS, fname), np.column_stack([sub, lps]))
    print(f"  [{sn_name}{'/bao' if use_bao else ''}{'/cmb' if use_cmb else ''}{'/' + tag if tag else ''}] "
          f"{model}: chi2_map={chi2_map:.2f} acc={acc:.2f} ({summ['seconds']}s) "
          + " ".join(f"{n}={m:.3f}+-{s:.3f}" for n, m, s in
                     zip(names, samples.mean(axis=0), samples.std(axis=0))), flush=True)
    return summ


def compare(res, combo_label, out, fname):
    for a, b, k in [("lcdm", "wcdm", 1), ("lcdm", "w0wa", 2), ("wcdm", "w0wa", 1)]:
        if a not in res or b not in res:
            continue
        d = res[b]["chi2_map"] - res[a]["chi2_map"]
        sig, p = sigma_from_dchi2(d, k)
        res[f"{b}_vs_{a}"] = dict(dchi2=float(d), sigma=sig, p=float(p))
        print(f"  >> {combo_label}: {b} vs {a}: dchi2={d:+.2f} -> {abs(sig):.2f} sigma "
              f"({'prefers ' + b if d < 0 else 'prefers ' + a})", flush=True)
    out.append({"combo": combo_label, "runs": res})
    with open(os.path.join(RESULTS, fname), "w") as f:
        json.dump(out, f, indent=2)


def sn_sys_mode():
    out = []
    for sn in ["dovekie", "pantheonplus", "desy5"]:
        for use_bao, use_cmb in [(True, True), (True, False)]:
            label = f"{sn}_bao{'_cmb' if use_cmb else ''}"
            res = {}
            for model in ["lcdm", "wcdm", "w0wa"]:
                res[model] = run(model, use_bao, sn, use_cmb, nsteps=18000, seed=5)
            compare(res, label, out, "sn_sys_results.json")
    print("SN_SYS DONE", flush=True)


def zcut_mode():
    out = []
    for zmin in [0.0, 0.1, 0.3]:
        global z_dov, mu_dov, invC_dov
        if zmin > 0:
            m = z_dov > zmin
            z_dov_m, mu_dov_m, invC_dov_m = z_dov[m], mu_dov[m], invC_dov[np.ix_(m, m)]
            SN["dovekie"] = (lambda bg, H0: sn_hd_chi2(bg, H0, z_dov_m, mu_dov_m, invC_dov_m),
                             f"Dovekie z>{zmin}")
        label = f"dovekie_zmin{zmin}_bao_cmb"
        res = {}
        for model in ["lcdm", "w0wa"]:
            res[model] = run(model, True, "dovekie", True, nsteps=18000, seed=6, tag=f"zmin{zmin}")
        compare(res, label, out, "zcut_results.json")
    print("ZCUT DONE", flush=True)


def seeds_mode():
    out = []
    for seed in [11, 12, 13]:
        res = {}
        for model in ["lcdm", "w0wa"]:
            res[model] = run(model, True, "dovekie", True, nsteps=24000, seed=seed, tag=f"seed{seed}")
        compare(res, f"bao_sn_cmb_seed{seed}", out, "seeds_results.json")
    print("SEEDS DONE", flush=True)


# ----------------------------------------------------------------------------
# Fast modes operating on saved chains
# ----------------------------------------------------------------------------

def _load_chain(model, combo, seed=1, tag=""):
    p = os.path.join(RESULTS, f"chain_{model}_{combo}{('_' + tag) if tag else ''}_s{seed}.npy")
    d = np.load(p)
    return d[:, :-1], d[:, -1]


def cov_stress_mode():
    """Robustness to BAO covariance inflation: chi2_BAO -> chi2_BAO/alpha^2.
    Profile over saved posterior samples (approximation, labelled as such)."""
    alphas = [1.0, 1.25, 1.5, 2.0, 3.0, 5.0]
    out = []
    for combo in ["bao_sn_cmb", "bao_sn", "bao_cmb", "bao"]:
        use_sn = "sn" in combo
        use_cmb = "cmb" in combo
        entry = {"combo": combo, "alphas": alphas, "curves": {}}
        chi2_bao = {}
        for model in ["lcdm", "w0wa"]:
            s, lp = _load_chain(model, combo)
            names = MODELS[model]
            cb = np.empty(len(s))
            for i, x in enumerate(s):
                d = dict(zip(names, x))
                w0 = d.get("w0", -1.0)
                wa = d.get("wa", 0.0)
                if "w" in d:
                    w0, wa = d["w"], 0.0
                bg = Background(d["Om"], w0, wa)
                cb[i] = bao_chi2(bg, 100.0 * d["h"], d["rd"], bao_all)
            chi2_bao[model] = cb
            chi2_tot = -2.0 * lp
            curve = []
            for a in alphas:
                c2 = chi2_tot - cb + cb / a ** 2
                curve.append(float(c2.min()))
            entry["curves"][model] = curve
        for j, a in enumerate(alphas):
            d = entry["curves"]["w0wa"][j] - entry["curves"]["lcdm"][j]
            sig, p = sigma_from_dchi2(d, 2)
            entry.setdefault("dchi2", {})[f"alpha={a}"] = dict(
                dchi2=float(d), sigma=sig, p=float(p))
            print(f"  [{combo}] alpha={a}: dchi2={d:+.2f} -> {abs(sig):.2f} sigma "
                  f"({'prefers w0wa' if d < 0 else 'prefers lcdm'})", flush=True)
        out.append(entry)
        with open(os.path.join(RESULTS, "cov_stress_results.json"), "w") as f:
            json.dump(out, f, indent=2)
    print("COV_STRESS DONE", flush=True)


def evidence_mode():
    """Bayes factors w0wa/lcdm and wcdm/lcdm:
    (a) Laplace approx: refine MAP with scipy, numerical Hessian, flat-prior volume;
    (b) Savage-Dickey: posterior density of (w0,wa) at (-1,0) via 2D KDE / prior density."""
    from scipy.optimize import minimize
    from scipy.stats import gaussian_kde
    out = []
    for combo in ["bao_sn_cmb", "bao_sn", "bao_cmb", "bao", "sn_cmb"]:
        use_bao, use_sn, use_cmb = "bao" in combo, "sn" in combo, "cmb" in combo
        entry = {"combo": combo, "models": {}}
        maps, hess = {}, {}
        sn_name = "dovekie" if use_sn else None
        for model in ["lcdm", "wcdm", "w0wa"]:
            logpost, names = make_logpost(model, use_bao, sn_name, use_cmb)
            s, lp = _load_chain(model, combo)
            x_map = s[int(np.argmax(lp))]
            # refine MAP
            res = minimize(lambda x: -logpost(x), x_map, method="Nelder-Mead",
                           options={"maxiter": 1500, "xatol": 1e-4, "fatol": 1e-5})
            xm = res.x if res.fun < -lp.max() else x_map
            chi2_map = -2.0 * logpost(xm)
            k = len(names)
            h = np.array([1e-4 if n in ("Om", "h") else 1e-5 if n == "Ob" else
                          0.05 if n == "rd" else 1e-3 for n in names])
            H = np.zeros((k, k))
            for i in range(k):
                for j in range(i, k):
                    xp = xm.copy(); xp[i] += h[i]; xp[j] += h[j]
                    xm_ = xm.copy(); xm_[i] -= h[i]; xm_[j] -= h[j]
                    xpm = xm.copy(); xpm[i] += h[i]; xpm[j] -= h[j]
                    xmp = xm.copy(); xmp[i] -= h[i]; xmp[j] += h[j]
                    def f(x):
                        return -logpost(x)
                    H[i, j] = H[j, i] = (f(xp) - f(xpm) - f(xmp) + f(xm_)) / (4 * h[i] * h[j])
            # symmetrize; H = Hessian of -logpost is PD at the MAP. Clip near-zero
            # eigenvalues (flat directions, e.g. Ob when no CMB) and flag them.
            logZ = float("nan")
            flat = True
            logprior = sum(np.log(PRIORS[n][1] - PRIORS[n][0]) for n in names)
            if np.all(np.isfinite(H)):
                try:
                    H = 0.5 * (H + H.T)
                    ev, V = np.linalg.eigh(H)
                    ev_raw_min = float(ev.min())
                    ev = np.clip(ev, max(ev.max() * 1e-12, 1e-12), None)
                    Sigma = (V / ev) @ V.T          # cov of posterior
                    logdet = float(-np.log(ev).sum())  # logdet(Sigma) = -sum(log ev)
                    flat = bool(ev_raw_min < 1e-6 * float(ev.max()))
                    logZ = float(logpost(xm) + 0.5 * k * np.log(2 * np.pi) + 0.5 * logdet + logprior)
                except np.linalg.LinAlgError:
                    pass
            maps[model] = dict(zip(names, map(float, xm)))
            maps[model]["chi2_map"] = chi2_map
            entry["models"][model] = {"map": maps[model], "logZ_laplace": logZ,
                                      "laplace_flat_direction": flat,
                                      "nsamples": int(len(s))}
            print(f"  [{combo}] {model}: chi2_map={chi2_map:.2f} logZ_laplace={logZ:.1f}"
                  f"{' [FLAT DIR]' if flat else ''}", flush=True)
        # Bayes factors from Laplace
        for a, b in [("lcdm", "wcdm"), ("lcdm", "w0wa")]:
            lz = entry["models"][b]["logZ_laplace"] - entry["models"][a]["logZ_laplace"]
            entry[f"BF_{b}_vs_{a}_laplace"] = float(np.exp(min(lz, 700)))
            entry[f"logBF_{b}_vs_{a}_laplace"] = float(lz)
            print(f"  >> {combo}: logBF({b}/{a}) = {lz:+.2f} (BF={np.exp(min(lz,700)):.3g})", flush=True)
        # Savage-Dickey from chains
        s_w0wa, _ = _load_chain("w0wa", combo)
        s_wcdm, _ = _load_chain("wcdm", combo)
        w0, wa = s_w0wa[:, 4], s_w0wa[:, 5]
        w = s_wcdm[:, 4]
        kde2 = gaussian_kde(np.vstack([w0, wa]))
        post_null = float(kde2(np.array([[-1.0], [0.0]]))[0])
        prior_null = 1.0 / ((PRIORS["w0"][1] - PRIORS["w0"][0]) * (PRIORS["wa"][1] - PRIORS["wa"][0]))
        entry["BF_w0wa_vs_lcdm_savagedickey"] = post_null / prior_null
        kde1 = gaussian_kde(w)
        post_null1 = float(kde1(np.array([-1.0]))[0])
        prior_null1 = 1.0 / (PRIORS["w"][1] - PRIORS["w"][0])
        entry["BF_wcdm_vs_lcdm_savagedickey"] = post_null1 / prior_null1
        print(f"  >> {combo}: SD BF w0wa/lcdm = {entry['BF_w0wa_vs_lcdm_savagedickey']:.3g}, "
              f"wcdm/lcdm = {entry['BF_wcdm_vs_lcdm_savagedickey']:.3g}", flush=True)
        out.append(entry)
        with open(os.path.join(RESULTS, "evidence_results.json"), "w") as f:
            json.dump(out, f, indent=2)
    print("EVIDENCE DONE", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "sn_sys"
    {"sn_sys": sn_sys_mode, "zcut": zcut_mode, "seeds": seeds_mode,
     "cov_stress": cov_stress_mode, "evidence": evidence_mode}[mode]()
