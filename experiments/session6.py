"""session6.py — NEW PHYSICS: fit candidate dark-energy models against the same
public data and let the data judge them.

Models (all flat, non-interacting backgrounds via lnF(z), see src/models.py):
  lcdm    : w = -1
  w0wa    : CPL w(z) = w0 + wa z/(1+z)              (baseline evolving-DE)
  ldpt    : Late Dark-sector Phase Transition — rho_DE(z) = rho_DE0 * g(z) with a
            logistic step of depth f at z=zc, width d (late vacuum metamorphosis)
  quintom : w(z) = wb + A exp(-((z-zc)/d)^2) — localized phantom-divide crossing
            (two-field dark sector)

Modes:
  fit <combo>     combo in {bao_sn_cmb, bao_sn_cmb_fs8}: MCMC all 4 models,
                  dchi2 vs lcdm, sigma, Laplace Bayes factor, H0 posterior.
                  -> results/session6_results.json, chain_s6_*.npy
  h0map <model>   fix h on a grid, MAP-fit (Nelder-Mead), chi2(H0) profile
                  (full SN sample) -> results/h0map_s6_<model>.json
  predict <model> <combo>: best-fit w(z), Omega_DE(z), fs8 prediction vs data
                  -> results/predict_s6_<model>_<combo>.json

Usage: python experiments/session6.py fit bao_sn_cmb
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import Background, bao_chi2, sn_hd_chi2, cmb_compressed_chi2, fs8_chi2, metropolis
from models import BackgroundFn, lnF_ldpt, make_lnF_quintom, w_quintom, w_eff_ldpt

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")
os.makedirs(RESULTS, exist_ok=True)

print("loading data ...", flush=True)
bao_all = cosmo.load_desi_bao(which="ALL")[0]
z_sn, mu_sn, invC_sn = cosmo.load_dovekie()
z_f, fs8v, invC_f = cosmo.load_fs8()
print(f"  BAO {len(bao_all['z'])} rows; SN {len(z_sn)}; fs8 {len(z_f)}", flush=True)

PRIORS = {
    "Om": (0.05, 0.9), "Ob": (0.005, 0.10), "h": (0.55, 0.85), "rd": (120.0, 165.0),
    "s8": (0.2, 1.4), "w0": (-1.8, -0.2), "wa": (-3.0, 3.0),
    "f": (0.0, 0.9), "zc": (0.2, 4.0), "d": (0.03, 1.5),
    "wb": (-1.3, -0.6), "A": (-1.0, 0.5), "zq": (0.2, 3.0), "dq": (0.05, 2.0),
}
X0 = {"Om": 0.31, "Ob": 0.049, "h": 0.68, "rd": 147.0, "s8": 0.81,
      "w0": -0.9, "wa": -0.3, "f": 0.3, "zc": 1.2, "d": 0.3,
      "wb": -0.95, "A": -0.3, "zq": 1.0, "dq": 0.6}

MODEL_PARAMS = {
    "lcdm": ["Om", "Ob", "h", "rd"],
    "w0wa": ["Om", "Ob", "h", "rd", "w0", "wa"],
    "ldpt": ["Om", "Ob", "h", "rd", "f", "zc", "d"],
    "quintom": ["Om", "Ob", "h", "rd", "wb", "A", "zq", "dq"],
}


def sigma_from_dchi2(dchi2, k):
    from scipy.stats import chi2 as chi2dist, norm
    p = chi2dist.sf(abs(dchi2), k)
    z = float(norm.isf(p)) if p < 0.5 else 0.0
    return z * float(np.sign(dchi2)) if dchi2 != 0 else 0.0, float(p)


def make_logpost(model, combo, h_fixed=None):
    use_fs8 = "fs8" in combo
    names = list(MODEL_PARAMS[model]) + (["s8"] if use_fs8 else [])
    if h_fixed is not None:
        names = [n for n in names if n != "h"]

    def logpost(x):
        d = dict(zip(names, x))
        for n in names:
            lo, hi = PRIORS[n]
            if not (lo <= d[n] <= hi):
                return -np.inf
        Om, Ob, rd = d["Om"], d["Ob"], d["rd"]
        h = h_fixed if h_fixed is not None else d["h"]
        if model == "lcdm":
            bg = Background(Om)
            lnF = None
        elif model == "w0wa":
            bg = Background(Om, d["w0"], d["wa"])
            lnF = None
        elif model == "ldpt":
            lnF = lambda z: lnF_ldpt(z, d["f"], d["zc"], d["d"])
            bg = BackgroundFn(Om, lnF)
        else:
            lnF = make_lnF_quintom(d["wb"], d["A"], d["zq"], d["dq"])
            bg = BackgroundFn(Om, lnF)
        H0 = 100.0 * h
        chi2 = bao_chi2(bg, H0, rd, bao_all) + sn_hd_chi2(bg, H0, z_sn, mu_sn, invC_sn) \
            + cmb_compressed_chi2(Om, Ob, h, bg)
        if use_fs8:
            chi2 += fs8_chi2(Om, d["s8"], z_f, fs8v, invC_f, lnF=lnF)
        return -0.5 * chi2

    return logpost, names


def run_mcmc(model, combo, nsteps=14000, seed=61):
    logpost, names = make_logpost(model, combo)
    x0 = np.array([X0[n] for n in names])
    step = np.array([{"Om": 0.01, "Ob": 0.005, "h": 0.01, "rd": 0.5, "s8": 0.02,
                      "w0": 0.05, "wa": 0.08, "f": 0.03, "zc": 0.1, "d": 0.05,
                      "wb": 0.03, "A": 0.05, "zq": 0.08, "dq": 0.08}[n] for n in names])
    t0 = time.time()
    samples, acc = metropolis(logpost, x0, nsteps=nsteps, burn=max(3500, nsteps // 4),
                              step=step, seed=seed)
    stride = max(1, len(samples) // 2500)
    sub = samples[::stride]
    lps = np.array([logpost(s) for s in sub])
    x_map = sub[int(np.argmax(lps))]
    chi2_map = float(-2.0 * lps.max())
    summ = {"model": model, "combo": combo, "chi2_map": chi2_map, "acceptance": float(acc),
            "map": dict(zip(names, map(float, x_map))),
            "mean": dict(zip(names, map(float, samples.mean(axis=0)))),
            "std": dict(zip(names, map(float, samples.std(axis=0)))),
            "seconds": round(time.time() - t0, 1)}
    np.save(os.path.join(RESULTS, f"chain_s6_{model}_{combo}_s{seed}.npy"),
            np.column_stack([sub, lps]))
    print(f"  [{combo}] {model}: chi2_map={chi2_map:.2f} acc={acc:.2f} "
          + " ".join(f"{n}={m:.3f}+-{s:.3f}" for n, m, s in
                     zip(names, samples.mean(axis=0), samples.std(axis=0)))
          + f" ({summ['seconds']}s)", flush=True)
    return summ, logpost, names


def laplace_logZ(logpost, names, x_map):
    from scipy.optimize import minimize
    res = minimize(lambda x: -logpost(x), x_map, method="Nelder-Mead",
                   options={"maxiter": 2000, "xatol": 1e-4, "fatol": 1e-6})
    xm = res.x if res.fun < -logpost(x_map) else x_map
    lp_map = logpost(xm)
    k = len(names)
    h = np.array([1e-4 if n in ("Om", "h") else 1e-5 if n == "Ob" else
                  0.05 if n == "rd" else 1e-3 for n in names])
    H = np.zeros((k, k))
    f = lambda x: -logpost(x)
    for i in range(k):
        for j in range(i, k):
            xp = xm.copy(); xp[i] += h[i]; xp[j] += h[j]
            xm_ = xm.copy(); xm_[i] -= h[i]; xm_[j] -= h[j]
            xpm = xm.copy(); xpm[i] += h[i]; xpm[j] -= h[j]
            xmp = xm.copy(); xmp[i] -= h[i]; xmp[j] += h[j]
            H[i, j] = H[j, i] = (f(xp) - f(xpm) - f(xmp) + f(xm_)) / (4 * h[i] * h[j])
    logZ = float("nan")
    flat = True
    if np.all(np.isfinite(H)):
        try:
            H = 0.5 * (H + H.T)
            ev, V = np.linalg.eigh(H)
            flat = bool(ev.min() < 1e-6 * ev.max())
            ev = np.clip(ev, max(ev.max() * 1e-12, 1e-12), None)
            logdet = float(-np.log(ev).sum())
            logprior = sum(np.log(PRIORS[n][1] - PRIORS[n][0]) for n in names)
            logZ = float(lp_map + 0.5 * k * np.log(2 * np.pi) + 0.5 * logdet + logprior)
        except np.linalg.LinAlgError:
            pass
    return logZ, flat, dict(zip(names, map(float, xm))), float(-2.0 * lp_map)


def fit_mode():
    combo = sys.argv[2] if len(sys.argv) > 2 else "bao_sn_cmb"
    try:
        with open(os.path.join(RESULTS, "session6_results.json")) as f:
            all_entries = json.load(f)
        out = [e for e in all_entries if e["combo"] != combo]
        prev = next((e for e in all_entries if e["combo"] == combo), None)
    except Exception:
        out, prev = [], None
    entry = prev if prev else {"combo": combo, "models": {}}
    entry.setdefault("models", {})
    logZ = {}
    for model in ["lcdm", "w0wa", "ldpt", "quintom"]:
        if model in entry["models"] and "logZ_laplace" in entry["models"][model]:
            logZ[model] = entry["models"][model]["logZ_laplace"]
            print(f"  [{combo}] {model}: already done, skipping", flush=True)
            continue
        summ, logpost, names = run_mcmc(model, combo)
        lz, flat, mapd, chi2r = laplace_logZ(logpost, names,
                                             np.array([summ["map"][n] for n in names]))
        summ["logZ_laplace"] = lz
        summ["laplace_flat"] = flat
        summ["chi2_map_refined"] = chi2r
        logZ[model] = lz
        entry["models"][model] = summ
        print(f"  >> [{combo}] {model}: logZ={lz:.2f}{' [FLAT]' if flat else ''} "
              f"chi2_refined={chi2r:.2f}", flush=True)
        with open(os.path.join(RESULTS, "session6_results.json"), "w") as f:
            json.dump(out + [entry], f, indent=2)
    base = entry["models"]["lcdm"]["chi2_map_refined"]
    for model in ["w0wa", "ldpt", "quintom"]:
        d = entry["models"][model]["chi2_map_refined"] - base
        k = len(MODEL_PARAMS[model]) - len(MODEL_PARAMS["lcdm"])
        sig, p = sigma_from_dchi2(d, k)
        dlz = logZ[model] - logZ["lcdm"]
        bf = float(np.exp(min(dlz, 700))) if np.isfinite(dlz) else None
        entry[f"{model}_vs_lcdm"] = {"dchi2": float(d), "sigma": float(sig), "p": float(p),
                                     "BF_laplace": bf}
        print(f"  >> [{combo}] {model} vs lcdm: dchi2={d:+.2f} ({abs(sig):.2f}s, k={k}) "
              f"BF_laplace={'n/a' if bf is None else format(bf, '.3g')}", flush=True)
    out.append(entry)
    with open(os.path.join(RESULTS, "session6_results.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("FIT DONE", flush=True)


def h0map_mode():
    model = sys.argv[2] if len(sys.argv) > 2 else "ldpt"
    from scipy.optimize import minimize
    hgrid = np.arange(0.660, 0.741, 0.01)
    out = []
    for h in hgrid:
        row = {"h": float(h), "H0": float(100 * h)}
        for m in [model, "lcdm"]:
            logpost, names = make_logpost(m, "bao_sn_cmb", h_fixed=float(h))
            starts = []
            base = [X0[n] for n in names if n != "h"]
            starts.append(base)
            if m == "ldpt":
                starts.append([0.31, 0.049, 143.5, 0.5, 1.5, 0.3])
            best = None
            for x0 in starts:
                res = minimize(lambda x: -logpost(x), np.array(x0, float), method="Nelder-Mead",
                               options={"maxiter": 5000, "xatol": 1e-5, "fatol": 1e-7})
                if res.fun < (best[0] if best else np.inf):
                    best = (res.fun, res.x)
            row[m] = {"chi2_map": float(-2.0 * (-best[0])),
                      "map": dict(zip(names, map(float, best[1])))}
            print(f"  h={h:.3f} {m}: chi2={row[m]['chi2_map']:.2f} "
                  + " ".join(f"{n}={v:.3f}" for n, v in zip(names, best[1])), flush=True)
        out.append(row)
        with open(os.path.join(RESULTS, f"h0map_s6_{model}.json"), "w") as f:
            json.dump(out, f, indent=2)
    print("H0MAP6 DONE", flush=True)


def predict_mode():
    model = sys.argv[2] if len(sys.argv) > 2 else "ldpt"
    combo = sys.argv[3] if len(sys.argv) > 3 else "bao_sn_cmb"
    with open(os.path.join(RESULTS, "session6_results.json")) as f:
        entry = next(e for e in json.load(f) if e["combo"] == combo)
    mp = entry["models"][model]["map"]
    Om = mp["Om"]
    zg = np.linspace(0.0, 3.0, 31)
    if model == "ldpt":
        f, zc, d = mp["f"], mp["zc"], mp["d"]
        lnF = lambda z: lnF_ldpt(z, f, zc, d)
        w = w_eff_ldpt(zg, f, zc, d)
        kind = "ldpt"
    elif model == "quintom":
        wb, A, zq, dq = mp["wb"], mp["A"], mp["zq"], mp["dq"]
        lnF = make_lnF_quintom(wb, A, zq, dq)
        w = w_quintom(zg, wb, A, zq, dq)
        kind = "quintom"
    elif model == "w0wa":
        lnF = None
        w = mp["w0"] + mp["wa"] * zg / (1 + zg)
        kind = "w0wa"
    else:
        lnF = None
        w = -np.ones_like(zg)
        kind = "lcdm"
    bg = BackgroundFn(Om, lnF) if lnF is not None else Background(Om, mp.get("w0", -1), mp.get("wa", 0))
    Ode = 1.0 - Om
    E2 = Om * (1 + zg) ** 3 + Ode * np.exp(lnF(zg)) if lnF is not None else \
        Om * (1 + zg) ** 3 + Ode
    Omde = Ode * np.exp(lnF(zg)) / E2 if lnF is not None else Ode / E2
    H0 = 100.0 * mp["h"]
    out = {"model": model, "combo": combo, "map": mp,
           "z": zg.tolist(), "w_eff": np.asarray(w, float).tolist(),
           "Omega_DE": np.asarray(Omde, float).tolist(), "H0": float(H0)}
    if "fs8" in combo and "s8" in mp:
        from cosmo import growth, fs8_model
        s8 = mp["s8"]
        pred = fs8_model(Om, s8, z_f, lnF=lnF) if lnF is not None else \
            fs8_model(Om, s8, z_f, mp.get("w0", -1), mp.get("wa", 0))
        out["fs8"] = {"z": z_f.tolist(), "data": fs8v.tolist(), "model": pred.tolist()}
        r = fs8v - pred
        out["fs8_chi2"] = float(r @ invC_f @ r)
    with open(os.path.join(RESULTS, f"predict_s6_{model}_{combo}.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"  [{model}/{combo}] H0={H0:.2f} w(z=0,0.5,1,1.5,2,3)="
          + ",".join(f"{wi:+.3f}" for wi in
                     [float(np.interp(zz, zg, w)) for zz in [0, 0.5, 1, 1.5, 2, 3]])
          + f"  fs8_chi2={out.get('fs8_chi2')}", flush=True)
    print("PREDICT DONE", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "fit"
    {"fit": fit_mode, "h0map": h0map_mode, "predict": predict_mode}[mode]()
