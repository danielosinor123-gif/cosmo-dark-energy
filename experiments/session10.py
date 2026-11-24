"""session10.py — session 10: running-vacuum (RVM) + 5-bin w(z) shape.

rvm <combo>  : fit lcdm / w0wa / RVM on BAO+Dovekie+CMB[(+fs8)].
               RVM: rho_DE(z) = Ode0 + nu*Om*(1+z)^3, Ode0 = 1 - Om*(1+nu) > 0.
               w_eff(z) = -1 + nu*Om*(1+z)^3 / (Ode0 + nu*Om*(1+z)^3).
               (vacuum density runs with matter — a 'channel' between vacuum and matter)
cert5        : 5-bin w(z), edges [0,0.3,0.75,1.5,2.5,inf], with f sigma8 in the joint,
               BAO+Dovekie+CMB+fs8. Resolves the dip shape (N2 crossing test).

Usage: python experiments/session10.py rvm bao_sn_cmb
       python experiments/session10.py cert5
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import Background, BackgroundBinned, lnF_binned, bao_chi2, sn_hd_chi2, \
    cmb_compressed_chi2, fs8_chi2, metropolis
from models import BackgroundFn

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

print("loading data ...", flush=True)
bao_all = cosmo.load_desi_bao(which="ALL")[0]
z_sn, mu_sn, invC_sn = cosmo.load_dovekie()
z_f, fs8v, invC_f = cosmo.load_fs8()

PRIORS = {"Om": (0.05, 0.9), "Ob": (0.005, 0.10), "h": (0.55, 0.85), "rd": (120.0, 165.0),
          "s8": (0.2, 1.4), "w0": (-1.8, -0.2), "wa": (-3.0, 3.0), "nu": (-0.5, 1.0)}
for i in range(5):
    PRIORS[f"w{i}"] = (-3.0, 1.0)
X0 = {"Om": 0.31, "Ob": 0.049, "h": 0.68, "rd": 147.0, "s8": 0.81, "w0": -0.9, "wa": -0.3,
      "nu": 0.05, "w1": -1.0, "w2": -1.0, "w3": -1.0, "w4": -1.0}


def sigma_from_dchi2(dchi2, k):
    from scipy.stats import chi2 as chi2dist, norm
    p = chi2dist.sf(abs(dchi2), k)
    z = float(norm.isf(p)) if p < 0.5 else 0.0
    return z * float(np.sign(dchi2)) if dchi2 != 0 else 0.0, float(p)


def make_logpost(model, combo):
    use_fs8 = "fs8" in combo
    if model == "lcdm":
        names = ["Om", "Ob", "h", "rd"] + (["s8"] if use_fs8 else [])
    elif model == "w0wa":
        names = ["Om", "Ob", "h", "rd", "w0", "wa"] + (["s8"] if use_fs8 else [])
    else:  # rvm
        names = ["Om", "Ob", "h", "rd", "nu"] + (["s8"] if use_fs8 else [])

    def logpost(x):
        d = dict(zip(names, x))
        for n in names:
            lo, hi = PRIORS[n]
            if not (lo <= d[n] <= hi):
                return -np.inf
        Om, Ob, h, rd = d["Om"], d["Ob"], d["h"], d["rd"]
        if model == "lcdm":
            bg = Background(Om); lnF = None
        elif model == "w0wa":
            bg = Background(Om, d["w0"], d["wa"]); lnF = None
        else:
            nu = d["nu"]
            Ode0 = 1.0 - Om * (1.0 + nu)
            if Ode0 <= 1e-3:
                return -np.inf
            lnF = lambda z: np.log((Ode0 + nu * Om * (1.0 + np.asarray(z, float)) ** 3) / Ode0)
            bg = BackgroundFn(Om, lnF)
        H0 = 100.0 * h
        chi2 = bao_chi2(bg, H0, rd, bao_all) + sn_hd_chi2(bg, H0, z_sn, mu_sn, invC_sn) \
            + cmb_compressed_chi2(Om, Ob, h, bg)
        if use_fs8:
            chi2 += fs8_chi2(Om, d["s8"], z_f, fs8v, invC_f, lnF=lnF)
        return -0.5 * chi2

    return logpost, names


def run_mcmc(model, combo, nsteps=16000, seed=101):
    logpost, names = make_logpost(model, combo)
    x0 = np.array([X0[n] for n in names])
    step = np.array([{"Om": 0.01, "Ob": 0.005, "h": 0.01, "rd": 0.5, "s8": 0.02,
                      "w0": 0.05, "wa": 0.08, "nu": 0.02}[n] for n in names])
    t0 = time.time()
    samples, acc = metropolis(logpost, x0, nsteps=nsteps, burn=max(4000, nsteps // 4),
                              step=step, seed=seed)
    stride = max(1, len(samples) // 2500)
    sub = samples[::stride]
    lps = np.array([logpost(s) for s in sub])
    chi2_map = float(-2.0 * lps.max())
    r = {"model": model, "combo": combo, "chi2_map": chi2_map, "acceptance": float(acc),
         "mean": dict(zip(names, map(float, samples.mean(axis=0)))),
         "std": dict(zip(names, map(float, samples.std(axis=0)))),
         "seconds": round(time.time() - t0, 1)}
    np.save(os.path.join(RESULTS, f"chain_s10_{model}_{combo}_s{seed}.npy"),
            np.column_stack([sub, lps]))
    print(f"  [{combo}] {model}: chi2_map={chi2_map:.2f} acc={acc:.2f} "
          + " ".join(f"{n}={m:.3f}+-{s:.3f}" for n, m, s in
                     zip(names, samples.mean(axis=0), samples.std(axis=0)))
          + f" ({r['seconds']}s)", flush=True)
    return r


def rvm_mode():
    combo = sys.argv[2] if len(sys.argv) > 2 else "bao_sn_cmb"
    try:
        with open(os.path.join(RESULTS, "rvm_results.json")) as f:
            out = [e for e in json.load(f) if e["combo"] != combo]
        prev = next(e for e in json.load(open(os.path.join(RESULTS, "rvm_results.json")))
                    if e["combo"] == combo)
    except Exception:
        out, prev = [], None
    entry = prev if prev else {"combo": combo, "models": {}}
    for model in ["lcdm", "w0wa", "rvm"]:
        if model in entry.get("models", {}):
            print(f"  {model} already done", flush=True)
            continue
        entry.setdefault("models", {})[model] = run_mcmc(model, combo)
        with open(os.path.join(RESULTS, "rvm_results.json"), "w") as f:
            json.dump(out + [entry], f, indent=2)
    base = entry["models"]["lcdm"]["chi2_map"]
    for model in ["w0wa", "rvm"]:
        d = entry["models"][model]["chi2_map"] - base
        k = 2 if model == "w0wa" else 1
        sig, p = sigma_from_dchi2(d, k)
        entry[f"{model}_vs_lcdm"] = {"dchi2": float(d), "sigma": float(sig), "p": float(p)}
        print(f"  >> [{combo}] {model} vs lcdm: dchi2={d:+.2f} ({abs(sig):.2f}s, k={k})",
              flush=True)
    with open(os.path.join(RESULTS, "rvm_results.json"), "w") as f:
        json.dump(out + [entry], f, indent=2)
    print("RVM DONE", flush=True)


EDGES5 = [0.0, 0.3, 0.75, 1.5, 2.5, np.inf]


def cert5_mode():
    names = ["Om", "Ob", "h", "rd", "s8"] + [f"w{i}" for i in range(5)]

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
        chi2 = bao_chi2(bg, H0, rd, bao_all) + sn_hd_chi2(bg, H0, z_sn, mu_sn, invC_sn) \
            + cmb_compressed_chi2(Om, Ob, h, bg) \
            + fs8_chi2(Om, s8, z_f, fs8v, invC_f, lnF=lnF)
        return -0.5 * chi2

    x0 = np.array([X0[n] for n in names])
    step = np.array([0.01 if n == "Om" else 0.005 if n == "Ob" else 0.01 if n == "h"
                     else 0.5 if n == "rd" else 0.02 if n == "s8" else 0.06 for n in names])
    t0 = time.time()
    samples, acc = metropolis(logpost, x0, nsteps=24000, burn=6000, step=step, seed=111)
    stride = max(1, len(samples) // 3000)
    sub = samples[::stride]
    lps = np.array([logpost(s) for s in sub])
    wi = samples[:, 5:10]
    entry = {"edges": EDGES5, "chi2_map": float(-2.0 * lps.max()), "acceptance": float(acc),
             "w_bins_mean": wi.mean(axis=0).tolist(), "w_bins_std": wi.std(axis=0).tolist(),
             "mean": dict(zip(names, map(float, samples.mean(axis=0)))),
             "std": dict(zip(names, map(float, samples.std(axis=0)))),
             "seconds": round(time.time() - t0, 1)}
    np.save(os.path.join(RESULTS, "chain_s10_cert5.npy"), np.column_stack([sub, lps]))
    edges_lbl = ["0-0.3", "0.3-0.75", "0.75-1.5", "1.5-2.5", ">2.5"]
    print("  cert5 w bins: " + " | ".join(
        f"{lbl}: {m:+.3f}+-{s:.3f}" for lbl, m, s in zip(edges_lbl, wi.mean(axis=0), wi.std(axis=0)))
        + f"  chi2_map={entry['chi2_map']:.2f} acc={acc:.2f} ({entry['seconds']}s)", flush=True)
    with open(os.path.join(RESULTS, "cert5_results.json"), "w") as f:
        json.dump(entry, f, indent=2)
    print("CERT5 DONE", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "rvm"
    {"rvm": rvm_mode, "cert5": cert5_mode}[mode]()
