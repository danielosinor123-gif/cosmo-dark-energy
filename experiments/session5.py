"""session5.py — NEW quantitative discoveries (not just stress tests):

step  : free low-z magnitude step dM for z<z_step — how big a calibration step
        absorbs the whole w0wa preference?  (lcdm vs w0wa, with/without dM,
        BAO+Dovekie+CMB, z_step in {0.05,0.1,0.2})
xres  : cross-sample low-z residual certificate — do Dovekie / Pantheon+ / DESY5
        disagree with each other in the z<0.1 bin where the signal lives?
        (geometry-only LCDM from BAO+CMB, per-sample offset marginalized,
        low-z vs high-z residual step d_i per sample, pairwise tensions)
h0map : H0-tension <-> w(z) equivalence map — fix h on a grid, MAP-fit w0wa and
        lcdm on BAO+Dovekie(z>0.1)+CMB; gives d(w0,wa)/dH0 and the chi2(h) profile.
bfbin : nonparametric Savage-Dickey BF for 3-bin free w(z) vs LCDM from the
        session-4 binned chains (3D KDE at w=(-1,-1,-1)).

Usage: python experiments/session5.py step [0.05,0.1 | 0.2]
       python experiments/session5.py xres
       python experiments/session5.py h0map
       python experiments/session5.py bfbin
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
os.makedirs(RESULTS, exist_ok=True)

print("loading data ...", flush=True)
bao_all = cosmo.load_desi_bao(which="ALL")[0]
z_d, mu_d, invC_d = cosmo.load_dovekie()
z_p, mB_p, invC_p = cosmo.load_pantheonplus()
z_y, mu_y, invC_y = cosmo.load_desy5()
print(f"  BAO {len(bao_all['z'])} rows; SN dovekie={len(z_d)} pantheon+={len(z_p)} desy5={len(z_y)}",
      flush=True)

PRIORS = {
    "Om": (0.05, 0.9), "Ob": (0.005, 0.10), "h": (0.55, 0.85), "rd": (120.0, 165.0),
    "w": (-1.6, -0.4), "w0": (-1.8, -0.2), "wa": (-3.0, 3.0), "dm": (-0.3, 0.3),
}
X0 = {"Om": 0.31, "Ob": 0.049, "h": 0.68, "rd": 147.0, "w0": -0.9, "wa": -0.3, "dm": 0.0}


def sigma_from_dchi2(dchi2, k):
    from scipy.stats import chi2 as chi2dist, norm
    p = chi2dist.sf(abs(dchi2), k)
    z = float(norm.isf(p)) if p < 0.5 else 0.0
    return z * float(np.sign(dchi2)) if dchi2 != 0 else 0.0, float(p)


# ----------------------------------------------------------------------------
def make_logpost(model, zmask, with_dm, z_step, h_fixed=None, sn="dovekie"):
    if sn == "dovekie":
        z_s, mu_s, invC_s = z_d, mu_d, invC_d
    elif sn == "pantheonplus":
        z_s, mu_s, invC_s = z_p, mB_p, invC_p
    else:
        z_s, mu_s, invC_s = z_y, mu_y, invC_y
    if zmask is not None:
        m = z_s > zmask
        z_s, mu_s, invC_s = z_s[m], mu_s[m], invC_s[np.ix_(m, m)]
    step_mask = (z_s < z_step).astype(float) if with_dm else None
    names = ["Om", "Ob", "rd"] + ([] if h_fixed is not None else ["h"]) \
        + (["w0", "wa"] if model == "w0wa" else (["w"] if model == "wcdm" else [])) \
        + (["dm"] if with_dm else [])

    def logpost(x):
        d = dict(zip(names, x))
        for n in names:
            lo, hi = PRIORS[n]
            if not (lo <= d[n] <= hi):
                return -np.inf
        Om, Ob, rd = d["Om"], d["Ob"], d["rd"]
        h = h_fixed if h_fixed is not None else d["h"]
        w0 = d.get("w0", d.get("w", -1.0))
        wa = d.get("wa", 0.0)
        bg = Background(Om, w0, wa)
        H0 = 100.0 * h
        r = mu_s - bg.mu(z_s, H0)
        if with_dm:
            r = r - d["dm"] * step_mask
        chi2 = bao_chi2(bg, H0, rd, bao_all) + float(r @ invC_s @ r) \
            + cmb_compressed_chi2(Om, Ob, h, bg)
        return -0.5 * chi2

    return logpost, names


def run_mcmc(model, tag, zmask, with_dm, z_step, nsteps=14000, seed=41, sn="dovekie"):
    logpost, names = make_logpost(model, zmask, with_dm, z_step, sn=sn)
    x0 = np.array([X0[n] for n in names])
    step = np.array([0.01 if n in ("Om",) else 0.005 if n == "Ob" else
                     0.01 if n == "h" else 0.5 if n == "rd" else
                     0.05 if n in ("w0", "wa", "w") else 0.01 for n in names])
    t0 = time.time()
    samples, acc = metropolis(logpost, x0, nsteps=nsteps, burn=max(3500, nsteps // 4),
                              step=step, seed=seed)
    stride = max(1, len(samples) // 2500)
    sub = samples[::stride]
    lps = np.array([logpost(s) for s in sub])
    chi2_map = float(-2.0 * lps.max())
    summ = {"tag": tag, "chi2_map": chi2_map, "acceptance": float(acc),
            "mean": dict(zip(names, map(float, samples.mean(axis=0)))),
            "std": dict(zip(names, map(float, samples.std(axis=0)))),
            "seconds": round(time.time() - t0, 1)}
    np.save(os.path.join(RESULTS, f"chain_s5_{tag}.npy"), np.column_stack([sub, lps]))
    print(f"  [{tag}] chi2_map={chi2_map:.2f} acc={acc:.2f} "
          + " ".join(f"{n}={m:.3f}+-{s:.3f}" for n, m, s in
                     zip(names, samples.mean(axis=0), samples.std(axis=0)))
          + f" ({summ['seconds']}s)", flush=True)
    return summ


def step_mode():
    zsteps = [float(s) for s in sys.argv[2].split(",")] if len(sys.argv) > 2 else [0.05, 0.1]
    try:
        with open(os.path.join(RESULTS, "step_results.json")) as f:
            out = [e for e in json.load(f) if e["z_step"] not in zsteps]
    except Exception:
        out = []
    for zs in zsteps:
        res = {}
        for model in ["lcdm", "w0wa"]:
            for dm in [False, True]:
                tag = f"step{zs}_{model}_dm{int(dm)}"
                res[tag] = run_mcmc(model, tag, None, dm, zs)
        d_plain = res[f"step{zs}_w0wa_dm0"]["chi2_map"] - res[f"step{zs}_lcdm_dm0"]["chi2_map"]
        d_dm = res[f"step{zs}_w0wa_dm1"]["chi2_map"] - res[f"step{zs}_lcdm_dm1"]["chi2_map"]
        s_plain, _ = sigma_from_dchi2(d_plain, 2)
        s_dm, _ = sigma_from_dchi2(d_dm, 2)
        dm_post = res["step{0}_w0wa_dm1".format(zs)]
        entry = {"z_step": zs,
                 "dchi2_no_dm": float(d_plain), "sigma_no_dm": float(s_plain),
                 "dchi2_with_dm": float(d_dm), "sigma_with_dm": float(s_dm),
                 "dm_w0wa_mean": dm_post["mean"].get("dm"), "dm_w0wa_std": dm_post["std"].get("dm"),
                 "dm_lcdm_mean": res[f"step{zs}_lcdm_dm1"]["mean"].get("dm"),
                 "dm_lcdm_std": res[f"step{zs}_lcdm_dm1"]["std"].get("dm"),
                 "runs": res}
        print(f"  >> z_step={zs}: w0wa vs lcdm  no-dm dchi2={d_plain:+.2f} ({abs(s_plain):.2f}s)"
              f"  with-dm dchi2={d_dm:+.2f} ({abs(s_dm):.2f}s)  "
              f"dm(w0wa)={entry['dm_w0wa_mean']:+.4f}+-{entry['dm_w0wa_std']:.4f} mag", flush=True)
        out.append(entry)
        with open(os.path.join(RESULTS, "step_results.json"), "w") as f:
            json.dump(out, f, indent=2)
    print("STEP DONE", flush=True)


# ----------------------------------------------------------------------------
def xres_mode():
    """Cross-sample low-z residual certificate.
    Geometry-only LCDM from BAO+CMB; per-sample constant offset marginalized
    analytically; compare (low-z bin mean - high-z bin mean) of residuals."""
    from scipy.optimize import minimize
    # geometry-only LCDM fit
    def lp_geo(x):
        Om, Ob, h, rd = x
        if not (0.05 < Om < 0.9 and 0.005 < Ob < 0.10 and 0.55 < h < 0.85 and 120 < rd < 165):
            return -np.inf
        bg = Background(Om)
        return -0.5 * (bao_chi2(bg, 100 * h, rd, bao_all)
                       + cmb_compressed_chi2(Om, Ob, h, bg))
    r0 = minimize(lambda x: -lp_geo(x), [0.31, 0.049, 0.68, 147.0], method="Nelder-Mead",
                  options={"maxiter": 3000})
    Om, Ob, h, rd = r0.x
    bg = Background(Om)
    H0 = 100 * h
    print(f"  geometry-only LCDM (BAO+CMB): Om={Om:.4f} h={h:.4f} H0={H0:.2f} rd={rd:.2f} "
          f"chi2={-2*lp_geo(r0.x):.2f}", flush=True)

    samples = {"dovekie": (z_d, mu_d, invC_d), "pantheonplus": (z_p, mB_p, invC_p),
               "desy5": (z_y, mu_y, invC_y)}
    bins = [("lowz", 0.0, 0.1), ("mid", 0.1, 0.3), ("high", 0.3, 99.0)]
    out = {"geometry": {"Om": float(Om), "h": float(h), "H0": float(H0), "rd": float(rd),
                        "chi2": float(-2 * lp_geo(r0.x))}, "samples": {}}
    for name, (z, data, invC) in samples.items():
        r = data - bg.mu(z, H0)
        one = np.ones(len(z))
        C1 = invC @ one
        # analytic marginalization over per-sample constant offset M
        r_off = r - (r @ C1) / (one @ C1) * one
        entry = {}
        for bname, zlo, zhi in bins:
            m = (z >= zlo) & (z < zhi)
            idx = np.where(m)[0]
            Cbb = invC[np.ix_(idx, idx)]
            rbb = r_off[idx]
            ob = np.ones(len(idx))
            w = Cbb @ ob
            mean = float(rbb @ w / (ob @ w))
            var = float(1.0 / (ob @ w))
            entry[bname] = {"mean_resid": mean, "sigma": np.sqrt(var), "n": int(m.sum())}
        # low-z step relative to high-z baseline, per sample
        d = entry["lowz"]["mean_resid"] - entry["high"]["mean_resid"]
        sd = np.hypot(entry["lowz"]["sigma"], entry["high"]["sigma"])
        entry["lowz_step_d"] = float(d)
        entry["lowz_step_sigma"] = float(sd)
        entry["lowz_step_tension_vs_zero"] = float(abs(d) / sd)
        out["samples"][name] = entry
        print(f"  [{name}] low-z residual mean={entry['lowz']['mean_resid']:+.4f}"
              f"+-{entry['lowz']['sigma']:.4f} (n={entry['lowz']['n']}); "
              f"lowz-minus-highz step d={d:+.4f}+-{sd:.4f} ({abs(d)/sd:.2f}s vs 0)", flush=True)
    names = list(out["samples"])
    out["pairwise_lowz_step_tensions"] = {}
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = out["samples"][names[i]], out["samples"][names[j]]
            dd = a["lowz_step_d"] - b["lowz_step_d"]
            ss = np.hypot(a["lowz_step_sigma"], b["lowz_step_sigma"])
            out["pairwise_lowz_step_tensions"][f"{names[i]}_vs_{names[j]}"] = \
                {"delta": float(dd), "sigma": float(ss), "tension": float(abs(dd) / ss)}
            print(f"  >> low-z step {names[i]} vs {names[j]}: dd={dd:+.4f}+-{ss:.4f} "
                  f"({abs(dd)/ss:.2f}s)", flush=True)
    with open(os.path.join(RESULTS, "xres_results.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("XRES DONE", flush=True)


# ----------------------------------------------------------------------------
def h0map_mode():
    """Fix h on a grid; MAP-fit lcdm and w0wa (Nelder-Mead) on BAO+Dovekie(z>0.1)+CMB."""
    from scipy.optimize import minimize
    hgrid = np.arange(0.660, 0.741, 0.01)
    out = []
    for h in hgrid:
        row = {"h": float(h), "H0": float(100 * h)}
        for model in ["lcdm", "w0wa"]:
            logpost, names = make_logpost(model, 0.1, False, 0.1, h_fixed=float(h))
            best = None
            starts = [[0.30, 0.045, 147.0]] if model == "lcdm" \
                else [[0.31, 0.046, 143.5, -0.85, -0.3], [0.30, 0.045, 147.0, -1.0, 0.0]]
            for x0 in starts:
                x0 = np.array(x0, float)
                res = minimize(lambda x: -logpost(x), x0, method="Nelder-Mead",
                               options={"maxiter": 4000, "xatol": 1e-5, "fatol": 1e-7})
                if res.fun < (best[0] if best else np.inf):
                    best = (res.fun, res.x)
            chi2 = -2.0 * (-best[0])
            row[model] = {"chi2_map": float(chi2), "map": dict(zip(names, map(float, best[1])))}
            print(f"  h={h:.3f} {model}: chi2={chi2:.2f} "
                  + " ".join(f"{n}={v:.3f}" for n, v in zip(names, best[1])), flush=True)
        d = row["w0wa"]["chi2_map"] - row["lcdm"]["chi2_map"]
        row["dchi2_w0wa_lcdm"] = float(d)
        out.append(row)
    with open(os.path.join(RESULTS, "h0map_results.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("H0MAP DONE", flush=True)


# ----------------------------------------------------------------------------
def bfbin_mode():
    """Savage-Dickey BF for 3-bin w(z) vs LCDM from session-4 binned chains."""
    from scipy.stats import gaussian_kde
    prior_dens = 1.0 / ((PRIORS["w0"][1] - PRIORS["w0"][0])
                        * (3.0 - 1.0) ** 0  # placeholder, fixed below
                        )
    # w1, w2 priors in session4 were (-3, 1); w0 (-1.8, -0.2)
    prior_dens = 1.0 / (1.6 * 4.0 * 4.0)
    out = {}
    for tag in ["joint_all", "joint_fs8", "bao_cmb_fs8", "bao_sn_fs8", "fs8_sn_cmb"]:
        f = os.path.join(RESULTS, f"chain_s4_binned_{tag}.npy")
        if not os.path.exists(f):
            continue
        ch = np.load(f)
        wb = ch[:, 5:8]
        # subsample for KDE speed
        if len(wb) > 4000:
            wb = wb[:: max(1, len(wb) // 4000)]
        kde = gaussian_kde(wb.T)
        post_null = float(kde(np.array([[-1.0], [-1.0], [-1.0]]))[0])
        bf = post_null / prior_dens
        out[tag] = {"BF_binned_vs_lcdm_SD": float(bf), "post_density_at_null": post_null,
                    "nsamples": int(len(wb))}
        print(f"  [{tag}] BF(3-bin w vs LCDM) SD = {bf:.3g} "
              f"(post density {post_null:.3g} vs prior {prior_dens:.3g})", flush=True)
    with open(os.path.join(RESULTS, "bfbin_results.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("BFBIN DONE", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "xres"
    {"step": step_mode, "xres": xres_mode, "h0map": h0map_mode, "bfbin": bfbin_mode}[mode]()
