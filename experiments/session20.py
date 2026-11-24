"""session20.py — anchor-clean refit: lcdm / w0wa / quintom / ldpt on BAO +
Pantheon+ NON-CALIBRATOR SNe + CMB (+fs8). Does N2 (quintom dip) survive without the
Cepheid calibrators?  Compare chi2_map and w(z) vs session6 (Dovekie).

Usage: python experiments/session20.py fit [bao_sn_cmb | bao_sn_cmb_fs8]
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import Background, bao_chi2, sn_mB_chi2, cmb_compressed_chi2, fs8_chi2, metropolis
from models import BackgroundFn, lnF_ldpt, make_lnF_quintom

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

print("loading data ...", flush=True)
bao_all = cosmo.load_desi_bao(which="ALL")[0]
z_p, mB_p, invC_p = cosmo.load_pantheonplus()
raw = np.genfromtxt(os.path.join(cosmo.DATA, "sn_data", "PantheonPlus", "Pantheon+SH0ES.dat"),
                    names=True, dtype=None, encoding="utf-8")
IS_CAL = np.asarray(raw["IS_CALIBRATOR"], float)
keep = IS_CAL <= 0.5
idx = np.where(keep)[0]
z_sn = z_p[keep]
mB_sn = mB_p[keep]
C_full = np.linalg.inv(invC_p)
invC_sn = np.linalg.inv(C_full[np.ix_(idx, idx)])
z_f, fs8v, invC_f = cosmo.load_fs8()
print(f"  anchor-clean Pantheon+: {keep.sum()} SNe", flush=True)

MODEL_PARAMS = {
    "lcdm": ["Om", "Ob", "h", "rd"],
    "w0wa": ["Om", "Ob", "h", "rd", "w0", "wa"],
    "ldpt": ["Om", "Ob", "h", "rd", "f", "zc", "d"],
    "quintom": ["Om", "Ob", "h", "rd", "wb", "A", "zq", "dq"],
}
PRIOR = {"Om": (0.05, 0.9), "Ob": (0.005, 0.10), "h": (0.55, 0.85), "rd": (120.0, 165.0),
         "w0": (-1.8, -0.2), "wa": (-3.0, 3.0), "f": (0.0, 0.6), "zc": (0.0, 2.0),
         "d": (0.02, 1.0), "wb": (-1.6, -0.4), "A": (-1.5, 1.5), "zq": (0.0, 3.0),
         "dq": (0.05, 1.5), "s8": (0.2, 1.4)}
X0 = {"Om": 0.31, "Ob": 0.049, "h": 0.68, "rd": 147.0, "w0": -0.9, "wa": -0.4,
      "f": 0.2, "zc": 0.5, "d": 0.2, "wb": -0.95, "A": -0.3, "zq": 1.2, "dq": 0.4, "s8": 0.8}


def make_logpost(model, combo):
    use_fs8 = "fs8" in combo
    names = list(MODEL_PARAMS[model]) + (["s8"] if use_fs8 else [])

    def logpost(x):
        d = dict(zip(names, x))
        for n in names:
            lo, hi = PRIOR[n]
            if not (lo <= d[n] <= hi):
                return -np.inf
        Om, Ob, h, rd = d["Om"], d["Ob"], d["h"], d["rd"]
        if model == "lcdm":
            bg = Background(Om); lnF = None
        elif model == "w0wa":
            bg = Background(Om, d["w0"], d["wa"]); lnF = None
        elif model == "ldpt":
            lnF = lambda z: lnF_ldpt(z, d["f"], d["zc"], d["d"])
            bg = BackgroundFn(Om, lnF)
        else:
            lnF = make_lnF_quintom(d["wb"], d["A"], d["zq"], d["dq"])
            bg = BackgroundFn(Om, lnF)
        H0 = 100.0 * h
        chi2 = bao_chi2(bg, H0, rd, bao_all) + sn_mB_chi2(bg, H0, z_sn, mB_sn, invC_sn) \
            + cmb_compressed_chi2(Om, Ob, h, bg)
        if use_fs8:
            chi2 += fs8_chi2(Om, d["s8"], z_f, fs8v, invC_f, lnF=lnF)
        return -0.5 * chi2

    return logpost, names


def fit_mode():
    combo = sys.argv[2] if len(sys.argv) > 2 else "bao_sn_cmb"
    out = {}
    entry = {"combo": combo + "_pantheon_nocal", "models": {}}
    for model in ["lcdm", "w0wa", "ldpt", "quintom"]:
        logpost, names = make_logpost(model, combo)
        x0 = np.array([X0[n] for n in names])
        step = np.array([{"Om": 0.01, "Ob": 0.005, "h": 0.01, "rd": 0.5, "w0": 0.05,
                          "wa": 0.08, "f": 0.05, "zc": 0.05, "d": 0.03, "wb": 0.04,
                          "A": 0.08, "zq": 0.08, "dq": 0.05, "s8": 0.02}[n] for n in names])
        t0 = time.time()
        samples, acc = metropolis(logpost, x0, nsteps=14000, burn=4000, step=step, seed=171)
        stride = max(1, len(samples) // 2500)
        sub = samples[::stride]
        lps = np.array([logpost(s) for s in sub])
        chi2_map = float(-2.0 * lps.max())
        # w(z) at MAP
        dmap = dict(zip(names, sub[int(np.argmax(lps))]))
        zg = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 3.0])
        if model == "lcdm":
            w = np.full(len(zg), -1.0)
        elif model == "w0wa":
            w = dmap["w0"] + dmap["wa"] * zg / (1 + zg)
        elif model == "ldpt":
            from models import w_eff_ldpt
            w = w_eff_ldpt(zg, dmap["f"], dmap["zc"], dmap["d"])
        else:
            from models import w_quintom
            w = w_quintom(zg, dmap["wb"], dmap["A"], dmap["zq"], dmap["dq"])
        entry["models"][model] = {
            "chi2_map": chi2_map, "acceptance": float(acc),
            "mean": dict(zip(names, map(float, samples.mean(axis=0)))),
            "std": dict(zip(names, map(float, samples.std(axis=0)))),
            "w_z": w.tolist(), "seconds": round(time.time() - t0, 1)}
        print(f"  [{combo}/nocal] {model}: chi2_map={chi2_map:.2f} acc={acc:.2f} "
              f"w(z=0,0.5,1,1.5,2,3)={np.round(w,3).tolist()} ({entry['models'][model]['seconds']}s)",
              flush=True)
        np.save(os.path.join(RESULTS, f"chain_s20_{model}_{combo}_nocal.npy"),
                np.column_stack([sub, lps]))
    base = entry["models"]["lcdm"]["chi2_map"]
    for model, k in [("w0wa", 2), ("ldpt", 3), ("quintom", 4)]:
        d = entry["models"][model]["chi2_map"] - base
        entry[f"{model}_vs_lcdm_dchi2"] = float(d)
        print(f"  >> {model} vs lcdm: dchi2={d:+.2f} (k={k})", flush=True)
    with open(os.path.join(RESULTS, "session20_results.json"), "w") as f:
        json.dump(entry, f, indent=2)
    print("S20 DONE", flush=True)


if __name__ == "__main__":
    fit_mode()
