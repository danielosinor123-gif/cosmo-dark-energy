"""session7.py — harnessing the H0 tension: which new-physics hypotheses absorb it?

For each model, chi2(h) profile on BAO + Dovekie(full SN) + CMB-compressed (R, wb),
with h FIXED on a grid and everything else MAP-fitted (Nelder-Mead, multi-start).
Model list:
  lcdm : baseline
  neff : + extra relativistic species dNeff (radiation density scaled)
  ede  : + early-dark-energy-like bump in E^2 near z~zc (relativistic scaling),
         A_x (1+z)^4 exp(-((z-zc)/d)^2)
Report: h at chi2 minimum, chi2 at H0=73.04 (SH0ES) minus chi2_min -> can the model
sit at SH0ES's H0 without a big chi2 penalty?  (the "harness" test)

Usage: python experiments/session7.py h0tension [model ...]
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import C_KMS, OGH2, bao_chi2, sn_hd_chi2, PLANCK18_COMPRESSED

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

print("loading data ...", flush=True)
bao_all = cosmo.load_desi_bao(which="ALL")[0]
z_sn, mu_sn, invC_sn = cosmo.load_dovekie()
ZG = np.concatenate([np.linspace(0.0, 5.0, 5001), np.geomspace(5.0, 5000.0, 4001)[1:]])
ZSTAR = PLANCK18_COMPRESSED["zstar"]


class BG7:
    """E^2 = Om(1+z)^3 + Or(1+z)^4 + Ode + extra(z); Or includes Neff."""

    def __init__(self, Om, h, Or, extra):
        self.Om, self.h, self.Or, self.extra = Om, h, Or, extra
        Ode = 1.0 - Om - Or * h * h  # Or here is Omega_r (physical); see Ez
        self.Ode = max(Ode, 1e-6)
        E = np.sqrt(Om * (1 + ZG) ** 3 + Or * h * h * (1 + ZG) ** 4 + self.Ode
                    + extra(ZG))
        self.zg = ZG
        self.Ig = np.concatenate([[0.0], np.cumsum(0.5 * (1 / E[1:] + 1 / E[:-1]) * np.diff(ZG))])

    def _E(self, z):
        z = np.asarray(z, float)
        return np.sqrt(self.Om * (1 + z) ** 3 + self.Or * self.h * self.h * (1 + z) ** 4
                       + self.Ode + self.extra(z))

    def DM(self, z, H0):
        return (C_KMS / H0) * np.interp(np.asarray(z, float), self.zg, self.Ig, left=0.0)

    def DH(self, z, H0):
        return C_KMS / (H0 * self._E(z))

    def DV(self, z, H0):
        z = np.asarray(z, float)
        return (z * self.DM(z, H0) ** 2 * self.DH(z, H0)) ** (1.0 / 3.0)

    def mu(self, z, H0):
        return 5.0 * np.log10((1 + np.asarray(z, float)) * self.DM(z, H0)) + 25.0


def make_extra(model, p):
    if model == "neff":
        dNeff = p
        return lambda z: 0.0 * z, dNeff  # handled in Or
    if model == "ede":
        Ax, zc, d = p
        return (lambda z: Ax * (1 + z) ** 4 * np.exp(-((z - zc) / d) ** 2)), None
    return lambda z: 0.0 * z, None


def chi2_at(h, model, x, use_sn=True):
    """x = [Om, Ob, rd] + model params. Returns total chi2 (BAO+SN+CMB-compressed)."""
    Om, Ob, rd = x[0], x[1], x[2]
    mp = x[3:]
    if not (0.05 < Om < 0.9 and 0.005 < Ob < 0.10 and 120 < rd < 165):
        return 1e12
    if model == "neff":
        dNeff = mp[0]
        if not (0.0 <= dNeff <= 3.0):
            return 1e12
    if model == "ede":
        Ax, zc, d = mp
        if not (0.0 <= Ax <= 0.5 and 300 < zc < 9000 and 50 < d < 4000):
            return 1e12
    Or = OGH2 * (1.0 + 0.2271 * (3.046 + (mp[0] if model == "neff" else 0.0)))
    extra, _ = make_extra(model, mp if model == "ede" else None)
    bg = BG7(Om, h, Or / h / h, extra)   # BG7 expects Omega_r = Or/h^2
    H0 = 100.0 * h
    chi2 = bao_chi2(bg, H0, rd, bao_all)
    if use_sn:
        chi2 += sn_hd_chi2(bg, H0, z_sn, mu_sn, invC_sn)
    DMstar = bg.DM(ZSTAR, H0)
    R = np.sqrt(Om) * H0 * DMstar / C_KMS
    wb = Ob * h * h
    chi2 += ((R - PLANCK18_COMPRESSED["R"][0]) / PLANCK18_COMPRESSED["R"][1]) ** 2
    chi2 += ((wb - PLANCK18_COMPRESSED["wb"][0]) / PLANCK18_COMPRESSED["wb"][1]) ** 2
    return float(chi2)


def h0tension_mode():
    models = sys.argv[2:] if len(sys.argv) > 2 else ["lcdm", "neff", "ede"]
    hgrid = np.arange(0.66, 0.761, 0.01)
    from scipy.optimize import minimize
    out = {}
    try:
        with open(os.path.join(RESULTS, "h0tension_results.json")) as f:
            out = {e["model"]: e for e in json.load(f)}
    except Exception:
        pass
    for model in models:
        if model in out and out[model].get("done"):
            print(f"  {model}: already done", flush=True)
            continue
        rows = []
        for h in hgrid:
            nmp = {"lcdm": 0, "neff": 1, "ede": 3}[model]
            starts = []
            base = [0.31, 0.049, 147.0]
            if model == "neff":
                starts = [base + [0.5], base + [1.5]]
            elif model == "ede":
                starts = [base + [0.05, 4000.0, 800.0], base + [0.1, 3000.0, 500.0],
                          base + [0.02, 6000.0, 1500.0]]
            else:
                starts = [base]
            best = None
            for x0 in starts:
                res = minimize(lambda x: chi2_at(h, model, x), np.array(x0, float),
                               method="Nelder-Mead",
                               options={"maxiter": 4000, "xatol": 1e-5, "fatol": 1e-6})
                if res.fun < (best[0] if best else np.inf):
                    best = (res.fun, res.x)
            rows.append({"h": float(h), "H0": float(100 * h), "chi2": float(best[0]),
                         "map": best[1].tolist()})
            print(f"  [{model}] h={h:.3f} chi2={best[0]:.2f} map="
                  + ",".join(f"{v:.4g}" for v in best[1]), flush=True)
        chi2s = np.array([r["chi2"] for r in rows])
        i_min = int(np.argmin(chi2s))
        chi2_7304 = float(np.interp(0.7304, [r["h"] for r in rows], chi2s))
        entry = {"model": model, "rows": rows, "h_min": rows[i_min]["H0"],
                 "chi2_min": float(chi2s.min()),
                 "chi2_at_H0_7304": chi2_7304,
                 "penalty_at_SH0ES": chi2_7304 - float(chi2s.min()), "done": True}
        out[model] = entry
        print(f"  >> {model}: H0_min={entry['h_min']:.2f} chi2_min={entry['chi2_min']:.2f} "
              f"chi2(H0=73.04)={chi2_7304:.2f} penalty={entry['penalty_at_SH0ES']:.2f}",
              flush=True)
        with open(os.path.join(RESULTS, "h0tension_results.json"), "w") as f:
            json.dump(list(out.values()), f, indent=2)
    print("H0TENSION DONE", flush=True)


if __name__ == "__main__":
    h0tension_mode()
