"""session13.py — data tests for N6-N10 (HARNESS_N6_N10.md):

  angular : N7 bubble-wall test — dipole+quadrupole fit to low-z Pantheon+ residuals
            on the sphere (RA/DEC) vs monopole-only. A local expanding bubble/void
            gives l>=1 structure; a global calibration step does not.
  dust    : N8 environment test — low-z Pantheon+ residual slope vs MWEBV
            (dust column; complements the host-mass test §19).
  cert5   : N9/N10 — load chain_s10_cert5.npy: (a) cosine-oscillation template
            (N9 axion stirring) vs the 5 w bins; (b) sign alternation of bin
            residuals (N10 time-crystal) -> period bound.

Usage: python experiments/session13.py angular|dust|cert5|all
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "src"))
import cosmo
from cosmo import Background, bao_chi2, cmb_compressed_chi2

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results")

print("loading data ...", flush=True)
bao_all = cosmo.load_desi_bao(which="ALL")[0]
z_p, mB_p, invC_p = cosmo.load_pantheonplus()
raw = np.genfromtxt(os.path.join(cosmo.DATA, "sn_data", "PantheonPlus", "Pantheon+SH0ES.dat"),
                    names=True, dtype=None, encoding="utf-8")
RA = np.deg2rad(np.asarray(raw["RA"], float))
DEC = np.deg2rad(np.asarray(raw["DEC"], float))
MWEBV = np.asarray(raw["MWEBV"], float)
print(f"  Pantheon+ {len(z_p)} rows", flush=True)

from scipy.optimize import minimize


def geo_fit():
    def lp(x):
        Om, Ob, h, rd = x
        if not (0.05 < Om < 0.9 and 0.005 < Ob < 0.10 and 0.55 < h < 0.85 and 120 < rd < 165):
            return -np.inf
        bg = Background(Om)
        return -0.5 * (bao_chi2(bg, 100 * h, rd, bao_all)
                       + cmb_compressed_chi2(Om, Ob, h, bg))
    r = minimize(lambda x: -lp(x), [0.31, 0.049, 0.68, 147.0], method="Nelder-Mead",
                 options={"maxiter": 3000})
    return Background(r.x[0]), 100 * r.x[2], float(-2 * lp(r.x))


bg, H0, chi2_geo = geo_fit()
print(f"  geometry LCDM: H0={H0:.2f} chi2={chi2_geo:.2f}", flush=True)

out = {"geometry": {"H0": float(H0), "chi2": float(chi2_geo)}}


def resid(zmask):
    m = zmask(z_p)
    r = mB_p[m] - bg.mu(z_p[m], H0)
    idx = np.where(m)[0]
    C = invC_p[np.ix_(idx, idx)]
    one = np.ones(len(idx))
    w = C @ one
    r = r - (r @ w) / (one @ w)  # marginalize monopole
    return idx, r, C


def angular_mode():
    for zcut, tag in [(0.1, "z<0.1"), (0.3, "z<0.3")]:
        idx, r, C = resid(lambda z: z < zcut)
        n = len(idx)
        ra, dec = RA[idx], DEC[idx]
        # real spherical harmonics l=1 (3) and l=2 (5)
        x = np.sin(dec) * np.cos(ra)
        y = np.sin(dec) * np.sin(ra)
        zc = np.cos(dec)
        l1 = np.column_stack([x, y, zc])
        l2 = np.column_stack([x * y, y * zc, zc * x, zc ** 2 - 1.0 / 3.0, x ** 2 - y ** 2])
        res = {"n": int(n)}
        for name, X in [("monopole-only", np.ones((n, 1))),
                        ("dipole l=1", np.column_stack([np.ones((n, 1)), l1])),
                        ("l<=2", np.column_stack([np.ones((n, 1)), l1, l2]))]:
            XtWX = X.T @ C @ X
            beta = np.linalg.solve(XtWX, X.T @ C @ r)
            chi2 = float(r @ C @ r - beta @ (X.T @ C @ r))
            res[name] = {"chi2": chi2, "beta": beta.tolist()}
        d1 = res["monopole-only"]["chi2"] - res["dipole l=1"]["chi2"]
        d2 = res["dipole l=1"]["chi2"] - res["l<=2"]["chi2"]
        res["dchi2_dipole_vs_mono"] = float(d1)
        res["dchi2_quad_vs_dipole"] = float(d2)
        dip = np.array(res["dipole l=1"]["beta"][1:4])
        res["dipole_amplitude_mag"] = float(np.linalg.norm(dip))
        # sigma of amplitude: marginal cov
        X = np.column_stack([np.ones((n, 1)), l1])
        cov = np.linalg.inv(X.T @ C @ X)
        res["dipole_amplitude_sigma"] = float(np.sqrt(np.trace(cov[1:, 1:])))
        out.setdefault("angular", {})[tag] = res
        print(f"  [angular {tag}] n={n} dchi2(dipole-mono)={d1:+.2f} "
              f"dchi2(quad-dip)={d2:+.2f} |dipole|={res['dipole_amplitude_mag']:.4f}"
              f"+-{res['dipole_amplitude_sigma']:.4f} mag", flush=True)
    print("ANGULAR DONE", flush=True)


def dust_mode():
    for zcut, tag in [(0.1, "z<0.1"), (0.3, "z<0.3")]:
        idx, r, C = resid(lambda z: z < zcut)
        ebv = MWEBV[idx]
        X = np.column_stack([np.ones(len(idx)), ebv])
        XtWX = X.T @ C @ X
        beta = np.linalg.solve(XtWX, X.T @ C @ r)
        cov = np.linalg.inv(XtWX)
        slope, sig = float(beta[1]), float(np.sqrt(cov[1, 1]))
        out.setdefault("dust", {})[tag] = {"n": int(len(idx)), "slope_dr_dEBV": slope,
                                            "sigma": sig, "tension": abs(slope) / sig}
        print(f"  [dust {tag}] n={len(idx)} slope dr/dMWEBV={slope:+.4f}+-{sig:.4f} "
              f"({abs(slope)/sig:.2f}s)", flush=True)
    print("DUST DONE", flush=True)


def cert5_mode():
    ch = np.load(os.path.join(RESULTS, "chain_s10_cert5.npy"))
    wi = ch[:, :10][:, 5:10]  # 5 w bins
    m = wi.mean(axis=0)
    s = wi.std(axis=0)
    edges = [0.0, 0.3, 0.75, 1.5, 2.5, 99.0]
    zc = np.array([0.15, 0.525, 1.125, 2.0, 3.5])
    # N9: cosine template w = -1 + A cos(2 pi (z - z0)/P); coarse: fit A,P,z0 on bin means
    from scipy.optimize import minimize as mn
    def cmodel(p):
        A, P, z0 = p
        return -1.0 + A * np.cos(2 * np.pi * (zc - z0) / P)
    def chi2c(p):
        return float(np.sum((m - cmodel(p)) ** 2 / s ** 2))
    best = None
    for P0 in [1.0, 1.5, 2.0, 3.0, 5.0]:
        for z00 in [0.5, 1.0, 1.5]:
            r = mn(chi2c, [0.3, P0, z00], method="Nelder-Mead")
            if best is None or r.fun < best.fun:
                best = r
    chi2_cos = float(best.fun)
    chi2_flat = float(np.sum((m + 1.0) ** 2 / s ** 2))
    # N10: sign alternation of (w_bin - median)
    dev = m - np.median(m)
    signs = np.sign(dev)
    n_alt = int(np.sum(signs[1:] * signs[:-1] < 0))
    out["cert5tests"] = {
        "w_bins": m.tolist(), "w_std": s.tolist(), "z_centers": zc.tolist(),
        "N9_cosine": {"best_params": best.x.tolist(), "chi2": chi2_cos,
                      "chi2_flat_w=-1": chi2_flat, "dchi2": chi2_flat - chi2_cos},
        "N10_alternation": {"deviations": dev.tolist(), "n_sign_changes": n_alt,
                            "max_bins": 5},
    }
    print(f"  [cert5] w={np.round(m,3).tolist()} sd={np.round(s,3).tolist()}")
    print(f"  [N9 cosine] chi2={chi2_cos:.2f} vs flat-w=-1 chi2={chi2_flat:.2f} "
          f"(dchi2={chi2_flat-chi2_cos:+.2f}); params A,P,z0={np.round(best.x,3).tolist()}")
    print(f"  [N10 alternation] dev signs={signs.astype(int).tolist()} "
          f"n_changes={n_alt}/4", flush=True)
    print("CERT5TESTS DONE", flush=True)


if __name__ == "__main__":
    modes = sys.argv[1:] if len(sys.argv) > 1 else ["angular", "dust", "cert5"]
    if modes == ["all"]:
        modes = ["angular", "dust", "cert5"]
    for m in modes:
        {"angular": angular_mode, "dust": dust_mode, "cert5": cert5_mode}[m]()
    with open(os.path.join(RESULTS, "session13_results.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("S13 DONE", flush=True)
