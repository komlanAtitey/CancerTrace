#!/usr/bin/env python3
# run_algorithm_2.py
# ------------------------------------------------------------
# CancerTrace — Algorithm 2 (Python port)
# ------------------------------------------------------------
import os
import sys
import argparse
import warnings
from pathlib import Path

# Use a non-interactive backend if no display is present
if not os.environ.get("DISPLAY"):
    os.environ["MPLBACKEND"] = "Agg"

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.linalg import toeplitz

# Optional: read RData
try:
    import pyreadr  # pip install pyreadr
except Exception:
    pyreadr = None


# ---------- helpers ----------
def near_psd(A: np.ndarray, eps: float = 1e-10) -> np.ndarray:
    """Project a symmetric matrix to the nearest PSD by flooring eigenvalues."""
    A = (A + A.T) / 2.0
    vals, vecs = np.linalg.eigh(A)
    vals[vals < eps] = eps
    return (vecs * vals) @ vecs.T


def ess(w: np.ndarray) -> float:
    w = np.asarray(w, float)
    s = w.sum()
    if not np.isfinite(s) or s <= 0:
        return np.nan
    p = w / s
    d = (p**2).sum()
    return np.nan if (not np.isfinite(d) or d <= 0) else 1.0 / d


def hdi_interval(x: np.ndarray, cred_mass: float = 0.95) -> tuple[float, float]:
    """HDI for a 1D sample/weights using the shortest interval."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return (0.0, 0.0)
    x = np.sort(x)
    m = int(np.floor(cred_mass * x.size))
    if m < 1:
        return (float(x.min()), float(x.max()))
    widths = x[m:] - x[:-m]
    j = int(np.argmin(widths))
    return (float(x[j]), float(x[j + m]))


def first_matching_column(df: pd.DataFrame, patterns: list[str], numeric_ok: bool = False) -> pd.Series:
    """Pick the first column whose name matches any pattern (regex). Fallback to first numeric column if allowed."""
    for pat in patterns:
        hits = df.filter(regex=pat, axis=1)
        if hits.shape[1] > 0:
            return hits.iloc[:, 0]
    if numeric_ok:
        numcols = df.select_dtypes(include=[np.number])
        if numcols.shape[1] > 0:
            return numcols.iloc[:, 0]
    raise ValueError(f"No matching columns found for patterns={patterns}.")


# ---------- core algorithm ----------
def cancertrace_algorithm_2(dataset) -> dict:
    """
    dataset: numeric vector (or 1-row/1-col DataFrame/ndarray)
    returns dict with:
      sir.est, obs.vec, hdi.range, driver.effect, density.new
    """
    # 0) sanitize
    if isinstance(dataset, (pd.DataFrame, np.ndarray, list, tuple)):
        arr = np.asarray(dataset, dtype=float)
        if arr.ndim == 2:
            if arr.shape[0] == 1 or arr.shape[1] == 1:
                dataset = arr.ravel()
            else:
                raise ValueError("`dataset` must be a numeric vector or 1-row/1-col table.")
        else:
            dataset = arr
    else:
        dataset = np.asarray(dataset, dtype=float)

    if dataset.size == 0 or not np.all(np.isfinite(dataset)):
        raise ValueError("`dataset` must be a finite numeric vector.")

    # 1) light augmentation (like R)
    arma_data = dataset * 10.0
    add_n = max(0, 17 - arma_data.size)
    if add_n > 0:
        rng_aug = np.random.default_rng(123)
        rnd = rng_aug.uniform(arma_data.min(), arma_data.max(), size=add_n)
        arma_data = np.concatenate([arma_data, rnd])

    # 2) hyperparameters
    coef = 1.15
    N = 1000
    n = arma_data.size
    H = 0.9
    nu2 = 1.78
    epsW = 1e-16
    minVar = 1e-12

    # 3) fractional-Brownian-like covariance
    j = np.arange(1, n + 1)
    cor_inc = 0.5 * (np.abs(j + 1) ** (2 * H) - 2 * np.abs(j) ** (2 * H) + np.abs(j - 1) ** (2 * H))
    Sigma = near_psd(nu2 * toeplitz(cor_inc))

    # 4) particle trajectories: n x N
    rng = np.random.default_rng(12345)
    noise_dist = rng.multivariate_normal(mean=np.zeros(n), cov=Sigma, size=N).T
    noise_dist[~np.isfinite(noise_dist)] = 0.0

    # 5) buffers
    obs_vec = arma_data.copy()
    particles = np.zeros((n, N))
    variance = np.zeros((n, N))
    inv_var = np.zeros((n, N))
    sir_est = np.zeros(n)
    w = np.ones(N) / N

    # 6) scale prior from observations
    gam_fun = np.abs(obs_vec) / max(coef, epsW)
    arma_coef_var = np.exp(gam_fun)
    arma_coef = float(np.mean(arma_coef_var))

    # 7) SIR loop
    for t in range(n):
        # propagate
        particles = arma_coef * noise_dist

        # obs variance at time t
        vt = np.exp(particles[t, :])
        vt[~np.isfinite(vt)] = minVar
        vt = np.maximum(vt, minVar)
        variance[t, :] = vt
        inv_var[t, :] = 1.0 / vt

        # log-likelihood under N(0, vt) for y_t
        ll = -0.5 * (np.log(2 * np.pi) + np.log(vt) + (obs_vec[t] ** 2) * inv_var[t, :])
        ll[~np.isfinite(ll)] = -1e6

        # weights
        lw = np.log(np.maximum(w, epsW)) + ll
        m = np.max(lw)
        w = np.exp(lw - m)
        s = w.sum()
        if (not np.isfinite(s)) or s <= 0:
            w = np.ones(N) / N
        else:
            w = w / s
        w[~np.isfinite(w)] = 1.0 / N
        w = w / w.sum()

        # SIR estimate
        sir_est[t] = float((w * particles[t, :]).sum())

        # ESS + resample
        cur_ess = ess(w)
        if not np.isfinite(cur_ess):
            cur_ess = N
        if cur_ess < 0.25 * N:
            idx = rng.choice(N, size=N, replace=True, p=np.maximum(w, epsW) / np.maximum(w.sum(), epsW))
            particles = particles[:, idx]
            noise_dist = noise_dist[:, idx]
            w = np.ones(N) / N

    # 8) post-processing
    sir_est[~np.isfinite(sir_est)] = 0.0
    sir_est = np.maximum(sir_est, 0.0)

    density_function = w.copy()
    density_function[~np.isfinite(density_function)] = 0.0
    density_function /= max(density_function.sum(), epsW)

    density_new = density_function[density_function > np.exp(-20)]
    if density_new.size == 0:
        density_new = density_function

    # 9) HDI
    lo, hi = hdi_interval(density_new, cred_mass=0.95)
    sel = (density_new >= lo) & (density_new <= hi)
    hdi_slice = density_new[sel]
    sd_slice = np.std(hdi_slice)
    if (not np.isfinite(sd_slice)) or sd_slice == 0:
        driver_effect = float(np.mean(hdi_slice))
    else:
        driver_effect = float(np.mean((hdi_slice / sd_slice) ** 2))

    return {
        "sir.est": sir_est,
        "obs.vec": obs_vec,
        "hdi.range": (float(lo), float(hi)),
        "driver.effect": driver_effect,
        "density.new": density_new,
    }


# ---------- data loading ----------
def load_inputs_from_rdata(p1="epithelial.level.time1.rdata",
                           p2="epithelial.level.time2.rdata",
                           p3="epithelial.level.time3.rdata"):
    if pyreadr is None:
        raise RuntimeError("pyreadr is not installed. `pip install pyreadr` or use CSV loader.")
    d1 = pyreadr.read_r(p1)
    d2 = pyreadr.read_r(p2)
    d3 = pyreadr.read_r(p3)
    df1 = next(iter(d1.values()))
    df2 = next(iter(d2.values()))
    df3 = next(iter(d3.values()))
    return df1, df2, df3


def load_inputs_from_csv(p1="epithelial.level.time1.csv",
                         p2="epithelial.level.time2.csv",
                         p3="epithelial.level.time3.csv"):
    return pd.read_csv(p1), pd.read_csv(p2), pd.read_csv(p3)


def build_input_frame(df1: pd.DataFrame, df2: pd.DataFrame, df3: pd.DataFrame) -> tuple[pd.DataFrame, pd.Series]:
    # Expect columns: level_1 / level_2 / level_3 and id.time1 (or fallbacks)
    level_1 = pd.to_numeric(first_matching_column(df1, ["^level_1$"], numeric_ok=True), errors="coerce")
    level_2 = pd.to_numeric(first_matching_column(df2, ["^level_2$"], numeric_ok=True), errors="coerce")
    level_3 = pd.to_numeric(first_matching_column(df3, ["^level_3$"], numeric_ok=True), errors="coerce")
    gene_id = first_matching_column(df1, [r"id\.time1$", r"id_time1$", r"^id$", r"^gene$"], numeric_ok=False).astype(str)

    # Align lengths safely
    m = min(len(level_1), len(level_2), len(level_3), len(gene_id))
    level_1, level_2, level_3, gene_id = level_1.iloc[:m], level_2.iloc[:m], level_3.iloc[:m], gene_id.iloc[:m]

    epithelial_gene_level = pd.DataFrame({"t1": level_1.values, "t2": level_2.values, "t3": level_3.values})
    return epithelial_gene_level, gene_id


# ---------- main ----------
def main():
    parser = argparse.ArgumentParser(description="Run CancerTrace Algorithm 2 (Python port).")
    parser.add_argument("--rdata", action="store_true", help="Load inputs from .rdata files (requires pyreadr).")
    parser.add_argument("--in1", default="epithelial.level.time1.rdata", help="Time1 file (.rdata or .csv).")
    parser.add_argument("--in2", default="epithelial.level.time2.rdata", help="Time2 file (.rdata or .csv).")
    parser.add_argument("--in3", default="epithelial.level.time3.rdata", help="Time3 file (.rdata or .csv).")
    parser.add_argument("--outpng", default="algorithm2_top20.png", help="Output PNG for the plot.")
    args = parser.parse_args()

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    # 1) Load data
    if args.rdata:
        df1, df2, df3 = load_inputs_from_rdata(args.in1, args.in2, args.in3)
    else:
        df1, df2, df3 = load_inputs_from_csv(
            args.in1.replace(".rdata", ".csv"),
            args.in2.replace(".rdata", ".csv"),
            args.in3.replace(".rdata", ".csv"),
        )

    # 2) Build per-gene frame
    epithelial_gene_level, gene_id = build_input_frame(df1, df2, df3)

    # 3) Run Algorithm 2 per gene (row-wise)
    effects = []
    for i in range(epithelial_gene_level.shape[0]):
        row = epithelial_gene_level.iloc[i, :].values
        res = cancertrace_algorithm_2(row)
        effects.append(res["driver.effect"])

    dr_coef = np.asarray(effects, float)
    dr_coef[np.isinf(dr_coef)] = 0.0

    gene_dr = pd.DataFrame({"gene.id": gene_id.values, "dr.coef": dr_coef})
    gene_dr = gene_dr.sort_values("dr.coef", ascending=False, kind="mergesort").reset_index(drop=True)
    print(gene_dr.head(10).to_string(index=False))

    # 4) Visualization (Top 20, known drivers flagged with “X”)
    top = gene_dr.head(20).copy()
    top["gene.id"] = pd.Categorical(top["gene.id"], categories=top["gene.id"], ordered=True)
    known_drivers = {"UHRF1", "CD82", "TRIM44", "APC", "CHEK1", "TP53INP1"}
    top["is_driver"] = top["gene.id"].isin(known_drivers)

    mx = float(np.nanmax(top["dr.coef"].values)) if len(top) else 1.0
    pad = mx * 0.08

    fig, ax = plt.subplots(figsize=(8, 8))
    y = np.arange(len(top))

    # left strip boxes
    ax.barh(y, [pad * 0] * len(top), height=0.9, left=-pad * 0.7, color="white", edgecolor="black")

    # X marks
    for yy, flag in enumerate(top["is_driver"].values):
        if flag:
            ax.scatter(-pad * 0.35, yy, marker="x", s=50, c="black", linewidths=1.5)

    # bars
    ax.barh(y, top["dr.coef"].values, color="darkred", height=0.8)

    # tiny right-end ticks
    eps = 1e-3
    for yy, v in enumerate(top["dr.coef"].values):
        ax.plot([v - eps, v + eps], [yy, yy], color="red", lw=1.2)

    ax.set_yticks(y)
    ax.set_yticklabels(top["gene.id"].astype(str).tolist())
    ax.set_xlim([-pad if pad > 0 else -0.1, mx * 1.05])
    ax.set_xlabel("dr.coef")
    ax.set_title("Top 20 genes (high → low) with known drivers marked by X")
    ax.grid(False)
    plt.tight_layout()

    outpath = Path(args.outpng)
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    print(f"[OK] Saved plot to: {outpath.resolve()}")

    # If a GUI is available, show the figure interactively
    if os.environ.get("DISPLAY"):
        plt.show()


if __name__ == "__main__":
    main()



