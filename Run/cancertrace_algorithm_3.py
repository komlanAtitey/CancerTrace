# cancertrace_algorithm_3.py
# ------------------------------------------------------------
# Python port of CancerTrace Algorithm 3 (R version)
# ------------------------------------------------------------
from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Sequence, Dict, Any, List, Tuple

# stats / ML
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from statsmodels.tsa.api import VAR
import warnings

# ------------------------------------------------------------
# Utilities: sanitation
# ------------------------------------------------------------
def _sanitize_numeric_vector(x, preferred_names=("level_1", "level_2", "level_3", "value")) -> np.ndarray:
    """
    Accepts a pandas Series/ DataFrame/ list-like.
    Returns a 1D numpy array of floats.
    """
    if isinstance(x, pd.DataFrame):
        have = [c for c in preferred_names if c in x.columns]
        if len(have) >= 1:
            x = x[have[0]].values
        else:
            num_cols = [c for c in x.columns if pd.api.types.is_numeric_dtype(x[c])]
            if not num_cols:
                raise ValueError("No numeric columns found in provided DataFrame.")
            x = x[num_cols[0]].values
    elif isinstance(x, pd.Series):
        x = x.values
    x = np.asarray(x, dtype=float).ravel()
    return x

def _sanitize_gene_vector(g) -> List[str]:
    """
    Accepts a list/Series/DataFrame; returns a list[str] of gene IDs.
    """
    if isinstance(g, pd.DataFrame):
        likely = [c for c in ("gene", "Gene", "symbol", "Symbol", "id", "ID") if c in g.columns]
        if not likely:
            raise ValueError("gene_vector is a DataFrame but no obvious gene column was found.")
        g = g[likely[0]]
    if isinstance(g, pd.Series):
        g = g.values
    return [str(s) for s in np.asarray(g).ravel()]

# ------------------------------------------------------------
# Time evolution with noise (row-stochastic Markov-like)
# ------------------------------------------------------------
def _evolve_markov_noise(v: np.ndarray, noise_sd: float = 0.05, seed: int | None = None) -> np.ndarray:
    v = np.asarray(v, dtype=float).ravel()
    n = v.shape[0]
    rng = np.random.default_rng(seed)
    P = np.eye(n) + rng.normal(loc=0.0, scale=noise_sd, size=(n, n))
    # row-stochastic
    row_sums = P.sum(axis=1, keepdims=True)
    # avoid division by zero
    row_sums[row_sums == 0] = 1.0
    P = P / row_sums
    evolved = v @ P
    in_min, in_max = float(np.nanmin(v)), float(np.nanmax(v))
    out_min, out_max = float(np.nanmin(evolved)), float(np.nanmax(evolved))
    if not np.isfinite(out_min) or not np.isfinite(out_max) or out_max == out_min:
        return np.full(n, (in_min + in_max) / 2.0)
    rescaled = (evolved - out_min) / (out_max - out_min) * (in_max - in_min) + in_min
    return rescaled

def _generate_evolved_matrix_with_causality(
    v1: np.ndarray,
    v2: np.ndarray,
    v3: np.ndarray,
    noise_sd: float = 0.1,
    seed: int = 42,
    inject_causality: bool = True,
    causal_weight: float = 0.8
) -> np.ndarray:
    # Early
    v1_1 = _evolve_markov_noise(v1, noise_sd, seed)
    v1_2 = _evolve_markov_noise(v1_1, noise_sd, seed)
    # Mid
    v2_1 = _evolve_markov_noise(v2, noise_sd, seed)
    v2_2 = _evolve_markov_noise(v2_1, noise_sd, seed)
    # Late
    v3_1 = _evolve_markov_noise(v3, noise_sd, seed)

    if inject_causality:
        rng = np.random.default_rng(seed)
        eps = rng.normal(0, 0.01, size=v3_1.shape[0])
        v3_2 = causal_weight * v1_2 + (1.0 - causal_weight) * v3_1 + eps
        # clamp to range of v3 for realism
        lo, hi = float(np.nanmin(v3)), float(np.nanmax(v3))
        v3_2 = np.clip(v3_2, lo, hi)
    else:
        v3_2 = _evolve_markov_noise(v3_1, noise_sd, seed)

    out = np.column_stack([v1, v1_1, v1_2, v2, v2_1, v2_2, v3, v3_1, v3_2])
    return out  # shape: (n_genes, 9)

# ------------------------------------------------------------
# Granger causality score (−log10 p)
# ------------------------------------------------------------
def _compute_granger_score(x: np.ndarray, y: np.ndarray, max_lag: int = 2) -> float:
    """
    Fit VAR on [x, y] and test if x -> y (Granger).
    Returns −log10(p) or 0 if fitting fails.
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    # build a T x 2 dataframe
    df = pd.DataFrame({"x": x, "y": y})
    # VAR needs enough observations
    if df.shape[0] < max(8, max_lag + 3):
        return 0.0
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = VAR(df)
            res = model.fit(maxlags=max_lag, trend="c")
            test = res.test_causality("y", ["x"], kind="f")  # does x cause y?
            pval = float(test.pvalue)
            return -np.log10(pval + 1e-10)
    except Exception:
        return 0.0

def _compute_CIS_matrix(
    expr_mat: pd.DataFrame, non_drivers: Sequence[str], drivers: Sequence[str], max_lag: int = 2
) -> pd.DataFrame:
    """
    expr_mat: DataFrame (genes x timepoints), index = gene names, columns = T1..T9
    Returns: DataFrame (non_drivers x drivers) of −log10 p
    """
    nd = [g for g in non_drivers if g in expr_mat.index]
    dr = [g for g in drivers if g in expr_mat.index]
    mat = np.zeros((len(nd), len(dr)), dtype=float)
    for i, ndg in enumerate(nd):
        x = expr_mat.loc[ndg].values
        for j, dg in enumerate(dr):
            y = expr_mat.loc[dg].values
            mat[i, j] = _compute_granger_score(x, y, max_lag=max_lag)
    return pd.DataFrame(mat, index=nd, columns=dr)

# ------------------------------------------------------------
# Top influencers per driver
# ------------------------------------------------------------
def _get_top_influencers_per_driver(CIS: pd.DataFrame, top_n: int = 5) -> Dict[str, pd.DataFrame]:
    out: Dict[str, pd.DataFrame] = {}
    for driver in CIS.columns:
        df = pd.DataFrame({
            "Driver_Gene": driver,
            "Non_Driver_Gene": CIS.index,
            "Influence_Score": CIS[driver].values
        }).sort_values("Influence_Score", ascending=False).head(top_n).reset_index(drop=True)
        out[driver] = df
    return out

# ------------------------------------------------------------
# Logistic “transformation likelihood” + CV AUC
# ------------------------------------------------------------
def _compute_transformation_likelihood(CIS: pd.DataFrame, driver_genes: Sequence[str], all_genes: Sequence[str]):
    influence = pd.Series(0.0, index=pd.Index(all_genes, name="Gene"))
    if CIS.shape[0] > 0:
        influence.loc[CIS.index] = CIS.sum(axis=1).values
    labels = pd.Series([1 if g in set(driver_genes) else 0 for g in all_genes], index=all_genes, name="Label")
    model_df = pd.DataFrame({"Gene": all_genes, "Label": labels.values, "Influence": influence.values})

    # Fit logistic regression (L2), then compute predicted prob
    if model_df["Influence"].var() == 0:
        # degenerate case
        model_df["Transformation_Likelihood"] = 0.5
        aucs, auc_mean = [0.5], 0.5
        model = None
        return {"model_df": model_df, "model": model, "aucs": aucs, "auc_mean": auc_mean}

    X = model_df[["Influence"]].values
    y = model_df["Label"].values
    lr = LogisticRegression(max_iter=1000, solver="lbfgs")
    lr.fit(X, y)
    model_df["Transformation_Likelihood"] = lr.predict_proba(X)[:, 1]

    # CV AUC
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    aucs: List[float] = []
    for tr, te in skf.split(X, y):
        y_te = y[te]
        # If test fold single class, skip
        if len(np.unique(y_te)) < 2:
            continue
        lr_fold = LogisticRegression(max_iter=1000, solver="lbfgs")
        lr_fold.fit(X[tr], y[tr])
        p = lr_fold.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y_te, p))
    auc_mean = float(np.mean(aucs)) if aucs else float("nan")

    return {"model_df": model_df, "model": lr, "aucs": aucs, "auc_mean": auc_mean}

# ------------------------------------------------------------
# Knockout GC analysis
# ------------------------------------------------------------
def _compute_gc_pair(non_driver_series: np.ndarray, driver_series: np.ndarray, max_lag: int = 2) -> Tuple[float, float]:
    """
    Returns (F_stat, p_value). If failure, returns (np.nan, np.nan).
    """
    df = pd.DataFrame({"x": np.asarray(non_driver_series, float), "y": np.asarray(driver_series, float)})
    if df.shape[0] < max(8, max_lag + 3):
        return (np.nan, np.nan)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = VAR(df).fit(maxlags=max_lag, trend="c")
            tst = res.test_causality("y", ["x"], kind="f")
            F = float(getattr(tst, "statistic", np.nan))
            p = float(getattr(tst, "pvalue", np.nan))
            return (F, p)
    except Exception:
        return (np.nan, np.nan)

def run_knockout_gc(
    expr_mat: pd.DataFrame,
    nd_vec: Sequence[str],
    d_vec: Sequence[str],
    max_lag: int = 2,
    knock_sd: float = 1e-4,
    seed: int = 42
) -> Dict[str, Any]:
    """
    expr_mat: genes x timepoints (index = genes)
    nd_vec:   non-driver gene names
    d_vec:    driver gene names
    Returns: dict with logp_orig, logp_knock, and a long results table
    """
    if expr_mat.index.isnull().any():
        raise ValueError("expr_mat must have valid gene rownames (index).")
    nd = [g for g in pd.unique(nd_vec) if g in expr_mat.index]
    dr = [g for g in pd.unique(d_vec) if g in expr_mat.index]
    if not nd or not dr:
        raise ValueError("No valid non-driver or driver genes found in expression matrix.")

    rng = np.random.default_rng(seed)
    rows = []
    for ndg in nd:
        for dg in dr:
            F0, p0 = _compute_gc_pair(expr_mat.loc[ndg].values, expr_mat.loc[dg].values, max_lag=max_lag)
            expr_knock = expr_mat.copy()
            expr_knock.loc[ndg] = rng.normal(0.0, knock_sd, size=expr_knock.shape[1])
            F1, p1 = _compute_gc_pair(expr_knock.loc[ndg].values, expr_knock.loc[dg].values, max_lag=max_lag)
            rows.append({
                "non_driver": ndg,
                "driver": dg,
                "F_orig": F0,
                "p_orig": p0,
                "F_knock": F1,
                "p_knock": p1,
                "logp_orig": -np.log10((p0 if np.isfinite(p0) else 1.0) + 1e-12),
                "logp_knock": -np.log10((p1 if np.isfinite(p1) else 1.0) + 1e-12),
            })
    df = pd.DataFrame(rows)
    return {"logp_orig": df["logp_orig"].values, "logp_knock": df["logp_knock"].values, "table": df}

# ------------------------------------------------------------
# Main entry point
# ------------------------------------------------------------
def cancertrace_algorithm_3(
    data_vector_1,   # Early  (array/Series/DataFrame with numeric column)
    data_vector_2,   # Mid
    data_vector_3,   # Late
    gene_vector,     # list/Series/DataFrame column of gene names (same length/order as above)
    driver_genes,    # list of known/putative drivers
    # ---- evolution / simulation controls ----
    noise_sd: float = 0.1,
    seed: int = 42,
    inject_causality: bool = True,
    causal_weight: float = 0.8,
    # ---- modeling controls ----
    max_lag: int = 2,
    top_n: int = 5,
    k_folds: int = 5,         # (kept for parity; we use internal 5-fold in likelihood)
    knock_sd: float = 1e-4
) -> Dict[str, Any]:
    # Sanitize inputs
    v1 = _sanitize_numeric_vector(data_vector_1)
    v2 = _sanitize_numeric_vector(data_vector_2)
    v3 = _sanitize_numeric_vector(data_vector_3)
    genes = _sanitize_gene_vector(gene_vector)

    n = len(genes)
    if not (len(v1) == n and len(v2) == n and len(v3) == n):
        raise ValueError(f"Length mismatch: gene_vector={n}, v1={len(v1)}, v2={len(v2)}, v3={len(v3)}")

    # Build time-augmented matrix (T1..T9)
    num_data = _generate_evolved_matrix_with_causality(
        v1, v2, v3, noise_sd=noise_sd, seed=seed,
        inject_causality=inject_causality, causal_weight=causal_weight
    )
    columns = [f"T{i}" for i in range(1, 10)]
    causality_data = pd.DataFrame(num_data, index=genes, columns=columns)

    # Filter driver/non-driver sets
    drivers = [g for g in driver_genes if g in causality_data.index]
    if not drivers:
        raise ValueError("None of the specified driver genes are present in the expression matrix rownames.")
    non_drivers = [g for g in causality_data.index if g not in drivers]
    if not non_drivers:
        raise ValueError("No non-driver genes remain after filtering.")

    # CIS + top influencers
    CIS = _compute_CIS_matrix(causality_data, non_drivers, drivers, max_lag=max_lag)
    top_influencers = _get_top_influencers_per_driver(CIS, top_n=top_n)

    # Likelihood + CV AUC
    like = _compute_transformation_likelihood(CIS, drivers, genes)
    likelihood_df = like["model_df"]
    auc_mean = like["auc_mean"]

    # Knockout analysis for ALL drivers
    per_driver = {}
    for drv in top_influencers.keys():
        nd_vec = top_influencers[drv]["Non_Driver_Gene"].tolist()
        ko = run_knockout_gc(causality_data, nd_vec, [drv], max_lag=max_lag, knock_sd=knock_sd, seed=seed)
        per_driver[drv] = ko["table"]

    knockout_table = pd.concat(
        [df.assign(Driver_Group=k) for k, df in per_driver.items()],
        axis=0, ignore_index=True
    )
    knockout_results = {
        "per_driver": per_driver,
        "table": knockout_table,
        "logp_orig": knockout_table["logp_orig"].values,
        "logp_knock": knockout_table["logp_knock"].values
    }

    return {
        "causality_data": causality_data,
        "CIS_matrix": CIS,
        "top_influencers": top_influencers,
        "likelihood_output": like,
        "likelihood_df": likelihood_df,
        "auc_mean": auc_mean,
        "knockout_results": knockout_results
    }

# ------------------------------------------------------------
# Example usage (uncomment and adapt to your environment)
# ------------------------------------------------------------
# if __name__ == "__main__":
#     # Example: data_vector_* could be pandas DataFrames with columns level_1/2/3
#     # and gene_vector a list/Series of gene names aligned to those rows.
#     result = cancertrace_algorithm_3(
#         data_vector_1 = df_time1["level_1"],
#         data_vector_2 = df_time2["level_2"],
#         data_vector_3 = df_time3["level_3"],
#         gene_vector   = df_time1["id.time1"],
#         driver_genes  = ["TP53INP1", "CCNL2", "VPS37D", "ATP11AUN", "FBXO6"]
#     )
#     print(result["auc_mean"])
#     print(result["CIS_matrix"].shape)
#     print(list(result["top_influencers"].keys()))

