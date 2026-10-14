"""Simple descriptive CPM summaries, ordination and exploratory marker screens."""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests


def counts_per_million(counts: pd.DataFrame) -> pd.DataFrame:
    """Rows are samples, columns are taxa. CPM is relative abundance scaled by 1e6."""
    totals = counts.sum(axis=1)
    return counts.div(totals.replace(0, np.nan), axis=0) * 1_000_000


def grouped_means(cpm: pd.DataFrame, frame) -> pd.DataFrame:
    meta = frame.df.loc[cpm.index].copy()
    meta["_sample"] = cpm.index.astype(str)
    data = meta.join(cpm)
    keys = ["group"] + (["time"] if frame.has_time else [])
    values = data.groupby(keys, observed=True)[list(cpm.columns)].mean().reset_index()
    sizes = data.groupby(keys, observed=True).agg(n_samples=("_sample", "nunique")).reset_index()
    if meta["upper"].notna().any():
        units = data.groupby(keys, observed=True)["upper"].nunique().rename("n_upper").reset_index()
    else:
        units = sizes.rename(columns={"n_samples": "n_upper"})
    sizes = sizes.merge(units, on=keys, how="left")
    return sizes.merge(values, on=keys, how="left")


def pca(log_cpm: pd.DataFrame, n_axes: int = 3) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Centered PCA of log2(CPM+1), without variance scaling of taxa."""
    X = np.log2(log_cpm.to_numpy(float) + 1.0)
    X -= X.mean(axis=0, keepdims=True)
    X[:, np.std(X, axis=0) == 0] = 0
    U, s, vt = np.linalg.svd(X, full_matrices=False)
    k = min(n_axes, int((s > np.finfo(float).eps * max(X.shape) * (s[0] if len(s) else 1)).sum()))
    if k == 0:
        return pd.DataFrame(index=log_cpm.index), pd.DataFrame(columns=["feature", "PC", "loading"])
    scores = U[:, :k] * s[:k]
    variance = s[:k] ** 2 / max(float((s ** 2).sum()), np.finfo(float).tiny)
    coords = pd.DataFrame(scores, index=log_cpm.index, columns=[f"PC{i + 1}" for i in range(k)])
    coords.attrs["variance_explained"] = variance
    loadings = pd.DataFrame(vt[:k].T, index=log_cpm.columns, columns=[f"PC{i + 1}" for i in range(k)])
    long = loadings.rename_axis("feature").reset_index().melt(
        id_vars="feature", var_name="PC", value_name="loading")
    return coords, long


def lda(log_cpm: pd.DataFrame, labels: pd.Series, n_axes: int = 2) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Regularized Fisher linear discriminant projection; descriptive, not a hypothesis test."""
    valid = labels.reindex(log_cpm.index).notna()
    X = np.log2(log_cpm.loc[valid].to_numpy(float) + 1.0)
    y = labels.reindex(log_cpm.index)[valid].astype(str).to_numpy()
    classes = pd.unique(y)
    if len(classes) < 2 or len(y) <= len(classes):
        return pd.DataFrame(index=log_cpm.index[valid]), pd.DataFrame(columns=["feature", "axis", "loading"])
    X -= X.mean(0, keepdims=True)
    grand = X.mean(0)
    within_rows = np.zeros_like(X)
    between_rows = np.zeros_like(X)
    for c in classes:
        mask = y == c
        xc = X[mask]
        mc = xc.mean(0)
        within_rows[mask] = xc - mc
        between_rows[mask] = mc - grand
    # Dual-space LDA avoids constructing a taxa x taxa covariance matrix.
    within = within_rows @ within_rows.T
    between = between_rows @ between_rows.T
    scale = float(np.trace(within) / max(within.shape[0], 1))
    ridge = max(scale * 1e-3, np.finfo(float).eps)
    within.flat[:: within.shape[0] + 1] += ridge
    eig, vec = np.linalg.eig(np.linalg.pinv(within) @ between)
    order = np.argsort(eig.real)[::-1]
    keep = [i for i in order if eig[i].real > 1e-10][:min(n_axes, len(classes) - 1)]
    if not keep:
        return pd.DataFrame(index=log_cpm.index[valid]), pd.DataFrame(columns=["feature", "axis", "loading"])
    dual = vec[:, keep].real
    W = X.T @ dual
    scores = X @ W
    cols = [f"LD{i + 1}" for i in range(W.shape[1])]
    coords = pd.DataFrame(scores, index=log_cpm.index[valid], columns=cols)
    load = pd.DataFrame(W, index=log_cpm.columns, columns=cols)
    return coords, load.rename_axis("feature").reset_index().melt(
        id_vars="feature", var_name="axis", value_name="loading")


def lefse_style(cpm: pd.DataFrame, frame, min_prevalence: float = 0.1,
                alpha: float = 0.05) -> pd.DataFrame:
    """KW + pairwise Wilcoxon-consistency screen and standardized LDA-like separation.

    Tests use one value per upper experimental unit in each time stratum, avoiding
    treating repeated or nested subsamples as independent treatment replicates.
    This transparent screen is LEfSe-inspired; it is not the canonical LEfSe
    implementation and should be treated as a discovery screen.
    """
    meta = frame.df.loc[cpm.index].copy()
    if meta["upper"].notna().any():
        keys = ["upper"] + (["time"] if frame.has_time else [])
        vals = cpm.join(meta[["group", "upper"] + (["time"] if frame.has_time else [])])
        vals = vals.groupby(keys, observed=True).mean(numeric_only=True)
        group_map = meta.groupby(keys, observed=True)["group"].first()
        vals["group"] = group_map
        if frame.has_time:
            vals["time"] = vals.index.get_level_values("time")
    else:
        vals = meta[["group"] + (["time"] if frame.has_time else [])].join(cpm)
    strata = frame.time_levels if frame.has_time else ["all"]
    out = []
    for time in strata:
        d = vals if time == "all" else vals[vals["time"] == time]
        groups = list(pd.unique(d["group"]))
        if len(groups) < 2:
            continue
        p_values, effects = [], []
        tested = []
        for feature in cpm.columns:
            prevalent = float((d[feature] > 0).mean())
            if prevalent < min_prevalence:
                continue
            arrays = [d.loc[d["group"] == g, feature].to_numpy(float) for g in groups]
            if min(map(len, arrays)) < 2:
                continue
            try:
                p = float(kruskal(*arrays).pvalue)
            except ValueError:
                continue
            log_values = [np.log10(a + 1.0) for a in arrays]
            means = np.array([a.mean() for a in log_values])
            variances = np.array([a.var(ddof=1) for a in log_values])
            pooled_sd = np.sqrt(np.average(variances, weights=[len(a) for a in log_values]))
            separation = float((means.max() - means.min()) / pooled_sd) if pooled_sd > 0 else 0.0
            # A conservative consistency gate: every pair must pass a two-sided
            # rank-sum test, with Holm adjustment within the feature.
            pair_p = [mannwhitneyu(log_values[i], log_values[j], alternative="two-sided").pvalue
                      for i, j in itertools.combinations(range(len(groups)), 2)]
            pair_ok = bool(pair_p) and bool(np.all(multipletests(pair_p, method="holm")[1] < alpha))
            tested.append((feature, prevalent, p, separation, means, pair_ok))
            p_values.append(p)
        if not tested:
            continue
        q_values = multipletests(p_values, method="fdr_bh")[1]
        for item, q in zip(tested, q_values):
            feature, prevalent, p, separation, means, pair_ok = item
            winner = groups[int(np.argmax(means))]
            out.append(dict(
                time=time, taxon=feature, prevalence=prevalent, kw_p=p, kw_q=float(q),
                pairwise_consistent=pair_ok, lda_like_separation=separation,
                highest_mean_group=winner, mean_log10_cpm_by_group="; ".join(
                    f"{g}={v:.5g}" for g, v in zip(groups, means)),
                selected=bool(q < alpha and pair_ok),
                note="Discovery screen; CPM-based KW/Wilcoxon plus standardized separation, not canonical LEfSe LDA score"))
    return pd.DataFrame(out)
