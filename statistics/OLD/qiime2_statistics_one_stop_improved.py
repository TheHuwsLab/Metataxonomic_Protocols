#!/usr/bin/env python3
"""
QIIME2 level-7 taxonomy statistics for the RUSITEC and Piglet Cecum studies.

Design-aware analyses implemented here:

RUSITEC
--------
16 vessels, 4 treatments, 2 repeated timepoints (0 h and 4 h).
The vessel is the experimental subject/random effect. Alpha-diversity is
modelled as:

    outcome ~ Treatment * Timepoint + (1 | Vessel)

Piglet Cecum
------------
4 treatments and 2 phases (Pre-starter / Starter). The statistical unit for
Treatment is the pen when treatment was assigned at pen level. Therefore the
recommended model is:

    outcome ~ Treatment * Phase + (1 | Pen)

If the same animal was sampled in both phases, Animal/Piglet can instead be
used as the repeated-measures grouping factor (and should be explicitly
provided in the metadata). The script never invents this pairing.

Beta diversity
--------------
Bray-Curtis and Jaccard distance matrices are generated with PCoA plots.
Treatment PERMANOVA is performed separately within each timepoint/phase,
which preserves the repeated-measures structure of RUSITEC rather than
permuting the 0 h and 4 h observations freely. An exploratory full-factorial
PERMANOVA is also available for independent-sample designs.

Taxon statistics (optional)
----------------------------
Each taxon can be tested using a CLR-transformed taxon abundance with the
same study-specific mixed model. These results are labelled exploratory: a
compositional DA method such as ANCOM-BC2/ALDEx2 should be used for a formal
microbiome differential-abundance analysis.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import re
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.spatial.distance import pdist, squareform
from scipy.stats import kruskal, mannwhitneyu, ttest_ind, wilcoxon
from statsmodels.formula.api import mixedlm, ols
from statsmodels.stats.anova import anova_lm
from statsmodels.stats.multitest import multipletests

try:
    from patsy import build_design_matrices
except ImportError:  # pragma: no cover
    build_design_matrices = None


# -----------------------------------------------------------------------------
# Input handling
# -----------------------------------------------------------------------------

def read_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    df = df.dropna(how="all")
    if df.empty:
        raise SystemExit(f"Feature/taxon table is empty: {path}")

    # QIIME exports normally contain Feature ID / Taxon as the first column.
    # Keep the first column as the row identifier and coerce only abundance
    # columns to numeric.
    first_col = df.columns[0]
    df = df.set_index(first_col)
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    return df


def read_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path, sep=None, engine="python")
    meta = meta.dropna(how="all")
    if meta.empty:
        raise SystemExit(f"Metadata table is empty: {path}")
    sample_col = meta.columns[0]
    meta = meta.rename(columns={sample_col: "sample_id"})
    meta["sample_id"] = meta["sample_id"].astype(str).str.strip()
    if meta["sample_id"].duplicated().any():
        dup = meta.loc[meta["sample_id"].duplicated(), "sample_id"].tolist()
        raise SystemExit(f"Duplicate sample IDs in metadata: {dup[:10]}")
    return meta.set_index("sample_id")


def infer_orientation(table: pd.DataFrame, sample_ids: pd.Index) -> pd.DataFrame:
    """Return canonical features x samples table."""
    sample_ids = pd.Index(sample_ids.astype(str))
    row_hits = table.index.intersection(sample_ids)
    col_hits = table.columns.intersection(sample_ids)

    if len(col_hits) >= len(row_hits) and len(col_hits) > 0:
        out = table.loc[:, table.columns.intersection(sample_ids)].copy()
    elif len(row_hits) > 0:
        out = table.loc[table.index.intersection(sample_ids)].T.copy()
    else:
        raise SystemExit(
            "No overlapping samples were found between the abundance table "
            "and metadata. Check sample IDs and whether the table is a QIIME2 "
            "export."
        )

    # Convert abundance values only after orientation. Non-numeric annotation
    # columns should not normally be present after the first column was used as
    # the index, but coercion is safer than silently retaining strings.
    out = out.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    out = out.loc[:, ~out.columns.duplicated()]
    out = out.groupby(out.index).sum()
    return out


def normalise(table: pd.DataFrame, method: str, pseudocount: float = 1.0) -> pd.DataFrame:
    method = method.lower()
    table = table.astype(float)
    if method == "none":
        return table.copy()
    if method == "relative":
        denom = table.sum(axis=0).replace(0, np.nan)
        return table.div(denom, axis=1).fillna(0.0)
    if method == "clr":
        x = table + pseudocount
        gm = np.exp(np.log(x).mean(axis=0))
        return np.log(x.div(gm, axis=1))
    raise ValueError(f"Unknown normalisation method: {method}")


def clr_samples(table: pd.DataFrame, pseudocount: float = 1.0) -> pd.DataFrame:
    x = table + pseudocount
    gm = np.exp(np.log(x).mean(axis=0))
    return np.log(x.div(gm, axis=1))


# -----------------------------------------------------------------------------
# Diversity calculations
# -----------------------------------------------------------------------------

def shannon(col: pd.Series) -> float:
    p = col[col > 0].astype(float)
    if p.sum() <= 0:
        return np.nan
    p = p / p.sum()
    return float(-(p * np.log(p)).sum())


def simpson(col: pd.Series) -> float:
    p = col[col > 0].astype(float)
    if p.sum() <= 0:
        return np.nan
    p = p / p.sum()
    return float(1.0 - np.sum(p ** 2))


def alpha_metrics(table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample in table.columns:
        col = table[sample]
        rows.append(
            {
                "sample_id": sample,
                "observed_features": int((col > 0).sum()),
                "shannon": shannon(col),
                "simpson": simpson(col),
                "read_count_or_total": float(col.sum()),
            }
        )
    return pd.DataFrame.from_records(rows).set_index("sample_id")


def beta_distance(table: pd.DataFrame, metric: str) -> pd.DataFrame:
    metric = metric.lower()
    x = table.T.values
    labels = table.columns.astype(str)
    if metric == "braycurtis":
        dist = pdist(x, metric="braycurtis")
    elif metric == "jaccard":
        dist = pdist(x > 0, metric="jaccard")
    elif metric == "euclidean":
        dist = pdist(x, metric="euclidean")
    else:
        raise ValueError(f"Unknown beta metric: {metric}")
    return pd.DataFrame(squareform(dist), index=labels, columns=labels)


def pcoa(dist: pd.DataFrame, n_axes: int = 2) -> pd.DataFrame:
    """Classical PCoA / principal coordinates using double-centred distances."""
    d = dist.values.astype(float)
    n = d.shape[0]
    if n < 3:
        raise ValueError("PCoA requires at least 3 samples")
    d2 = d ** 2
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ d2 @ H
    eigvals, eigvecs = np.linalg.eigh(B)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    positive = eigvals > 1e-12
    eigvals = eigvals[positive]
    eigvecs = eigvecs[:, positive]
    k = min(n_axes, len(eigvals))
    coords = eigvecs[:, :k] * np.sqrt(eigvals[:k])
    cols = [f"PC{i+1}" for i in range(k)]
    return pd.DataFrame(coords, index=dist.index, columns=cols)


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def save_alpha_plots(alpha: pd.DataFrame, meta: pd.DataFrame, treatment_col: str, time_col: str, outdir: Path):
    merged = alpha.join(meta[[treatment_col, time_col]], how="inner")
    for metric in ["shannon", "simpson", "observed_features"]:
        plt.figure(figsize=(8, 6))
        # Treatment x time interaction displayed as separate boxes by timepoint.
        merged.boxplot(column=metric, by=[time_col, treatment_col], rot=45)
        plt.title(f"{metric.title()} by {treatment_col} and {time_col}")
        plt.suptitle("")
        plt.tight_layout()
        plt.savefig(outdir / f"alpha_{metric}_by_{treatment_col}_x_{time_col}.png", dpi=300)
        plt.close()

    # Time course: mean +/- SEM by treatment.
    for metric in ["shannon", "simpson", "observed_features"]:
        plot_df = merged.reset_index()
        plt.figure(figsize=(8, 6))
        for treatment, sub in plot_df.groupby(treatment_col):
            # Preserve original category order where possible.
            ordered = sub.groupby(time_col)[metric].agg(["mean", "sem"]).reset_index()
            plt.errorbar(
                np.arange(len(ordered)),
                ordered["mean"],
                yerr=ordered["sem"],
                marker="o",
                label=str(treatment),
            )
            plt.xticks(np.arange(len(ordered)), [str(x) for x in ordered[time_col]])
        plt.xlabel(time_col)
        plt.ylabel(metric.replace("_", " ").title())
        plt.legend(title=treatment_col)
        plt.tight_layout()
        plt.savefig(outdir / f"alpha_timecourse_{metric}.png", dpi=300)
        plt.close()


def save_pcoa_plot(coords: pd.DataFrame, meta: pd.DataFrame, treatment_col: str, time_col: str, outdir: Path):
    # Coordinates may already contain the grouping columns. Only join metadata
    # when they are not present to avoid pandas duplicate-column errors.
    if treatment_col in coords.columns and time_col in coords.columns:
        plot_df = coords.copy()
    else:
        plot_df = coords.join(meta[[treatment_col, time_col]], how="left")
    plt.figure(figsize=(8, 7))
    for (time, treatment), sub in plot_df.groupby([time_col, treatment_col], dropna=False):
        plt.scatter(sub["PC1"], sub["PC2"], s=60, label=f"{time} | {treatment}")
    plt.xlabel("PCoA1")
    plt.ylabel("PCoA2")
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    plt.savefig(outdir / f"pcoa_{treatment_col}_x_{time_col}.png", dpi=300, bbox_inches="tight")
    plt.close()


# -----------------------------------------------------------------------------
# Experimental-design handling
# -----------------------------------------------------------------------------

def choose_columns(args, meta: pd.DataFrame) -> dict:
    study = args.study.lower()
    treatment_col = args.treatment_col or "Treatment"

    if args.time_col:
        time_col = args.time_col
    elif study == "rusitec":
        time_col = "Timepoint"
    elif study == "piglet":
        time_col = "Phase"
    else:
        raise SystemExit("For --study custom you must supply --time-col.")

    if args.random_col:
        random_col = args.random_col
    elif study == "rusitec":
        random_col = "Vessel"
    elif study == "piglet" and "Pen" in meta.columns:
        random_col = "Pen"
    else:
        random_col = None

    required = [treatment_col, time_col]
    missing = [c for c in required if c not in meta.columns]
    if missing:
        raise SystemExit(f"Missing required metadata columns: {', '.join(missing)}")
    if random_col and random_col not in meta.columns:
        raise SystemExit(f"Random-effects column not found: {random_col}")

    return {
        "study": study,
        "treatment_col": treatment_col,
        "time_col": time_col,
        "random_col": random_col,
    }


def validate_design(meta: pd.DataFrame, design: dict) -> pd.DataFrame:
    """Check the uploaded metadata against the study design described in the package."""
    rows = []
    study = design["study"]
    tr = design["treatment_col"]
    tm = design["time_col"]
    rnd = design["random_col"]

    def add(check, status, detail):
        rows.append({"check": check, "status": status, "detail": detail})

    add("sample_count", "INFO", f"{len(meta)} samples available after metadata alignment")
    add("treatments", "INFO", f"{meta[tr].nunique(dropna=True)} levels: {', '.join(map(str, meta[tr].dropna().unique()))}")
    add("time_or_phase", "INFO", f"{meta[tm].nunique(dropna=True)} levels: {', '.join(map(str, meta[tm].dropna().unique()))}")

    if study == "rusitec":
        vessels = meta[rnd].nunique(dropna=True) if rnd else 0
        add("RUSITEC_vessels", "PASS" if vessels == 16 else "CHECK", f"{vessels} unique vessels; expected 16 from design document")
        reps = meta.groupby(rnd).size() if rnd else pd.Series(dtype=int)
        bad = reps[reps != 2]
        add("RUSITEC_repeated_timepoints", "PASS" if len(bad) == 0 and len(reps) == 16 else "CHECK",
            "Each vessel has exactly two observations" if len(bad) == 0 and len(reps) == 16 else f"Unexpected vessel replication: {bad.to_dict()}")
        tr_by_vessel = meta.groupby(rnd)[tr].nunique() if rnd else pd.Series(dtype=int)
        add("RUSITEC_treatment_constant_within_vessel", "PASS" if len(tr_by_vessel) and (tr_by_vessel == 1).all() else "CHECK",
            "Treatment is constant within vessel" if len(tr_by_vessel) and (tr_by_vessel == 1).all() else "At least one vessel has multiple treatment labels")
        n_vessels_by_treatment = meta.groupby(tr)[rnd].nunique().to_dict() if rnd else {}
        add("RUSITEC_treatment_replication", "PASS" if set(n_vessels_by_treatment.values()) == {4} else "CHECK",
            f"Vessels per treatment: {n_vessels_by_treatment}; expected 4 each")

    elif study == "piglet":
        add("Piglet_sample_count", "PASS" if len(meta) == 40 else "CHECK", f"{len(meta)} samples; design package describes 40 cecal samples")
        counts = meta.groupby([tr, tm]).size().to_dict()
        add("Piglet_treatment_x_phase_balance", "PASS" if len(set(counts.values())) == 1 and set(counts.values()) == {5} else "CHECK",
            f"Samples per treatment x phase: {counts}; expected 5 per combination in the supplied 40-sample design")
        if rnd:
            pens_by_treatment = meta.groupby(tr)[rnd].nunique().to_dict()
            add("Piglet_pen_replication", "PASS" if set(pens_by_treatment.values()) == {5} else "CHECK",
                f"Pens per treatment: {pens_by_treatment}; design describes 5 pens per treatment")
        else:
            add("Piglet_pen_random_effect", "CHECK",
                "No Pen column supplied. Treatment was assigned to 5 pens per treatment in the design, so pen should be included as a cluster/random effect whenever the sample-to-pen mapping is available.")

    return pd.DataFrame(rows)


def metadata_cleanup(meta: pd.DataFrame, design: dict) -> pd.DataFrame:
    meta = meta.copy()
    for col in [design["treatment_col"], design["time_col"], design["random_col"]]:
        if col:
            meta[col] = meta[col].astype(str).str.strip()
            meta.loc[meta[col].isin(["", "nan", "None"]), col] = np.nan
    return meta


# -----------------------------------------------------------------------------
# Alpha-diversity mixed models
# -----------------------------------------------------------------------------

def fit_mixed_model(data: pd.DataFrame, metric: str, treatment_col: str, time_col: str, random_col: str | None):
    d = data[[metric, treatment_col, time_col] + ([random_col] if random_col else [])].dropna().copy()
    if d.empty:
        return None, None, d

    d["__y"] = pd.to_numeric(d[metric], errors="coerce")
    d = d.dropna(subset=["__y"])
    d["__treatment"] = d[treatment_col].astype(str)
    d["__time"] = d[time_col].astype(str)
    if random_col:
        d["__random"] = d[random_col].astype(str)

    formula = "__y ~ C(__treatment) * C(__time)"
    if random_col:
        # Primary model: random-intercept linear mixed model.
        try:
            model = mixedlm(formula, d, groups=d["__random"], re_formula="1")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result = model.fit(reml=False, method="lbfgs", maxiter=2000, disp=False)
            if not result.converged:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    result = model.fit(reml=False, method="powell", maxiter=2000, disp=False)
            if result.converged and np.isfinite(result.params.to_numpy(dtype=float)).all():
                return result, "mixedlm", d
        except Exception:
            pass

        # Cluster-aware fallback: GEE keeps the repeated-measures grouping
        # without pretending the observations are independent.
        try:
            gee = smf.gee(
                formula,
                groups="__random",
                data=d,
                cov_struct=sm.cov_struct.Exchangeable(),
                family=sm.families.Gaussian(),
            )
            result = gee.fit()
            return result, "gee", d
        except Exception as exc:
            warnings.warn(
                f"MixedLM and GEE failed for {metric} ({exc}); falling back to OLS. "
                "The output will be labelled as an independence-assumption fallback."
            )

    # OLS is only used when no clustering variable was supplied or when both
    # cluster-aware models fail.
    model = ols(formula, data=d).fit()
    return model, "ols", d


def likelihood_ratio_pvalue(full_model, reduced_model):
    if full_model is None or reduced_model is None:
        return np.nan, np.nan
    lr = 2.0 * (full_model.llf - reduced_model.llf)
    df = max(int(full_model.df_modelwc - reduced_model.df_modelwc), 1)
    if lr < 0:
        return np.nan, np.nan
    from scipy.stats import chi2
    p = float(chi2.sf(lr, df))
    return float(lr), p


def _wald_term_test(result, term_name: str):
    design_info = result.model.data.design_info
    term_slices = design_info.term_name_slices
    if term_name not in term_slices:
        return np.nan, np.nan, np.nan

    sl = term_slices[term_name]
    params = result.fe_params if hasattr(result, "fe_params") else result.params
    cov = result.cov_params()
    # MixedLM cov_params contains fixed effects followed by variance terms.
    k = len(params)
    params = np.asarray(params, dtype=float).reshape(-1)
    cov = np.asarray(cov, dtype=float)[:k, :k]
    idx = np.arange(sl.start, sl.stop)
    b = params[idx]
    V = cov[np.ix_(idx, idx)]
    try:
        stat = float(b.T @ np.linalg.pinv(V) @ b)
        df = int(np.linalg.matrix_rank(V))
        from scipy.stats import chi2
        p = float(chi2.sf(stat, df))
        return stat, df, p
    except Exception:
        return np.nan, np.nan, np.nan


def alpha_omnibus(data: pd.DataFrame, metric: str, treatment_col: str, time_col: str, random_col: str | None):
    result, model_type, d = fit_mixed_model(data, metric, treatment_col, time_col, random_col)
    if result is None or d.empty:
        return pd.DataFrame()

    term_map = {
        "treatment": "C(__treatment)",
        "time": "C(__time)",
        "treatment:time": "C(__treatment):C(__time)",
    }
    rows = []
    for label, term_name in term_map.items():
        stat, df, p = _wald_term_test(result, term_name)
        rows.append({
            "metric": metric,
            "term": label,
            "model": {"mixedlm": "MixedLM random-intercept", "gee": "GEE exchangeable", "ols": "OLS"}.get(model_type, model_type),
            "test_statistic": stat,
            "df": df,
            "p_value": p,
        })

    out = pd.DataFrame(rows)
    if not out.empty:
        valid = out["p_value"].notna()
        if valid.any():
            out.loc[valid, "p_fdr"] = multipletests(out.loc[valid, "p_value"], method="fdr_bh")[1]
    return out


def _fixed_effect_contrast(result, rows_a: dict, rows_b: dict):
    """Wald contrast for fixed effects of a fitted statsmodels model."""
    if build_design_matrices is None:
        return np.nan, np.nan, np.nan
    design_info = result.model.data.design_info
    Xa = build_design_matrices([design_info], pd.DataFrame([rows_a]))[0]
    Xb = build_design_matrices([design_info], pd.DataFrame([rows_b]))[0]
    params = result.fe_params if hasattr(result, "fe_params") else result.params
    fe_names = list(params.index)
    Xa = np.asarray(Xa)[:, :len(fe_names)].ravel()
    Xb = np.asarray(Xb)[:, :len(fe_names)].ravel()
    c = Xa - Xb
    cov = np.asarray(result.cov_params())[:len(fe_names), :len(fe_names)]
    est = float(c @ params.values)
    var = float(c @ cov @ c)
    if var <= 0:
        return est, np.nan, np.nan
    se = math.sqrt(var)
    z = est / se
    from scipy.stats import norm
    p = 2 * norm.sf(abs(z))
    return est, se, float(p)


def alpha_pairwise(data: pd.DataFrame, metric: str, treatment_col: str, time_col: str, random_col: str | None):
    """Model-based marginal mean contrasts for MixedLM; robust fallback tests for OLS."""
    result, model_type, d = fit_mixed_model(data, metric, treatment_col, time_col, random_col)
    if result is None or d.empty:
        return pd.DataFrame()

    treatments = list(pd.unique(d["__treatment"]))
    times = list(pd.unique(d["__time"]))
    rows = []

    if model_type in {"mixedlm", "gee"}:
        for time in times:
            for a, b in itertools.combinations(treatments, 2):
                ra = {"__treatment": a, "__time": time}
                rb = {"__treatment": b, "__time": time}
                est, se, p = _fixed_effect_contrast(result, ra, rb)
                rows.append({
                    "metric": metric,
                    "contrast_type": "treatment_within_time",
                    "time": time,
                    "group_a": a,
                    "group_b": b,
                    "estimate_a_minus_b": est,
                    "SE": se,
                    "p_value": p,
                    "model": "MixedLM fixed-effect Wald contrast",
                })
        for treatment in treatments:
            for a, b in itertools.combinations(times, 2):
                ra = {"__treatment": treatment, "__time": a}
                rb = {"__treatment": treatment, "__time": b}
                est, se, p = _fixed_effect_contrast(result, ra, rb)
                rows.append({
                    "metric": metric,
                    "contrast_type": "time_within_treatment",
                    "time": f"{a} vs {b}",
                    "group_a": a,
                    "group_b": b,
                    "treatment": treatment,
                    "estimate_a_minus_b": est,
                    "SE": se,
                    "p_value": p,
                    "model": "MixedLM fixed-effect Wald contrast",
                })
    else:
        # For OLS fallback, use Welch tests for treatment-within-time and
        # independent tests for time-within-treatment. This is deliberately
        # labelled as a fallback because the intended study models are mixed.
        for time in times:
            sub = d[d["__time"] == time]
            for a, b in itertools.combinations(treatments, 2):
                xa = sub.loc[sub["__treatment"] == a, "__y"]
                xb = sub.loc[sub["__treatment"] == b, "__y"]
                if len(xa) and len(xb):
                    _, p = ttest_ind(xa, xb, equal_var=False, nan_policy="omit")
                    rows.append({
                        "metric": metric, "contrast_type": "treatment_within_time", "time": time,
                        "group_a": a, "group_b": b, "estimate_a_minus_b": xa.mean() - xb.mean(),
                        "SE": np.nan, "p_value": float(p), "model": "Welch t-test fallback",
                    })
        for treatment in treatments:
            sub = d[d["__treatment"] == treatment]
            for a, b in itertools.combinations(times, 2):
                xa = sub.loc[sub["__time"] == a, "__y"]
                xb = sub.loc[sub["__time"] == b, "__y"]
                if len(xa) and len(xb):
                    _, p = ttest_ind(xa, xb, equal_var=False, nan_policy="omit")
                    rows.append({
                        "metric": metric, "contrast_type": "time_within_treatment", "time": f"{a} vs {b}",
                        "group_a": a, "group_b": b, "treatment": treatment,
                        "estimate_a_minus_b": xa.mean() - xb.mean(), "SE": np.nan,
                        "p_value": float(p), "model": "Welch t-test fallback",
                    })

    out = pd.DataFrame(rows)
    if not out.empty:
        valid = out["p_value"].notna()
        if valid.any():
            out.loc[valid, "p_fdr"] = multipletests(out.loc[valid, "p_value"], method="fdr_bh")[1]
    return out


# -----------------------------------------------------------------------------
# PERMANOVA (single-factor; used within fixed timepoint/phase strata)
# -----------------------------------------------------------------------------

def permanova_pseudo_f(dist: pd.DataFrame, labels: pd.Series) -> tuple[float, float, float]:
    """One-factor PERMANOVA pseudo-F, no permutations performed here."""
    ids = list(dist.index)
    lab = pd.Series(labels, index=ids).astype(str)
    groups = list(pd.unique(lab))
    n = len(ids)
    g = len(groups)
    if g < 2 or n <= g:
        return np.nan, np.nan, np.nan

    # Gower-centred matrix.
    D2 = dist.loc[ids, ids].values ** 2
    H = np.eye(n) - np.ones((n, n)) / n
    G = -0.5 * H @ D2 @ H
    SST = np.trace(G)

    group_ss = 0.0
    for gr in groups:
        idx = np.flatnonzero(lab.values == gr)
        Hg = np.zeros((n, n))
        if len(idx):
            Hg[np.ix_(idx, idx)] = 1.0 / len(idx)
        group_ss += np.trace(Hg @ G)  # equivalent to centroid SS contribution

    ss_between = float(group_ss)
    ss_within = float(SST - ss_between)
    df_between = g - 1
    df_within = n - g
    ms_b = ss_between / df_between if df_between else np.nan
    ms_w = ss_within / df_within if df_within else np.nan
    F = ms_b / ms_w if ms_w > 0 else np.nan
    return float(F), ss_between, ss_within


def permuted_pvalue(observed: float, values: list[float]) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if not len(arr):
        return np.nan
    return float((1 + np.sum(arr >= observed)) / (1 + len(arr)))


def permanova_time_stratified(dist: pd.DataFrame, meta: pd.DataFrame, treatment_col: str, time_col: str, permutations: int, seed: int, outdir: Path):
    """Treatment PERMANOVA separately at each phase/timepoint.

    For RUSITEC this is design-aware: only one observation per vessel enters
    each test, so permutations are unrestricted among vessels *within the
    selected timepoint*. No 0 h / 4 h labels are permuted across one another.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for time in pd.unique(meta[time_col].dropna().astype(str)):
        ids = meta.index[meta[time_col].astype(str) == time].intersection(dist.index)
        subdist = dist.loc[ids, ids]
        labels = meta.loc[ids, treatment_col].astype(str)
        observed, ssb, ssw = permanova_pseudo_f(subdist, labels)
        if not np.isfinite(observed):
            continue
        null = []
        vals = labels.values.copy()
        for _ in range(permutations):
            perm = pd.Series(rng.permutation(vals), index=ids)
            f, _, _ = permanova_pseudo_f(subdist, perm)
            null.append(f)
        p = permuted_pvalue(observed, null)
        rows.append({
            "time_or_phase": time,
            "term": treatment_col,
            "n": len(ids),
            "groups": labels.nunique(),
            "pseudo_F": observed,
            "SS_between": ssb,
            "SS_within": ssw,
            "permutations": permutations,
            "p_value": p,
        })

    out = pd.DataFrame(rows)
    if not out.empty:
        valid = out["p_value"].notna()
        out.loc[valid, "p_fdr"] = multipletests(out.loc[valid, "p_value"], method="fdr_bh")[1]
    if not out.empty:
        out.to_csv(outdir / "beta_permanova_treatment_within_phase.tsv", sep="\t", index=False)
    return out


def paired_time_permanova(dist: pd.DataFrame, meta: pd.DataFrame, treatment_col: str, time_col: str, block_col: str | None, permutations: int, seed: int, outdir: Path, metric: str):
    """Permutation PERMANOVA for the time/phase effect within treatment.

    When each Vessel/Pen has exactly one sample at each of the two timepoints,
    labels are permuted only by swapping the two time labels within each block.
    This is the appropriate paired permutation structure for a two-timepoint
    repeated-measures design.
    """
    if not block_col or block_col not in meta.columns:
        return pd.DataFrame()
    rng = np.random.default_rng(seed)
    rows = []
    for treatment in pd.unique(meta[treatment_col].dropna().astype(str)):
        ids = meta.index[meta[treatment_col].astype(str) == treatment].intersection(dist.index)
        m = meta.loc[ids]
        good_blocks = []
        for block, sub in m.groupby(block_col):
            counts = sub[time_col].astype(str).value_counts()
            if len(counts) == 2 and counts.max() == 1:
                good_blocks.append(block)
        if len(good_blocks) < 2:
            continue
        keep = m[block_col].isin(good_blocks)
        ids = m.index[keep]
        subdist = dist.loc[ids, ids]
        labels = meta.loc[ids, time_col].astype(str)
        observed, ssb, ssw = permanova_pseudo_f(subdist, labels)
        if not np.isfinite(observed):
            continue
        null = []
        blocks = list(meta.loc[ids, block_col].astype(str).unique())
        for _ in range(permutations):
            perm = labels.copy()
            for block in blocks:
                bi = [i for i in ids if str(meta.loc[i, block_col]) == block]
                if len(bi) == 2 and rng.random() < 0.5:
                    a, b = perm.loc[bi[0]], perm.loc[bi[1]]
                    perm.loc[bi[0]], perm.loc[bi[1]] = b, a
            f, _, _ = permanova_pseudo_f(subdist, perm)
            null.append(f)
        p = permuted_pvalue(observed, null)
        rows.append({
            "treatment": treatment,
            "n_blocks": len(good_blocks),
            "term": time_col,
            "pseudo_F": observed,
            "SS_between": ssb,
            "SS_within": ssw,
            "permutations": permutations,
            "p_value": p,
            "note": f"Paired within-{block_col} permutation; labels swapped within block only.",
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        valid = out["p_value"].notna()
        out.loc[valid, "p_fdr"] = multipletests(out.loc[valid, "p_value"], method="fdr_bh")[1]
    out.to_csv(outdir / f"beta_permanova_{time_col}_within_treatment_{metric}.tsv", sep="\t", index=False)
    return out


def paired_change_magnitude(dist: pd.DataFrame, meta: pd.DataFrame, treatment_col: str, time_col: str, block_col: str | None, permutations: int, seed: int, outdir: Path, metric: str):
    """Test whether treatment changes the magnitude of the paired phase/time shift.

    This is a supplementary interaction-oriented analysis: it compares the
    within-block distance between the two phases across treatments. It tests
    magnitude of change, not its direction in community space.
    """
    if not block_col or block_col not in meta.columns:
        return pd.DataFrame()
    pairs = []
    for block, sub in meta.groupby(block_col):
        sub = sub[sub.index.isin(dist.index)]
        if sub[time_col].nunique() != 2 or len(sub) != 2 or sub[treatment_col].nunique() != 1:
            continue
        ids = list(sub.index)
        d = float(dist.loc[ids[0], ids[1]])
        pairs.append({"block": str(block), "treatment": str(sub[treatment_col].iloc[0]), "change_distance": d})
    pair_df = pd.DataFrame(pairs)
    if pair_df.empty or pair_df["treatment"].nunique() < 2:
        return pair_df

    observed_stat = float(kruskal(*[g["change_distance"].to_numpy() for _, g in pair_df.groupby("treatment")]).statistic)
    rng = np.random.default_rng(seed)
    labels = pair_df["treatment"].to_numpy().copy()
    values = pair_df["change_distance"].to_numpy()
    null = []
    for _ in range(permutations):
        perm = rng.permutation(labels)
        groups = [values[perm == lab] for lab in np.unique(labels)]
        h = float(kruskal(*groups).statistic)
        null.append(h)
    p = permuted_pvalue(observed_stat, null)
    global_out = pd.DataFrame([{
        "metric": metric,
        "test": "Treatment effect on paired phase/time change magnitude",
        "n_blocks": len(pair_df),
        "Kruskal_H": observed_stat,
        "permutations": permutations,
        "p_value": p,
        "note": f"Within-{block_col} {time_col}-pair distance; this is a magnitude-of-change test, not a full directional multivariate interaction test.",
    }])
    # Pairwise treatment contrasts on change magnitude.
    pair_rows = []
    treatments = list(pd.unique(pair_df["treatment"]))
    for a, b in itertools.combinations(treatments, 2):
        xa = pair_df.loc[pair_df["treatment"] == a, "change_distance"]
        xb = pair_df.loc[pair_df["treatment"] == b, "change_distance"]
        if len(xa) and len(xb):
            u, pp = mannwhitneyu(xa, xb, alternative="two-sided")
            pair_rows.append({"metric": metric, "group_a": a, "group_b": b, "U": u, "p_value": pp})
    pair_out = pd.DataFrame(pair_rows)
    if not pair_out.empty:
        pair_out["p_fdr"] = multipletests(pair_out["p_value"], method="fdr_bh")[1]
    global_out.to_csv(outdir / f"beta_change_magnitude_global_{metric}.tsv", sep="\t", index=False)
    pair_out.to_csv(outdir / f"beta_change_magnitude_pairwise_{metric}.tsv", sep="\t", index=False)
    pair_df.to_csv(outdir / f"beta_paired_change_distances_{metric}.tsv", sep="\t", index=False)
    return global_out


def exploratory_full_factorial_permanova(dist: pd.DataFrame, meta: pd.DataFrame, treatment_col: str, time_col: str, permutations: int, seed: int, outdir: Path):
    """Exploratory full-factorial PERMANOVA via composite groups.

    This does not model repeated vessels/pens. It is retained as a visual/
    exploratory result and is not the primary RUSITEC inferential analysis.
    """
    ids = dist.index.intersection(meta.index)
    labels = meta.loc[ids, [treatment_col, time_col]].astype(str).agg(" | ".join, axis=1)
    observed, ssb, ssw = permanova_pseudo_f(dist.loc[ids, ids], labels)
    if not np.isfinite(observed):
        return pd.DataFrame()
    rng = np.random.default_rng(seed)
    vals = labels.values.copy()
    null = []
    for _ in range(permutations):
        perm = pd.Series(rng.permutation(vals), index=ids)
        f, _, _ = permanova_pseudo_f(dist.loc[ids, ids], perm)
        null.append(f)
    p = permuted_pvalue(observed, null)
    out = pd.DataFrame([{
        "analysis": "exploratory composite Treatment | Time/Phase PERMANOVA",
        "n": len(ids),
        "groups": labels.nunique(),
        "pseudo_F": observed,
        "SS_between": ssb,
        "SS_within": ssw,
        "permutations": permutations,
        "p_value": p,
        "note": "Does not account for repeated Vessel/Pen/Animal structure; use as exploratory only.",
    }])
    out.to_csv(outdir / "beta_permanova_exploratory_composite.tsv", sep="\t", index=False)
    return out


# -----------------------------------------------------------------------------
# Optional taxon-level CLR mixed models
# -----------------------------------------------------------------------------

def taxon_statistics(table: pd.DataFrame, meta: pd.DataFrame, design: dict, outdir: Path, min_prevalence: float, min_mean_relative: float):
    """Exploratory per-taxon CLR mixed/cluster-aware models with BH correction."""
    clr = clr_samples(table)
    rel = table.div(table.sum(axis=0).replace(0, np.nan), axis=1).fillna(0)
    prevalence = (table > 0).mean(axis=1)
    mean_rel = rel.mean(axis=1)
    keep = (prevalence >= min_prevalence) & (mean_rel >= min_mean_relative)
    taxa = list(table.index[keep])

    rows = []
    for taxon in taxa:
        d = pd.DataFrame({
            "__y": clr.loc[taxon].values,
            "__sample": clr.columns,
        }).set_index("__sample").join(meta, how="inner")
        d = d.dropna(subset=[design["treatment_col"], design["time_col"]])
        if len(d) < 8:
            continue
        d["__treatment"] = d[design["treatment_col"]].astype(str)
        d["__time"] = d[design["time_col"]].astype(str)

        model_type = None
        result = None
        try:
            if design["random_col"]:
                d["__random"] = d[design["random_col"]].astype(str)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    mm = mixedlm("__y ~ C(__treatment) * C(__time)", d, groups=d["__random"], re_formula="1")
                    result = mm.fit(reml=False, method="lbfgs", maxiter=500, disp=False)
                    if not result.converged:
                        result = mm.fit(reml=False, method="powell", maxiter=500, disp=False)
                if result.converged and np.isfinite(result.params.to_numpy(dtype=float)).all():
                    model_type = "MixedLM"
            if result is None or model_type is None:
                if design["random_col"]:
                    gee = smf.gee(
                        "__y ~ C(__treatment) * C(__time)",
                        groups="__random",
                        data=d,
                        cov_struct=sm.cov_struct.Exchangeable(),
                        family=sm.families.Gaussian(),
                    )
                    result = gee.fit()
                    model_type = "GEE"
                else:
                    result = ols("__y ~ C(__treatment) * C(__time)", d).fit()
                    model_type = "OLS"

            tests = {}
            for label, term_name in {
                "treatment": "C(__treatment)",
                "time": "C(__time)",
                "treatment_x_time": "C(__treatment):C(__time)",
            }.items():
                stat, dfree, pval = _wald_term_test(result, term_name)
                tests[label] = pval

            rows.append({
                "taxon": taxon,
                "prevalence": float(prevalence[taxon]),
                "mean_relative_abundance": float(mean_rel[taxon]),
                "model": model_type,
                "p_treatment": tests["treatment"],
                "p_time_or_phase": tests["time"],
                "p_treatment_x_time": tests["treatment_x_time"],
            })
        except Exception:
            continue

    columns = [
        "taxon", "prevalence", "mean_relative_abundance", "model",
        "p_treatment", "p_time_or_phase", "p_treatment_x_time",
        "q_treatment", "q_time_or_phase", "q_treatment_x_time",
    ]
    out = pd.DataFrame(rows, columns=columns)
    if not out.empty:
        for col in ["p_treatment", "p_time_or_phase", "p_treatment_x_time"]:
            qcol = col.replace("p_", "q_")
            valid = out[col].notna()
            if valid.any():
                out.loc[valid, qcol] = multipletests(out.loc[valid, col], method="fdr_bh")[1]
        out = out.sort_values(["q_treatment_x_time", "q_treatment", "q_time_or_phase"], na_position="last")
    out.to_csv(outdir / "taxon_clr_mixed_model_results.tsv", sep="\t", index=False)
    return out


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(
        description="Design-aware statistics for QIIME2 level-7 taxonomy tables from the RUSITEC and Piglet Cecum studies."
    )
    p.add_argument("--table", required=True, help="Level-7 taxonomy abundance CSV/TSV exported from QIIME2.")
    p.add_argument("--metadata", required=True, help="Sample metadata CSV/TSV; first column must be sample ID.")
    p.add_argument("--study", required=True, choices=["rusitec", "piglet", "custom"], help="Study design.")
    p.add_argument("--treatment-col", default="Treatment", help="Treatment metadata column (default: Treatment).")
    p.add_argument("--time-col", default=None, help="Time/phase metadata column. RUSITEC defaults to Timepoint; Piglet defaults to Phase.")
    p.add_argument("--random-col", default=None, help="Random/repeated-measures grouping column. RUSITEC defaults to Vessel; Piglet auto-detects Pen if present.")
    p.add_argument("--normalization", default="relative", choices=["none", "relative", "clr"], help="Normalisation used for diversity calculations; alpha is calculated from relative proportions.")
    p.add_argument("--beta-metrics", default="braycurtis,jaccard", help="Comma-separated beta metrics.")
    p.add_argument("--permutations", type=int, default=999, help="Number of PERMANOVA permutations (default: 999).")
    p.add_argument("--seed", type=int, default=42, help="Random seed for permutations and reproducibility.")
    p.add_argument("--taxa-stats", action="store_true", help="Run exploratory per-taxon CLR mixed-model tests.")
    p.add_argument("--min-prevalence", type=float, default=0.10, help="Taxon prevalence threshold for --taxa-stats (default: 0.10).")
    p.add_argument("--min-mean-relative", type=float, default=1e-5, help="Minimum mean relative abundance for --taxa-stats.")
    p.add_argument("--outdir", default="statistics_output", help="Output directory.")
    p.add_argument("--no-exploratory-factorial-permanova", action="store_true", help="Skip the composite-factor exploratory PERMANOVA.")
    args = p.parse_args()

    if args.permutations < 99:
        raise SystemExit("Use at least 99 permutations; 999 or 9999 is recommended for final results.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    table = read_table(Path(args.table))
    meta = read_metadata(Path(args.metadata))
    table = infer_orientation(table, meta.index)
    meta.index = meta.index.astype(str).str.strip()

    # Align exactly to abundance-table sample order.
    samples = table.columns.intersection(meta.index)
    if len(samples) == 0:
        raise SystemExit("No samples remain after aligning abundance table and metadata.")
    table = table.loc[:, samples]
    meta = meta.loc[samples].copy()

    design = choose_columns(args, meta)
    meta = metadata_cleanup(meta, design)
    validation = validate_design(meta, design)
    validation.to_csv(outdir / "design_validation.tsv", sep="\t", index=False)
    for _, row in validation.iterrows():
        if row["status"] == "CHECK":
            print(f"DESIGN CHECK: {row['check']}: {row['detail']}", file=sys.stderr)

    # Remove samples with missing design information.
    required = [design["treatment_col"], design["time_col"]] + ([design["random_col"]] if design["random_col"] else [])
    complete = meta[required].notna().all(axis=1)
    if (~complete).any():
        removed = meta.index[~complete].tolist()
        pd.DataFrame({"sample_id": removed, "reason": "missing design metadata"}).to_csv(
            outdir / "excluded_samples.tsv", sep="\t", index=False
        )
    meta = meta.loc[complete]
    table = table.loc[:, meta.index]

    if table.empty:
        raise SystemExit("No samples remain after design metadata filtering.")

    # Ensure samples have non-zero total abundance.
    nonzero = table.sum(axis=0) > 0
    if (~nonzero).any():
        removed = table.columns[~nonzero].tolist()
        pd.DataFrame({"sample_id": removed, "reason": "zero total abundance"}).to_csv(
            outdir / "excluded_zero_abundance_samples.tsv", sep="\t", index=False
        )
        table = table.loc[:, nonzero]
        meta = meta.loc[table.columns]

    # Record the exact design and command-line settings.
    params = vars(args).copy()
    params["resolved_design"] = design
    params["n_samples"] = int(table.shape[1])
    params["n_features_or_taxa"] = int(table.shape[0])
    with open(outdir / "analysis_parameters.json", "w") as fh:
        json.dump(params, fh, indent=2, default=str)

    design_summary = pd.DataFrame([
        {
            "study": design["study"],
            "n_samples": table.shape[1],
            "n_features_or_taxa": table.shape[0],
            "treatments": "; ".join(sorted(meta[design["treatment_col"]].dropna().astype(str).unique())),
            "time_or_phase": "; ".join(sorted(meta[design["time_col"]].dropna().astype(str).unique())),
            "random_effect": design["random_col"] or "NONE (OLS/fixed-effect fallback)",
        }
    ])
    design_summary.to_csv(outdir / "design_summary.tsv", sep="\t", index=False)
    meta.to_csv(outdir / "analysis_metadata.tsv", sep="\t")

    # Relative abundance is appropriate for Shannon/Simpson and Bray-Curtis.
    rel = normalise(table, "relative")
    rel.to_csv(outdir / "relative_abundance.tsv", sep="\t")

    alpha = alpha_metrics(table)
    alpha = alpha.join(meta, how="inner")
    alpha.to_csv(outdir / "alpha_diversity.tsv", sep="\t")
    save_alpha_plots(alpha[["shannon", "simpson", "observed_features", "read_count_or_total"]], meta, design["treatment_col"], design["time_col"], outdir)

    # Alpha inferential analysis.
    omnibus = pd.concat(
        [alpha_omnibus(alpha, metric, design["treatment_col"], design["time_col"], design["random_col"])
         for metric in ["shannon", "simpson", "observed_features"]],
        ignore_index=True,
    )
    if not omnibus.empty:
        omnibus.to_csv(outdir / "alpha_omnibus_statistics.tsv", sep="\t", index=False)

    pairwise = pd.concat(
        [alpha_pairwise(alpha, metric, design["treatment_col"], design["time_col"], design["random_col"])
         for metric in ["shannon", "simpson", "observed_features"]],
        ignore_index=True,
    )
    if not pairwise.empty:
        pairwise.to_csv(outdir / "alpha_pairwise_contrasts.tsv", sep="\t", index=False)

    # Beta diversity.
    beta_input = rel
    for metric in [m.strip() for m in args.beta_metrics.split(",") if m.strip()]:
        dist = beta_distance(beta_input, metric)
        dist.to_csv(outdir / f"beta_distance_{metric}.tsv", sep="\t")
        coords = pcoa(dist)
        coords = coords.join(meta[[design["treatment_col"], design["time_col"]]], how="left")
        coords.to_csv(outdir / f"beta_pcoa_coordinates_{metric}.tsv", sep="\t")
        if {"PC1", "PC2"}.issubset(coords.columns):
            save_pcoa_plot(coords, meta, design["treatment_col"], design["time_col"], outdir / ".")
            # Avoid overwriting different metric plots.
            src = outdir / f"pcoa_{design['treatment_col']}_x_{design['time_col']}.png"
            dst = outdir / f"pcoa_{metric}_{design['treatment_col']}_x_{design['time_col']}.png"
            if src.exists():
                src.replace(dst)

        permanova_time_stratified(
            dist, meta, design["treatment_col"], design["time_col"], args.permutations, args.seed, outdir
        ).to_csv(outdir / f"beta_permanova_treatment_within_phase_{metric}.tsv", sep="\t", index=False)

        paired_time_permanova(
            dist, meta, design["treatment_col"], design["time_col"], design["random_col"],
            args.permutations, args.seed, outdir, metric
        )
        paired_change_magnitude(
            dist, meta, design["treatment_col"], design["time_col"], design["random_col"],
            args.permutations, args.seed, outdir, metric
        )

        if not args.no_exploratory_factorial_permanova:
            exploratory_full_factorial_permanova(
                dist, meta, design["treatment_col"], design["time_col"], args.permutations, args.seed, outdir
            ).to_csv(outdir / f"beta_permanova_exploratory_composite_{metric}.tsv", sep="\t", index=False)

    # Optional taxon-level association tests.
    if args.taxa_stats:
        taxon_statistics(table, meta, design, outdir, args.min_prevalence, args.min_mean_relative)

    print("Analysis complete.")
    print(f"Study: {design['study']}")
    print(f"Samples analysed: {table.shape[1]}")
    print(f"Treatment column: {design['treatment_col']}")
    print(f"Time/phase column: {design['time_col']}")
    print(f"Random/repeated-measures column: {design['random_col'] or 'NONE'}")
    print(f"Results: {outdir.resolve()}")


if __name__ == "__main__":
    main()
