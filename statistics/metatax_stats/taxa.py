"""Differential abundance of taxa: CLR transform + the same repeated-measures REML model as alpha diversity."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .contrasts import adjust_p, contrast_family
from .models import fit_model, wald_term

LN2 = np.log(2.0)


def clr_transform(counts: pd.DataFrame, pseudocount: float, zero_method: str = "multiplicative") -> pd.DataFrame:
    """CLR transform relative abundances using multiplicative zero replacement by default.

    Pseudocount mode remains available for sensitivity comparisons, but adds
    count-scale constants and therefore depends more directly on library size.
    """
    x = counts.to_numpy(float)
    totals = x.sum(axis=1, keepdims=True)
    if (totals <= 0).any():
        raise ValueError("CLR transform received a zero-depth sample.")
    p = x / totals
    if zero_method == "multiplicative":
        replaced = p.copy()
        for i, row in enumerate(p):
            zero = row == 0
            if zero.any():
                positive = row[~zero]
                delta = min(float(positive.min()) / 2.0, 0.5 / int(zero.sum()))
                replaced[i, zero] = delta
                replaced[i, ~zero] *= (1.0 - delta * zero.sum()) / positive.sum()
        xlog = np.log(replaced)
    elif zero_method == "pseudocount":
        if pseudocount <= 0:
            raise ValueError("CLR pseudocount must be positive.")
        xlog = np.log(x + pseudocount)
    else:
        raise ValueError(f"Unknown CLR zero method '{zero_method}'.")
    return pd.DataFrame(xlog - xlog.mean(1, keepdims=True), index=counts.index, columns=counts.columns)


def run_taxa(counts: pd.DataFrame, unresolved: pd.Series, frame, fixed, rank, pseudocount: float,
             min_prevalence: float, min_mean_rel: float, include_unclassified: bool, notes, scope: str,
             zero_method: str = "multiplicative", min_unit_prevalence: float = 0.0) -> dict:
    """counts: samples x taxa (unrarefied, aggregated to the level). unresolved: bool Series by taxon."""
    idx = frame.df.index
    counts = counts.loc[idx]
    rel = counts.div(counts.sum(1).replace(0, np.nan), axis=0).fillna(0)
    clr = clr_transform(counts, pseudocount, zero_method)
    prevalence = (counts > 0).mean(0)
    mean_rel = rel.mean(0)
    detected = counts > 0
    support_meta = frame.df.loc[counts.index, ["group", "upper"]].copy()
    support_meta["unit"] = support_meta["upper"].where(support_meta["upper"].notna(), pd.Series(counts.index, index=counts.index))
    support = detected.copy()
    support.insert(0, "unit", support_meta["unit"].to_numpy())
    support.insert(0, "group", support_meta["group"].to_numpy())
    unit_detection = support.groupby(["group", "unit"], sort=False).max()
    group_support = {}
    for group in frame.group_levels:
        group_units = unit_detection.loc[group] if group in unit_detection.index.get_level_values("group") else pd.DataFrame()
        n_units = int(support_meta.loc[support_meta["group"] == group, "unit"].nunique())
        n_present = group_units.sum(0).astype(int) if n_units else pd.Series(0, index=counts.columns, dtype=int)
        group_support[group] = (n_present, n_units)
    all_units_present = unit_detection.sum(0).astype(int)
    all_units_total = int(sum(total for _, total in group_support.values()))
    unit_prev_by_group = {
        group: (present / total if total else pd.Series(0.0, index=counts.columns))
        for group, (present, total) in group_support.items()
    }
    max_group_unit_prevalence = pd.concat(unit_prev_by_group, axis=1).max(axis=1)
    support_details = []
    for taxon in counts.columns:
        parts = []
        for group, (present, total) in group_support.items():
            fraction = float(present[taxon] / total) if total else 0.0
            parts.append(f"{group}={int(present[taxon])}/{total} ({fraction:.3f})")
        support_details.append("; ".join(parts))
    keep = ((prevalence >= min_prevalence) & (mean_rel >= min_mean_rel) &
            (max_group_unit_prevalence >= min_unit_prevalence))
    if not include_unclassified:
        keep &= ~unresolved.reindex(counts.columns).fillna(False).to_numpy()
    filter_reasons = []
    for tx in counts.columns:
        failed = []
        if prevalence[tx] < min_prevalence:
            failed.append("sample prevalence")
        if mean_rel[tx] < min_mean_rel:
            failed.append("mean relative abundance")
        if max_group_unit_prevalence[tx] < min_unit_prevalence:
            failed.append("unit prevalence")
        if not include_unclassified and bool(unresolved.reindex([tx]).fillna(False).iloc[0]):
            failed.append("unclassified")
        filter_reasons.append("; ".join(failed))
    summary = pd.DataFrame({
        "taxon": counts.columns,
        "prevalence": prevalence.values,
        "mean_relative_abundance": mean_rel.values,
        "units_present": all_units_present.reindex(counts.columns).to_numpy(dtype=int),
        "units_total": all_units_total,
        "unit_prevalence": (all_units_present / all_units_total if all_units_total else 0).reindex(counts.columns).to_numpy(),
        "max_group_unit_prevalence": max_group_unit_prevalence.reindex(counts.columns).to_numpy(),
        "unit_support_by_group": support_details,
        "unresolved": unresolved.reindex(counts.columns).fillna(False).values,
        "tested": keep.values,
        "filter_reason": filter_reasons,
    })
    for group, (present, total) in group_support.items():
        summary[f"units_present[{group}]"] = present.reindex(counts.columns).to_numpy(dtype=int)
        summary[f"units_total[{group}]"] = total
        summary[f"unit_prevalence[{group}]"] = unit_prev_by_group[group].reindex(counts.columns).to_numpy()
    tested = list(counts.columns[keep.to_numpy()])
    out = dict(summary=summary, clr=clr, rel=rel)
    if not tested:
        notes.warn(scope, "No taxa passed the prevalence/abundance filters; no taxon models were fitted.")
        out.update(omnibus=pd.DataFrame(), contrasts=pd.DataFrame(), failures=pd.DataFrame())
        return out
    om, ct, fails = [], [], []
    for tx in tested:
        fit = fit_model(clr[tx].to_numpy(float), fixed, rank)
        if not fit.ok:
            fails.append(dict(taxon=tx, reason=fit.note))
            continue
        for term in fixed.terms:
            w = wald_term(fit, term)
            w["taxon"] = tx
            om.append(w)
        fam = contrast_family(fit, frame, mode="reference")
        if not fam.empty:
            fam.insert(0, "taxon", tx)
            ct.append(fam)
    omd = pd.DataFrame(om)
    ctd = pd.concat(ct, ignore_index=True) if ct else pd.DataFrame()
    support_cols = ["taxon", "units_present", "units_total", "unit_prevalence",
                    "max_group_unit_prevalence", "unit_support_by_group"]
    if not omd.empty:
        omd = omd.merge(summary[["taxon", "prevalence", "mean_relative_abundance", *support_cols[1:]]], on="taxon")
        omd["q_value"] = np.nan
        for _, ix in omd.groupby("term").groups.items():
            omd.loc[ix, "q_value"] = adjust_p(omd.loc[ix, "p_value"].to_numpy(float), "fdr_bh")
    if not ctd.empty:
        for col in ("estimate", "se", "ci_low", "ci_high"):
            ctd[f"log2_{col}"] = ctd[col] / LN2
        ctd["q_value"] = np.nan
        for _, ix in ctd.groupby(["family", "level_a", "level_b", "time", "group"], dropna=False).groups.items():
            ctd.loc[ix, "q_value"] = adjust_p(ctd.loc[ix, "p_value"].to_numpy(float), "fdr_bh")
        ctd = ctd.merge(summary[["taxon", "prevalence", "mean_relative_abundance", *support_cols[1:]]], on="taxon")
        method = (f"multiplicative zero replacement (delta=min positive proportion/2)" if zero_method == "multiplicative"
                  else f"pseudocount {pseudocount:g}")
        ctd["effect_scale"] = f"log2 CLR difference ({method}); level_a minus level_b"
    out.update(omnibus=omd, contrasts=ctd, failures=pd.DataFrame(fails))
    return out
