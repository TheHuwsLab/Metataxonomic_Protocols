#!/usr/bin/env python3
"""
QIIME2 taxonomy statistics for the RUSITEC and Piglet Cecum studies.

The input is normally a QIIME level-7 taxonomy table. The toolkit can aggregate
that table to kingdom, phylum, class, order, family, genus, or species, or run
the complete analysis at every recoverable rank with --taxonomic-level all.

v4 changes (normalisation/provenance fixes): sum-coded main effects (hand-built design matrix, no patsy),
REML + Wald-F with containment df, GEE/OLS fallbacks reported in a `note` column,
unclassified ranks labelled by lineage instead of "__", robust recognition of unknown/uncultured/unidentified taxonomy,
explicit resolved-vs-unresolved audit, exact within-vessel swap PERMANOVA with min attainable p reported,
validation run after exclusions,
rarefaction (--rarefy-depth), no silent taxon drops, metric-suffixed output files.

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

Normalisation / transformation outputs
-------------------------------------
The toolkit explicitly names every abundance representation so the input to each
analysis is unambiguous:

* `rarefied_relative_abundance.tsv`: sample-wise relative abundance after the
  optional common-depth rarefaction; used for Bray-Curtis and abundance plots.
* `rarefied_presence_absence.tsv`: 0/1 presence-absence derived from the same
  rarefied table; used for Jaccard.
* `unrarefied_relative_abundance.tsv`: descriptive relative abundance from the
  original, unrarefied counts.
* `unrarefied_clr.tsv`: centred log-ratio transformed unrarefied counts used for
  exploratory per-taxon mixed-model statistics.

Shannon, Simpson and observed features are calculated from the rarefied count
table; Shannon/Simpson internally convert each sample to proportions, while
observed features counts non-zero taxa. A `normalisation_manifest.tsv` records
these transformations for every taxonomic level.
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
from scipy.spatial.distance import pdist, squareform
from scipy.stats import chi2, f as f_dist, kruskal, mannwhitneyu, norm, t as t_dist, ttest_ind
from statsmodels.stats.multitest import multipletests


# -----------------------------------------------------------------------------
# Taxonomic resolution handling
# -----------------------------------------------------------------------------

TAXONOMIC_LEVELS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
TAXON_PREFIXES = {
    "k": "kingdom",
    "d": "kingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}


_BLANK = re.compile(r"^_*$")
# Taxonomy classifiers often use many variants of unresolved labels. Exact-match
# handling alone is unsafe: labels such as "uncultured_bacterium" must not become
# headline biological taxa.
_UNRESOLVED_PATTERNS = (
    re.compile(r"^(?:unknown|unclassified|uncultured|unidentified|unassigned|nan)$", re.I),
    re.compile(r"^(?:unknown|unclassified|uncultured|unidentified|unassigned|metagenome|environmental|no[_ -]?hit|ambiguous)(?:[_ -].*)$", re.I),
    re.compile(r"^sp\.?(?:$|[_ -].*)", re.I),
)
UNCLASSIFIED_SUFFIX = "_unclassified"


def unresolved_reason(value: str | None) -> str | None:
    if value is None:
        return "missing_or_blank"
    text = str(value).strip()
    if not text or _BLANK.match(text):
        return "missing_or_blank"
    low = text.lower()
    for pattern in _UNRESOLVED_PATTERNS:
        if pattern.match(low):
            return "unresolved_label"
    return None


def is_unresolved_value(value: str | None) -> bool:
    return unresolved_reason(value) is not None


def parse_taxonomy_string(taxon: str) -> dict[str, str | None]:
    """Parse a semicolon-delimited QIIME/GTDB taxonomy string into rank labels.

    Empty positional fields are preserved for legacy non-prefixed taxonomy strings.
    Explicit unresolved labels (e.g. unknown, unclassified, uncultured_bacterium)
    are converted to None so they cannot be mistaken for named biological taxa.
    """
    parts = str(taxon).split(";")
    parsed = {rank: None for rank in TAXONOMIC_LEVELS}
    for i, raw_part in enumerate(parts):
        part = raw_part.strip()
        if i >= len(TAXONOMIC_LEVELS) and not part:
            continue
        m = re.match(r"^([dkpcofgs])__(.*)$", part, flags=re.IGNORECASE) if part else None
        if m:
            rank = TAXON_PREFIXES[m.group(1).lower()]
            value = m.group(2).strip()
        elif i < len(TAXONOMIC_LEVELS):
            rank, value = TAXONOMIC_LEVELS[i], part
        else:
            continue
        if unresolved_reason(value) is not None:
            parsed[rank] = None
        else:
            parsed[rank] = value
    return parsed


def is_unclassified_label(label: str) -> bool:
    text = str(label).strip()
    if text == "Unassigned" or text.endswith(UNCLASSIFIED_SUFFIX):
        return True
    # This also catches legacy output labels such as uncultured_bacterium if a
    # table generated by an earlier toolkit version is re-used as input.
    return is_unresolved_value(text)


def taxonomy_rank_mapping(table: pd.DataFrame, level: str) -> pd.DataFrame:
    """Build an auditable mapping from each input taxonomy to the selected rank.

    Unresolved entities are represented as lineage-specific *_unclassified bins;
    they remain in community-level compositional analyses as unresolved signal but
    are explicitly flagged and are excluded from headline "top taxa" and taxon DA
    results by default.
    """
    level = level.lower()
    if level not in TAXONOMIC_LEVELS:
        raise ValueError(f"Unsupported taxonomic level: {level}")
    shallower = TAXONOMIC_LEVELS[: TAXONOMIC_LEVELS.index(level)]
    rows = []
    for source in table.index.astype(str):
        parsed = parse_taxonomy_string(source)
        resolved = parsed[level]
        if resolved is not None:
            label, status, reason = resolved, "resolved", ""
            deepest = resolved
        else:
            deepest = next((parsed[r] for r in reversed(shallower) if parsed[r]), None)
            label = f"{deepest}{UNCLASSIFIED_SUFFIX}" if deepest else "Unassigned"
            status = "unresolved"
            reason = f"{level}:{unresolved_reason(None) if parsed[level] is None else 'unresolved_label'}"
        rows.append({
            "source_taxonomy": source,
            "taxonomic_level": level,
            "aggregated_taxon": label,
            "taxon_status": status,
            "unresolved_reason": reason,
            "deepest_classified_ancestor": deepest if deepest is not None else "",
        })
    return pd.DataFrame(rows)


def aggregate_taxonomic_level(table: pd.DataFrame, level: str) -> pd.DataFrame:
    """Aggregate a level-7 table to one standard rank without merging unrelated unknowns."""
    mapping = taxonomy_rank_mapping(table, level)
    labels = mapping["aggregated_taxon"].tolist()
    out = table.copy()
    out.index = labels
    return out.groupby(out.index).sum()


def taxonomy_levels_available(table: pd.DataFrame) -> list[str]:
    availability = {rank: False for rank in TAXONOMIC_LEVELS}
    for taxon in table.index:
        for rank, value in parse_taxonomy_string(taxon).items():
            if value is not None:
                availability[rank] = True
    return [rank for rank in TAXONOMIC_LEVELS if availability[rank]]


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

    # Convert abundance values only after orientation. Do not silently convert
    # malformed/non-numeric abundance cells to zero: that could manufacture
    # absences and materially change diversity and differential-abundance results.
    numeric = out.apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any():
        bad = []
        for col in numeric.columns:
            mask = numeric[col].isna() & out[col].notna()
            for idx in out.index[mask][:3]:
                bad.append((str(idx), str(col), str(out.loc[idx, col])))
                if len(bad) >= 10:
                    break
            if len(bad) >= 10:
                break
        raise SystemExit(f"Non-numeric abundance values detected; refusing to convert them to zero. Examples: {bad}")
    out = numeric.astype(float)
    if out.columns.duplicated().any():
        dup = out.columns[out.columns.duplicated()].astype(str).tolist()
        raise SystemExit(f"Duplicate sample IDs/columns in abundance table: {dup[:10]}")
    # Duplicate feature/taxonomy rows are safe to combine because they represent
    # the same labelled feature in the exported count table.
    out = out.groupby(out.index).sum()
    return out


def validate_count_table(table: pd.DataFrame) -> None:
    """Fail loudly on data that are not non-negative integer counts."""
    arr = table.to_numpy(dtype=float)
    if not np.isfinite(arr).all():
        raise SystemExit("Abundance table contains non-finite values (NaN/inf) after parsing.")
    if (arr < 0).any():
        bad = int((arr < 0).sum())
        raise SystemExit(f"Abundance table contains {bad} negative values; count tables must be non-negative.")
    fractional = np.abs(arr - np.rint(arr)) > 1e-8
    if fractional.any():
        examples = arr[fractional][:5].tolist()
        raise SystemExit(
            "Input abundance table contains non-integer values. This toolkit expects raw integer counts "
            f"for rarefaction/CLR; example values: {examples}."
        )


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


def alpha_metrics(table: pd.DataFrame, unresolved_labels: set[str] | None = None) -> pd.DataFrame:
    rows = []
    unresolved_labels = unresolved_labels or set()
    for sample in table.columns:
        col = table[sample]
        unresolved_mask = col.index.astype(str).isin(unresolved_labels)
        unresolved_reads = float(col.loc[unresolved_mask].sum())
        total = float(col.sum())
        rows.append(
            {
                "sample_id": sample,
                "observed_features": int((col > 0).sum()),
                "observed_resolved_taxa": int(((col > 0) & (~unresolved_mask)).sum()),
                "observed_unresolved_bins": int(((col > 0) & unresolved_mask).sum()),
                "unresolved_relative_abundance": unresolved_reads / total if total > 0 else np.nan,
                "shannon": shannon(col),
                "simpson": simpson(col),
                "read_count_or_total": total,
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


def pcoa(dist: pd.DataFrame, n_axes: int = 2) -> tuple[pd.DataFrame, pd.Series]:
    """Classical PCoA with variance explained from positive eigenvalues.

    Returns
    -------
    coords
        Sample coordinates on the requested principal coordinates.
    variance_explained
        Fraction of the positive eigenvalue sum represented by each axis.

    Negative eigenvalues can occur with non-Euclidean dissimilarities such as
    Bray-Curtis. They are not plotted here; axis percentages therefore describe
    the proportion of the positive-eigenvalue solution represented by each axis.
    """
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
    pos_eigvals = eigvals[positive]
    pos_eigvecs = eigvecs[:, positive]
    if len(pos_eigvals) == 0:
        raise ValueError("PCoA produced no positive eigenvalues.")

    total_positive = float(pos_eigvals.sum())
    k = min(n_axes, len(pos_eigvals))
    coords = pos_eigvecs[:, :k] * np.sqrt(pos_eigvals[:k])
    cols = [f"PC{i+1}" for i in range(k)]
    variance = pd.Series(
        pos_eigvals[:k] / total_positive,
        index=cols,
        name="variance_explained",
    )
    return pd.DataFrame(coords, index=dist.index, columns=cols), variance


def rarefy(table: pd.DataFrame, depth: int, seed: int) -> pd.DataFrame:
    """Subsample every sample (column) to `depth` reads without replacement."""
    rng = np.random.default_rng(seed)
    validate_count_table(table)
    counts = table.astype(np.int64)
    cols = {}
    for s in counts.columns:
        cols[s] = rng.multivariate_hypergeometric(counts[s].to_numpy(), int(depth))
    return pd.DataFrame(cols, index=table.index)


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def save_alpha_plots(alpha: pd.DataFrame, meta: pd.DataFrame, treatment_col: str, time_col: str, outdir: Path):
    merged = alpha.join(meta[[treatment_col, time_col]], how="inner")
    for metric in ["shannon", "simpson", "observed_features", "observed_resolved_taxa"]:
        plt.figure(figsize=(8, 6))
        # Treatment x time interaction displayed as separate boxes by timepoint.
        merged.boxplot(column=metric, by=[time_col, treatment_col], rot=45)
        plt.title(f"{metric.title()} by {treatment_col} and {time_col}")
        plt.suptitle("")
        plt.tight_layout()
        plt.savefig(outdir / f"alpha_{metric}_by_{treatment_col}_x_{time_col}.png", dpi=300)
        plt.close()

    # Time course: mean +/- SEM by treatment.
    for metric in ["shannon", "simpson", "observed_features", "observed_resolved_taxa"]:
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


def save_pcoa_plot(
    coords: pd.DataFrame,
    meta: pd.DataFrame,
    treatment_col: str,
    time_col: str,
    outdir: Path,
    metric: str,
    variance_explained: pd.Series | None = None,
    random_col: str | None = None,
    taxonomic_level: str | None = None,
) -> None:
    """Plot PCoA with treatment encoded by colour and time/phase by marker.

    For a repeated-measures design such as RUSITEC, samples belonging to the
    same Vessel are connected when there is exactly one observation at each of
    the two timepoints. This visualises the within-vessel trajectory without
    implying that the line itself is an inferential test.
    """
    if treatment_col in coords.columns and time_col in coords.columns:
        plot_df = coords.copy()
    else:
        plot_df = coords.join(meta[[treatment_col, time_col]], how="left")

    plot_df = plot_df.join(meta[[random_col]] if random_col and random_col not in plot_df.columns else pd.DataFrame(index=plot_df.index), how="left")

    treatments = list(pd.unique(plot_df[treatment_col].dropna().astype(str)))
    times = list(pd.unique(plot_df[time_col].dropna().astype(str)))

    cmap = plt.get_cmap("tab10")
    treatment_colours = {t: cmap(i % 10) for i, t in enumerate(treatments)}
    marker_cycle = ["o", "^", "s", "D", "P", "X", "v", "<", ">", "*", "h", "p"]
    time_markers = {t: marker_cycle[i % len(marker_cycle)] for i, t in enumerate(times)}

    fig, ax = plt.subplots(figsize=(9.5, 7.5))

    # Connect paired observations within the experimental unit.
    if random_col and random_col in plot_df.columns:
        for block, sub in plot_df.dropna(subset=[random_col, time_col]).groupby(random_col):
            sub = sub.copy()
            sub["__time"] = sub[time_col].astype(str)
            if sub.shape[0] == 2 and sub["__time"].nunique() == 2 and sub[treatment_col].nunique() == 1:
                sub = sub.sort_values("__time")
                ax.plot(
                    sub["PC1"],
                    sub["PC2"],
                    color="0.65",
                    linewidth=1.1,
                    alpha=0.65,
                    zorder=1,
                )

    # Plot observations by treatment colour and time/phase marker.
    for time in times:
        for treatment in treatments:
            sub = plot_df[
                (plot_df[time_col].astype(str) == time)
                & (plot_df[treatment_col].astype(str) == treatment)
            ]
            if sub.empty:
                continue
            ax.scatter(
                sub["PC1"],
                sub["PC2"],
                s=80,
                marker=time_markers[time],
                color=treatment_colours[treatment],
                edgecolor="black",
                linewidth=0.45,
                alpha=0.9,
                zorder=3,
            )

    # Separate legends avoid the unreadable 8-way treatment x time legend.
    treatment_handles = [
        plt.Line2D(
            [0], [0], marker="o", linestyle="", markersize=8,
            markerfacecolor=treatment_colours[t], markeredgecolor="black",
            markeredgewidth=0.45, label=t,
        )
        for t in treatments
    ]
    time_handles = [
        plt.Line2D(
            [0], [0], marker=time_markers[t], linestyle="", markersize=8,
            markerfacecolor="0.65", markeredgecolor="black", label=t,
        )
        for t in times
    ]

    leg1 = ax.legend(
        handles=treatment_handles,
        title=treatment_col,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=True,
    )
    ax.add_artist(leg1)
    ax.legend(
        handles=time_handles,
        title=time_col,
        loc="upper left",
        bbox_to_anchor=(1.02, 0.58),
        frameon=True,
    )

    xlab = "PCoA1"
    ylab = "PCoA2"
    if variance_explained is not None:
        if "PC1" in variance_explained.index:
            xlab += f" ({variance_explained['PC1'] * 100:.1f}%)"
        if "PC2" in variance_explained.index:
            ylab += f" ({variance_explained['PC2'] * 100:.1f}%)"
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

    level_text = f" — {taxonomic_level} level" if taxonomic_level else ""
    ax.set_title(f"{metric.replace('_', '-').title()} PCoA{level_text}")
    ax.axhline(0, color="0.90", linewidth=0.7, zorder=0)
    ax.axvline(0, color="0.90", linewidth=0.7, zorder=0)
    ax.grid(False)
    fig.tight_layout()
    suffix = f"_{taxonomic_level}" if taxonomic_level else ""
    fig.savefig(
        outdir / f"pcoa_{metric}_{treatment_col}_x_{time_col}{suffix}.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


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


def validate_design(meta: pd.DataFrame, design: dict, excluded: list | None = None) -> pd.DataFrame:
    """Check the metadata AFTER sample exclusions against the expected study design."""
    rows = []
    study = design["study"]
    tr = design["treatment_col"]
    tm = design["time_col"]
    rnd = design["random_col"]

    def add(check, status, detail):
        rows.append({"check": check, "status": status, "detail": detail})

    add("sample_count", "INFO", f"{len(meta)} samples retained for analysis")
    if excluded:
        add("excluded_samples", "CHECK", "; ".join(f"{s} ({why})" for s, why in excluded))
    add("treatments", "INFO", f"{meta[tr].nunique(dropna=True)} levels: {', '.join(map(str, meta[tr].dropna().unique()))}")
    add("time_or_phase", "INFO", f"{meta[tm].nunique(dropna=True)} levels: {', '.join(map(str, meta[tm].dropna().unique()))}")

    if study == "rusitec":
        vessels = meta[rnd].nunique(dropna=True) if rnd else 0
        add("RUSITEC_vessels", "PASS" if vessels == 16 else "CHECK", f"{vessels} unique vessels; expected 16 from design document")
        reps = meta.groupby(rnd).size() if rnd else pd.Series(dtype=int)
        bad = reps[reps != 2]
        ok = len(bad) == 0 and len(reps) == 16
        add("RUSITEC_repeated_timepoints", "PASS" if ok else "CHECK",
            "Each vessel has exactly two observations" if ok else f"Unexpected vessel replication (observations per vessel): {bad.to_dict()}")
        tr_by_vessel = meta.groupby(rnd)[tr].nunique() if rnd else pd.Series(dtype=int)
        const = len(tr_by_vessel) > 0 and (tr_by_vessel == 1).all()
        add("RUSITEC_treatment_constant_within_vessel", "PASS" if const else "CHECK",
            "Treatment is constant within vessel" if const else "At least one vessel has multiple treatment labels")
        n_vessels_by_treatment = meta.groupby(tr)[rnd].nunique().to_dict() if rnd else {}
        add("RUSITEC_treatment_replication", "PASS" if set(n_vessels_by_treatment.values()) == {4} else "CHECK",
            f"Vessels per treatment: {n_vessels_by_treatment}; expected 4 each")

    elif study == "piglet":
        add("Piglet_sample_count", "PASS" if len(meta) == 40 else "CHECK", f"{len(meta)} samples; design package describes 40 cecal samples")
        counts = meta.groupby([tr, tm]).size().to_dict()
        add("Piglet_treatment_x_phase_balance", "PASS" if set(counts.values()) == {5} else "CHECK",
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
# Alpha-diversity / per-taxon models:  y ~ Treatment * Time  (+ random intercept)
# -----------------------------------------------------------------------------
#
# The design matrix is built by hand with SUM (effect) coding so that
#   * the 'treatment' and 'time' terms are true main effects (averaged over the
#     other factor), not simple effects at the reference level, and
#   * nothing depends on patsy/formulaic (statsmodels >= 0.15 dropped patsy).

MODEL_LABELS = {
    "mixedlm": "MixedLM random-intercept (REML)",
    "gee": "GEE exchangeable (bias-reduced)",
    "ols_cluster": "OLS with cluster-robust covariance (fallback)",
    "ols": "OLS",
}


class Design:
    def __init__(self, treatment: pd.Series, time: pd.Series):
        tr = treatment.astype(str)
        tm = time.astype(str)
        self.tlevels = sorted(tr.unique())
        self.mlevels = sorted(tm.unique())
        self.X, self.terms = self._build(tr.to_numpy(), tm.to_numpy(), treatment.index)
        self.columns = list(self.X.columns)

    def _build(self, tr, tm, index):
        cols = {"Intercept": np.ones(len(tr))}
        tcols, mcols = {}, {}
        for lv in self.tlevels[:-1]:
            tcols[f"T[{lv}]"] = (tr == lv).astype(float) - (tr == self.tlevels[-1]).astype(float)
        for lv in self.mlevels[:-1]:
            mcols[f"M[{lv}]"] = (tm == lv).astype(float) - (tm == self.mlevels[-1]).astype(float)
        terms = {"treatment": list(tcols), "time": list(mcols), "treatment:time": []}
        cols.update(tcols)
        cols.update(mcols)
        for tn, tc in tcols.items():
            for mn, mc in mcols.items():
                cols[f"{tn}:{mn}"] = tc * mc
                terms["treatment:time"].append(f"{tn}:{mn}")
        return pd.DataFrame(cols, index=index), terms

    def row(self, treatment, time) -> np.ndarray:
        X, _ = self._build(np.array([str(treatment)]), np.array([str(time)]), pd.RangeIndex(1))
        return X.iloc[0].to_numpy(float)


ICC_BOUNDARY = 0.05


def _anova_icc(y, X: pd.DataFrame, groups) -> float:
    """Rough ICC of OLS residuals by cluster (one-way ANOVA estimator)."""
    try:
        r = pd.Series(np.asarray(sm.OLS(y, X).fit().resid, float))
        g = pd.Series(np.asarray(groups))
        n, k = len(r), g.nunique()
        if k < 2 or n <= k:
            return np.nan
        gm = r.groupby(g).transform("mean")
        msw = float(((r - gm) ** 2).sum()) / (n - k)
        msb = float((gm ** 2).sum()) / (k - 1)
        n0 = (n - (g.value_counts() ** 2).sum() / n) / (k - 1)
        var_b = (msb - msw) / n0
        return float(var_b / (var_b + msw)) if (var_b + msw) > 0 else np.nan
    except Exception:
        return np.nan


def fit_design(d: pd.DataFrame, has_random: bool, maxiter: int = 2000):
    """Fit __y ~ Treatment*Time. MixedLM (REML) -> GEE -> OLS fallback chain.

    Returns (result, model_type, Design, note). `note` records why a fallback
    happened so it is never silent.
    """
    design = Design(d["__treatment"], d["__time"])
    X = design.X
    y = d["__y"].to_numpy(float)
    note = ""
    if not np.isfinite(y).all():
        return None, "constant", design, "response contains non-finite values"
    if np.nanstd(y) <= 1e-12:
        return None, "constant", design, "response has zero variance; treatment/time effects are not estimable"
    if has_random:
        groups = d["__random"].to_numpy()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mod = sm.MixedLM(y, X, groups=groups)
                res = mod.fit(reml=True, method="lbfgs", maxiter=maxiter, disp=False)
                if not res.converged:
                    res = mod.fit(reml=True, method="powell", maxiter=maxiter, disp=False)
            k = len(res.fe_params)
            cov_ok = np.isfinite(np.asarray(res.cov_params(), float)[:k, :k]).all()
            if res.converged and cov_ok and np.isfinite(np.asarray(res.fe_params, float)).all():
                return res, "mixedlm", design, note
            note = "MixedLM did not converge"
        except Exception as exc:  # singular Hessian, rank deficiency, ...
            note = f"MixedLM failed ({type(exc).__name__}: {exc})"
        # Singular Hessian usually means the random-intercept variance is at the
        # boundary (~0), where the mixed model equals OLS. Use OLS only when the
        # ANOVA estimate of the intraclass correlation is negligible; else GEE.
        icc = _anova_icc(y, X, groups)
        if np.isfinite(icc) and icc < ICC_BOUNDARY:
            base = sm.OLS(y, X).fit()
            try:
                res = base.get_robustcov_results(cov_type="cluster", groups=groups, use_correction=True, df_correction=True)
                return res, "ols_cluster", design, (
                    note + f"; cluster variance ~0 (ANOVA ICC={icc:.3f}); used cluster-robust OLS rather than independence-assumption OLS")
            except Exception as exc:
                note += f"; cluster-robust covariance failed ({type(exc).__name__}: {exc})"
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                gee = sm.GEE(y, X, groups=groups, cov_struct=sm.cov_struct.Exchangeable(),
                             family=sm.families.Gaussian())
                try:
                    res = gee.fit(cov_type="bias_reduced")   # small-cluster-corrected sandwich
                except Exception:
                    res = gee.fit()
            return res, "gee", design, note + (
                "; fell back to GEE (bias-reduced covariance, Wald F with containment df; "
                "approximate - treat with caution)")
        except Exception as exc:
            note += f"; GEE failed ({type(exc).__name__}: {exc})"
            # Last resort with a supplied clustering variable: retain the cluster
            # dependence in the covariance rather than silently reverting to iid OLS.
            try:
                base = sm.OLS(y, X).fit()
                res = base.get_robustcov_results(cov_type="cluster", groups=groups, use_correction=True, df_correction=True)
                return res, "ols_cluster", design, note + "; used cluster-robust OLS fallback"
            except Exception as robust_exc:
                note += f"; cluster-robust OLS failed ({type(robust_exc).__name__}: {robust_exc})"
    res = sm.OLS(y, X).fit()
    return res, "ols", design, note + ("; no clustering variable supplied" if not has_random else "; WARNING: independence-assumption OLS fallback")


def prepare_model_data(data, metric, treatment_col, time_col, random_col):
    cols = [metric, treatment_col, time_col] + ([random_col] if random_col else [])
    d = data[cols].dropna().copy()
    d["__y"] = pd.to_numeric(d[metric], errors="coerce")
    d = d.dropna(subset=["__y"])
    d["__treatment"] = d[treatment_col].astype(str)
    d["__time"] = d[time_col].astype(str)
    if random_col:
        d["__random"] = d[random_col].astype(str)
    return d


def fit_mixed_model(data, metric, treatment_col, time_col, random_col):
    d = prepare_model_data(data, metric, treatment_col, time_col, random_col)
    if d.empty:
        return None, None, d, None, ""
    result, model_type, design, note = fit_design(d, bool(random_col))
    return result, model_type, d, design, note


def _fixed(result):
    """Fixed-effect estimates and their covariance block."""
    raw = result.fe_params if hasattr(result, "fe_params") else result.params
    params = np.asarray(raw, float).reshape(-1)
    cov = np.asarray(result.cov_params(), float)[: len(params), : len(params)]
    return params, cov


def containment_ddf(d: pd.DataFrame) -> dict:
    """Containment denominator df for a random-intercept split-plot design.

    Between-cluster term (treatment, when constant within cluster): G - T.
    Within-cluster terms (time, treatment:time): N - G - [(M-1) + (T-1)(M-1)].
    """
    n = len(d)
    g = d["__random"].nunique()
    t = d["__treatment"].nunique()
    m = d["__time"].nunique()
    within = max(n - g - ((m - 1) + (t - 1) * (m - 1)), 1)
    constant = bool((d.groupby("__random")["__treatment"].nunique() == 1).all())
    between = max(g - t, 1) if constant else within
    return {"treatment": between, "time": within, "treatment:time": within}


def wald_term_test(result, design: Design, term: str, ddf=None):
    """Wald test of one model term. Returns (stat, df, ddf, p, label).

    With a denominator df the statistic is reported as F (= chi2/df, referred to
    F(df, ddf)); otherwise as a chi-square.
    """
    cols = design.terms.get(term, [])
    if not cols:
        return np.nan, np.nan, np.nan, np.nan, "n/a"
    params, cov = _fixed(result)
    idx = [design.columns.index(c) for c in cols]
    b = params[idx]
    V = cov[np.ix_(idx, idx)]
    try:
        stat = float(b @ np.linalg.pinv(V) @ b)
        df = int(np.linalg.matrix_rank(V))
        if df == 0 or not np.isfinite(stat):
            return np.nan, np.nan, np.nan, np.nan, "n/a"
        if ddf and ddf > 0:
            F = stat / df
            return F, df, float(ddf), float(f_dist.sf(F, df, ddf)), "Wald F"
        return stat, df, np.nan, float(chi2.sf(stat, df)), "Wald chi2"
    except Exception:
        return np.nan, np.nan, np.nan, np.nan, "n/a"


def _ddf_for(result, model_type, d):
    if model_type in {"mixedlm", "gee"}:
        return containment_ddf(d)
    if model_type == "ols_cluster":
        # Robust cluster covariance is approximately referenced to the number of
        # independent experimental units, not the raw number of observations.
        g = d["__random"].nunique() if "__random" in d.columns else len(d)
        r = float(max(g - 1, 1))
        return {"treatment": r, "time": r, "treatment:time": r}
    if model_type == "ols":
        r = float(result.df_resid)
        return {"treatment": r, "time": r, "treatment:time": r}
    return {}


def alpha_omnibus(data, metric, treatment_col, time_col, random_col):
    result, model_type, d, design, note = fit_mixed_model(data, metric, treatment_col, time_col, random_col)
    if d.empty:
        return pd.DataFrame()
    if result is None and model_type == "constant":
        return pd.DataFrame([
            {"metric": metric, "term": term, "model": "Not estimable", "test": "not estimable",
             "test_statistic": np.nan, "df": np.nan, "ddf": np.nan, "p_value": np.nan, "note": note}
            for term in ["treatment", "time", "treatment:time"]
        ])
    if result is None:
        return pd.DataFrame()
    ddfs = _ddf_for(result, model_type, d)
    rows = []
    for term in ["treatment", "time", "treatment:time"]:
        stat, df, ddf, p, label = wald_term_test(result, design, term, ddfs.get(term))
        rows.append({
            "metric": metric, "term": term, "model": MODEL_LABELS.get(model_type, model_type),
            "test": label, "test_statistic": stat, "df": df, "ddf": ddf, "p_value": p, "note": note,
        })
    return pd.DataFrame(rows)


def fixed_contrast(result, design: Design, a: tuple, b: tuple, ddf=None):
    """Wald contrast between two (treatment, time) cells."""
    c = design.row(*a) - design.row(*b)
    params, cov = _fixed(result)
    est = float(c @ params)
    var = float(c @ cov @ c)
    if not np.isfinite(var) or var <= 0:
        return est, np.nan, np.nan
    se = math.sqrt(var)
    z = est / se
    p = 2 * t_dist.sf(abs(z), ddf) if ddf and ddf > 0 else 2 * norm.sf(abs(z))
    return est, se, float(p)


def _bh_within(out: pd.DataFrame, by: str, p: str = "p_value", q: str = "p_fdr") -> pd.DataFrame:
    out[q] = np.nan
    for _, idx in out.groupby(by).groups.items():
        valid = out.loc[idx, p].notna()
        if valid.any():
            vidx = out.loc[idx].index[valid.to_numpy()]
            out.loc[vidx, q] = multipletests(out.loc[vidx, p], method="fdr_bh")[1]
    return out


def alpha_pairwise(data, metric, treatment_col, time_col, random_col):
    """Model-based contrasts (MixedLM/GEE) with Welch-test fallback for OLS."""
    result, model_type, d, design, note = fit_mixed_model(data, metric, treatment_col, time_col, random_col)
    if result is None or d.empty:
        return pd.DataFrame()

    treatments = list(pd.unique(d["__treatment"]))
    times = list(pd.unique(d["__time"]))
    rows = []

    if model_type in {"mixedlm", "gee"}:
        ddfs = _ddf_for(result, model_type, d)
        label = f"{MODEL_LABELS[model_type]} Wald contrast"
        for time in times:
            for a, b in itertools.combinations(treatments, 2):
                est, se, p = fixed_contrast(result, design, (a, time), (b, time), ddfs.get("treatment"))
                rows.append({"metric": metric, "contrast_type": "treatment_within_time", "time": time,
                             "group_a": a, "group_b": b, "estimate_a_minus_b": est, "SE": se,
                             "p_value": p, "model": label, "note": note})
        for treatment in treatments:
            for a, b in itertools.combinations(times, 2):
                est, se, p = fixed_contrast(result, design, (treatment, a), (treatment, b), ddfs.get("time"))
                rows.append({"metric": metric, "contrast_type": "time_within_treatment", "time": f"{a} vs {b}",
                             "group_a": a, "group_b": b, "treatment": treatment, "estimate_a_minus_b": est,
                             "SE": se, "p_value": p, "model": label, "note": note})
    else:
        for time in times:
            sub = d[d["__time"] == time]
            for a, b in itertools.combinations(treatments, 2):
                xa = sub.loc[sub["__treatment"] == a, "__y"]
                xb = sub.loc[sub["__treatment"] == b, "__y"]
                if len(xa) > 1 and len(xb) > 1:
                    _, p = ttest_ind(xa, xb, equal_var=False, nan_policy="omit")
                    rows.append({"metric": metric, "contrast_type": "treatment_within_time", "time": time,
                                 "group_a": a, "group_b": b, "estimate_a_minus_b": xa.mean() - xb.mean(),
                                 "SE": np.nan, "p_value": float(p), "model": "Welch t-test fallback", "note": note})
        for treatment in treatments:
            sub = d[d["__treatment"] == treatment]
            for a, b in itertools.combinations(times, 2):
                xa = sub.loc[sub["__time"] == a, "__y"]
                xb = sub.loc[sub["__time"] == b, "__y"]
                if len(xa) > 1 and len(xb) > 1:
                    _, p = ttest_ind(xa, xb, equal_var=False, nan_policy="omit")
                    rows.append({"metric": metric, "contrast_type": "time_within_treatment", "time": f"{a} vs {b}",
                                 "group_a": a, "group_b": b, "treatment": treatment,
                                 "estimate_a_minus_b": xa.mean() - xb.mean(), "SE": np.nan,
                                 "p_value": float(p), "model": "Welch t-test fallback", "note": note})

    out = pd.DataFrame(rows)
    if not out.empty:
        out = _bh_within(out, "contrast_type")
    return out


# -----------------------------------------------------------------------------
# PERMANOVA (single-factor; used within fixed timepoint/phase strata)
# -----------------------------------------------------------------------------

EXACT_MAX_BLOCKS = 15  # enumerate all 2^n within-block swaps up to this many blocks


def _gower(dist_values: np.ndarray) -> np.ndarray:
    n = dist_values.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    return -0.5 * H @ (dist_values ** 2) @ H


def _pseudo_f(G: np.ndarray, codes: np.ndarray, g: int):
    """One-factor pseudo-F from a Gower-centred matrix and integer group codes."""
    n = G.shape[0]
    if g < 2 or n <= g:
        return np.nan, np.nan, np.nan
    Z = np.zeros((n, g))
    Z[np.arange(n), codes] = 1.0
    ssb = float(((Z * (G @ Z)).sum(axis=0) / Z.sum(axis=0)).sum())
    ssw = float(np.trace(G)) - ssb
    F = (ssb / (g - 1)) / (ssw / (n - g)) if ssw > 0 else np.nan
    return float(F), ssb, ssw


def permanova_pseudo_f(dist: pd.DataFrame, labels: pd.Series) -> tuple[float, float, float]:
    """One-factor PERMANOVA pseudo-F, no permutations performed here."""
    ids = list(dist.index)
    codes, uniq = pd.factorize(pd.Series(labels, index=ids).astype(str))
    return _pseudo_f(_gower(dist.loc[ids, ids].to_numpy(float)), codes, len(uniq))


def permuted_pvalue(observed: float, values) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if not len(arr):
        return np.nan
    tol = 1e-10 * max(abs(observed), 1.0)
    return float((1 + np.sum(arr >= observed - tol)) / (1 + len(arr)))


def permanova_time_stratified(dist, meta, treatment_col, time_col, permutations, seed, outdir=None):
    """Treatment PERMANOVA separately at each phase/timepoint.

    One observation per vessel enters each test, so unrestricted permutation
    among vessels within the selected timepoint is valid.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for time in pd.unique(meta[time_col].dropna().astype(str)):
        ids = meta.index[meta[time_col].astype(str) == time].intersection(dist.index)
        codes, uniq = pd.factorize(meta.loc[ids, treatment_col].astype(str))
        G = _gower(dist.loc[ids, ids].to_numpy(float))
        g = len(uniq)
        observed, ssb, ssw = _pseudo_f(G, codes, g)
        if not np.isfinite(observed):
            continue
        null = [_pseudo_f(G, rng.permutation(codes), g)[0] for _ in range(permutations)]
        rows.append({
            "time_or_phase": time, "term": treatment_col, "n": len(ids), "groups": g,
            "pseudo_F": observed, "SS_between": ssb, "SS_within": ssw,
            "permutations": permutations, "p_value": permuted_pvalue(observed, null),
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out = _bh_within(out.assign(_all=1), "_all").drop(columns="_all")
    return out


def _swap_masks(n_blocks: int, permutations: int, rng):
    """All 2^n within-block swap patterns when feasible, else random ones."""
    if n_blocks <= EXACT_MAX_BLOCKS:
        bits = (np.arange(2 ** n_blocks)[:, None] >> np.arange(n_blocks)) & 1
        return bits.astype(bool), True
    return rng.random((permutations, n_blocks)) < 0.5, False


def _paired_time_test(dist, meta, ids, time_col, block_col, permutations, rng):
    """Time-label PERMANOVA with labels swapped only within blocks (vessels/pens)."""
    ids = list(ids)
    pos = {i: k for k, i in enumerate(ids)}
    pairs = []
    for _, sub in meta.loc[ids].groupby(block_col):
        pairs.append([pos[i] for i in sub.index])
    pairs = np.asarray(pairs)
    codes, uniq = pd.factorize(meta.loc[ids, time_col].astype(str))
    G = _gower(dist.loc[ids, ids].to_numpy(float))
    observed, ssb, ssw = _pseudo_f(G, codes, len(uniq))
    if not np.isfinite(observed):
        return None
    masks, exact = _swap_masks(len(pairs), permutations, rng)
    null = np.empty(len(masks))
    for k, mask in enumerate(masks):
        perm = codes.copy()
        sw = pairs[mask]
        if len(sw):
            perm[sw[:, 0]] = codes[sw[:, 1]]
            perm[sw[:, 1]] = codes[sw[:, 0]]
        null[k] = _pseudo_f(G, perm, len(uniq))[0]
    tol = 1e-10 * max(abs(observed), 1.0)
    if exact:
        p = float(np.mean(null >= observed - tol))      # identity is one of the masks
        min_p = 2.0 / len(masks)                         # identity + all-swapped give the same F
    else:
        p = float((1 + np.sum(null >= observed - tol)) / (1 + len(null)))
        min_p = 1.0 / (1 + len(null))
    return {"n_blocks": len(pairs), "pseudo_F": observed, "SS_between": ssb, "SS_within": ssw,
            "arrangements": len(masks), "exact": exact, "min_attainable_p": min_p, "p_value": p}


def paired_time_permanova(dist, meta, treatment_col, time_col, block_col, permutations, seed, outdir, metric):
    """Time/phase PERMANOVA within each treatment, plus one pooled test.

    Permutation is exact (all 2^n swaps) for <= EXACT_MAX_BLOCKS blocks. With only
    4 vessels per treatment the smallest attainable p-value is 2/2^4 = 0.125, so
    per-treatment tests cannot be significant; the pooled row (all vessels,
    swaps within vessel) is the adequately powered time test.
    """
    if not block_col or block_col not in meta.columns:
        return pd.DataFrame()
    rng = np.random.default_rng(seed)

    def good_ids(sub_ids):
        m = meta.loc[sub_ids]
        keep = []
        for block, sub in m.groupby(block_col):
            counts = sub[time_col].astype(str).value_counts()
            if len(counts) == 2 and counts.max() == 1:
                keep.extend(sub.index)
        return keep

    rows = []
    for treatment in pd.unique(meta[treatment_col].dropna().astype(str)):
        ids = meta.index[meta[treatment_col].astype(str) == treatment].intersection(dist.index)
        keep = good_ids(ids)
        if len(keep) < 4:
            continue
        res = _paired_time_test(dist, meta, keep, time_col, block_col, permutations, rng)
        if res:
            rows.append({"treatment": treatment, "term": time_col, **res})
    out = pd.DataFrame(rows)
    if not out.empty:
        out = _bh_within(out.assign(_all=1), "_all").drop(columns="_all")

    all_ids = good_ids(meta.index.intersection(dist.index))
    if len(all_ids) >= 4:
        res = _paired_time_test(dist, meta, all_ids, time_col, block_col, permutations, rng)
        if res:
            pooled = pd.DataFrame([{"treatment": "ALL (pooled)", "term": time_col, **res}])
            out = pd.concat([out, pooled], ignore_index=True)
    if not out.empty:
        out["note"] = f"Labels swapped within {block_col} only; p_fdr covers the per-treatment rows only."
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


def exploratory_full_factorial_permanova(dist, meta, treatment_col, time_col, permutations, seed, outdir=None):
    """Exploratory composite Treatment|Time PERMANOVA (ignores repeated measures)."""
    ids = dist.index.intersection(meta.index)
    labels = meta.loc[ids, [treatment_col, time_col]].astype(str).agg(" | ".join, axis=1)
    codes, uniq = pd.factorize(labels)
    G = _gower(dist.loc[ids, ids].to_numpy(float))
    g = len(uniq)
    observed, ssb, ssw = _pseudo_f(G, codes, g)
    if not np.isfinite(observed):
        return pd.DataFrame()
    rng = np.random.default_rng(seed)
    null = [_pseudo_f(G, rng.permutation(codes), g)[0] for _ in range(permutations)]
    return pd.DataFrame([{
        "analysis": "exploratory composite Treatment | Time/Phase PERMANOVA",
        "n": len(ids), "groups": g, "pseudo_F": observed, "SS_between": ssb, "SS_within": ssw,
        "permutations": permutations, "p_value": permuted_pvalue(observed, null),
        "note": "Does not account for repeated Vessel/Pen/Animal structure; use as exploratory only.",
    }])


# -----------------------------------------------------------------------------
# Optional taxon-level CLR mixed models
# -----------------------------------------------------------------------------

def taxon_statistics(table, meta, design, outdir, min_prevalence, min_mean_relative, include_unclassified=False):
    """Exploratory per-taxon CLR models (same fallback chain as alpha) with BH correction.

    Failures are recorded in taxon_clr_failed.tsv instead of being dropped silently.
    Unclassified pseudo-taxa ("X_unclassified", "Unassigned") are excluded by default.
    """
    clr = clr_samples(table)
    rel = table.div(table.sum(axis=0).replace(0, np.nan), axis=1).fillna(0)
    prevalence = (table > 0).mean(axis=1)
    mean_rel = rel.mean(axis=1)
    keep = (prevalence >= min_prevalence) & (mean_rel >= min_mean_relative)
    if not include_unclassified:
        keep &= ~pd.Series([is_unclassified_label(i) for i in table.index], index=table.index)
    taxa = list(table.index[keep])

    tr, tm, rnd = design["treatment_col"], design["time_col"], design["random_col"]
    rows, failed = [], []
    for taxon in taxa:
        d = pd.DataFrame({"__y": clr.loc[taxon].to_numpy(float)}, index=clr.columns).join(meta, how="inner")
        d = d.dropna(subset=[tr, tm] + ([rnd] if rnd else []))
        if len(d) < 8:
            failed.append({"taxon": taxon, "reason": f"only {len(d)} usable samples"})
            continue
        d["__treatment"] = d[tr].astype(str)
        d["__time"] = d[tm].astype(str)
        if rnd:
            d["__random"] = d[rnd].astype(str)
        try:
            result, model_type, dsn, note = fit_design(d, bool(rnd), maxiter=500)
            ddfs = _ddf_for(result, model_type, d)
            p = {t: wald_term_test(result, dsn, t, ddfs.get(t))[3] for t in ["treatment", "time", "treatment:time"]}
        except Exception as exc:
            failed.append({"taxon": taxon, "reason": f"{type(exc).__name__}: {exc}"})
            continue
        rows.append({
            "taxon": taxon,
            "taxon_status": "unresolved" if is_unclassified_label(taxon) else "resolved",
            "prevalence_in_unrarefied_counts": float(prevalence[taxon]),
            "mean_unrarefied_relative_abundance": float(mean_rel[taxon]),
            "model": MODEL_LABELS.get(model_type, model_type),
            "p_treatment": p["treatment"], "p_time_or_phase": p["time"], "p_treatment_x_time": p["treatment:time"],
            "response_transformation": "unrarefied CLR (pseudocount=1)",
            "note": note,
        })

    columns = ["taxon", "taxon_status", "prevalence_in_unrarefied_counts", "mean_unrarefied_relative_abundance", "model",
               "p_treatment", "p_time_or_phase", "p_treatment_x_time",
               "q_treatment", "q_time_or_phase", "q_treatment_x_time", "response_transformation", "note"]
    out = pd.DataFrame(rows, columns=columns)
    if not out.empty:
        for col in ["p_treatment", "p_time_or_phase", "p_treatment_x_time"]:
            qcol = col.replace("p_", "q_")
            valid = out[col].notna()
            if valid.any():
                out.loc[valid, qcol] = multipletests(out.loc[valid, col], method="fdr_bh")[1]
        out = out.sort_values(["q_treatment_x_time", "q_treatment", "q_time_or_phase"], na_position="last")
    out.to_csv(outdir / "taxon_clr_mixed_model_results.tsv", sep="\t", index=False)
    if failed:
        pd.DataFrame(failed).to_csv(outdir / "taxon_clr_failed.tsv", sep="\t", index=False)
        print(f"taxon stats: {len(failed)} of {len(taxa)} taxa could not be fitted (see taxon_clr_failed.tsv)", file=sys.stderr)
    return out


# -----------------------------------------------------------------------------
# Taxon-level output: normalised tables, summaries and plots
# -----------------------------------------------------------------------------

def _safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "group"


def taxon_abundance_outputs(
    level_raw: pd.DataFrame,
    level_analysis: pd.DataFrame,
    meta: pd.DataFrame,
    design: dict,
    outdir: Path,
    top_n: int = 20,
    taxonomy_status: pd.Series | None = None,
) -> dict:
    """Write taxon-level relative-abundance tables and publication-oriented plots.

    `level_raw` is the unrarefied count table used for taxon DA/CLR modelling.
    `level_analysis` is the table actually used for alpha/beta analysis (usually
    rarefied to a common depth). Relative-abundance plots use this same rarefied
    relative-abundance representation so the visualisation corresponds to the
    diversity/ordination input.
    """
    tr, tm = design["treatment_col"], design["time_col"]

    # Explicitly named abundance representations.
    # `level_analysis` is normally the rarefied table used for community-level
    # analysis; `level_raw` is always the unrarefied count table.
    rel = normalise(level_analysis, "relative")
    rel.to_csv(outdir / "rarefied_relative_abundance.tsv", sep="\t")
    rel.T.to_csv(outdir / "rarefied_relative_abundance_samples_x_taxa.tsv", sep="\t")

    presence_absence = (level_analysis > 0).astype(np.int8)
    presence_absence.to_csv(outdir / "rarefied_presence_absence.tsv", sep="\t")
    presence_absence.T.to_csv(outdir / "rarefied_presence_absence_samples_x_taxa.tsv", sep="\t")

    raw_counts = level_raw.copy()
    raw_counts.to_csv(outdir / "unrarefied_counts.tsv", sep="\t")
    raw_rel = normalise(level_raw, "relative")
    raw_rel.to_csv(outdir / "unrarefied_relative_abundance.tsv", sep="\t")

    # CLR is only for the exploratory taxon-level inferential analysis. Keep the
    # exact transformed matrix so users can inspect/reproduce those results.
    unrarefied_clr = clr_samples(level_raw)
    unrarefied_clr.to_csv(outdir / "unrarefied_clr.tsv", sep="\t")
    unrarefied_clr.T.to_csv(outdir / "unrarefied_clr_samples_x_taxa.tsv", sep="\t")

    # Machine-readable provenance describing which representation feeds which
    # analysis. This avoids ambiguous terms such as simply "normalised table".
    manifest = pd.DataFrame([
        {"output": "unrarefied_counts.tsv", "source": "original input after sample QC/exclusions", "transformation": "none", "used_for": "reference/raw counts; CLR source"},
        {"output": "rarefied_relative_abundance.tsv", "source": "rarefied count table", "transformation": "sample-wise total-sum scaling", "used_for": "Bray-Curtis; relative-abundance plots"},
        {"output": "rarefied_presence_absence.tsv", "source": "rarefied count table", "transformation": "count > 0 -> 1 else 0", "used_for": "Jaccard; observed-features presence"},
        {"output": "unrarefied_relative_abundance.tsv", "source": "unrarefied count table", "transformation": "sample-wise total-sum scaling", "used_for": "descriptive unrarefied abundance summaries"},
        {"output": "unrarefied_clr.tsv", "source": "unrarefied count table", "transformation": "pseudocount=1 followed by centred log-ratio", "used_for": "exploratory per-taxon mixed-model statistics"},
        {"output": "alpha_diversity.tsv", "source": "rarefied count table", "transformation": "Shannon/Simpson use within-sample proportions; observed_features uses non-zero taxa", "used_for": "alpha diversity"},
    ])
    manifest.to_csv(outdir / "normalisation_manifest.tsv", sep="\t", index=False)

    # Human-readable companion file so a user can see exactly what is normalised
    # and what is not without having to inspect the code.
    with open(outdir / "normalisation_README.txt", "w", encoding="utf-8") as fh:
        fh.write(
            "NORMALISATION / TRANSFORMATION GUIDE\n"
            "==================================\n\n"
            "unrarefied_counts.tsv: original count table after sample QC/exclusions.\n"
            "rarefied_relative_abundance.tsv: rarefied counts converted to sample-wise proportions; used for Bray-Curtis and abundance plots.\n"
            "rarefied_presence_absence.tsv: rarefied counts converted to 0/1 presence-absence; used for Jaccard.\n"
            "unrarefied_relative_abundance.tsv: original unrarefied counts converted to sample-wise proportions; descriptive only.\n"
            "unrarefied_clr.tsv: unrarefied counts + pseudocount 1, followed by CLR; used for exploratory per-taxon models.\n"
            "alpha_diversity.tsv: Shannon/Simpson are calculated from the rarefied count table after within-sample conversion to proportions; observed_features counts all non-zero bins, while observed_resolved_taxa excludes unresolved bins.\n"
            "taxonomy_rank_mapping.tsv / taxonomy_rank_status.tsv: audit which input lineages were resolved at the requested rank. Unresolved bins remain in community-level composition but are never treated as named taxa in headline top-taxa outputs.\n"
            "taxonomy_qc_summary.tsv and unresolved_abundance_by_group.tsv: quantify how much community abundance is unresolved at the selected rank.\n"
        )

    prevalence = (level_analysis > 0).mean(axis=1)
    mean_rel = rel.mean(axis=1)
    if taxonomy_status is None:
        status = pd.Series({tax: ("unresolved" if is_unclassified_label(tax) else "resolved") for tax in rel.index})
    else:
        status = taxonomy_status.reindex(rel.index).fillna("resolved").astype(str)
    prevalence_table = pd.DataFrame({
        "taxon": rel.index,
        "taxon_status": status.reindex(rel.index).values,
        "prevalence_in_analysis_table": prevalence.reindex(rel.index).values,
        "mean_rarefied_relative_abundance": mean_rel.reindex(rel.index).values,
    }).sort_values(["taxon_status", "mean_rarefied_relative_abundance", "prevalence_in_analysis_table"], ascending=[True, False, False])
    prevalence_table.to_csv(outdir / "taxa_prevalence_and_mean_abundance.tsv", sep="\t", index=False)

    # Explicit unresolved-taxonomy burden. This is deliberately kept separate
    # from named biological taxa so an unresolved category can never silently
    # become the headline "largest taxon".
    unresolved = prevalence_table[prevalence_table["taxon_status"] == "unresolved"].copy()
    unresolved.to_csv(outdir / "unresolved_taxa_summary.tsv", sep="\t", index=False)
    unresolved_mean_by_sample = rel.loc[rel.index.intersection(unresolved["taxon"])] .sum(axis=0) if not unresolved.empty else pd.Series(0.0, index=rel.columns)
    taxonomy_qc = pd.DataFrame([{
        "n_resolved_taxa_bins": int((status == "resolved").sum()),
        "n_unresolved_taxa_bins": int((status == "unresolved").sum()),
        "mean_unresolved_relative_abundance_across_samples": float(unresolved_mean_by_sample.mean()),
        "max_unresolved_relative_abundance_single_sample": float(unresolved_mean_by_sample.max()),
        "median_unresolved_relative_abundance_across_samples": float(unresolved_mean_by_sample.median()),
        "interpretation": "Unresolved bins are retained in community composition but are not counted as named taxa in headline top-taxa or taxon-CLR results unless explicitly requested.",
    }])
    taxonomy_qc.to_csv(outdir / "taxonomy_qc_summary.tsv", sep="\t", index=False)

    # Sample metadata / group labels.
    plot_meta = meta[[tr, tm]].copy()
    plot_meta["group"] = plot_meta[tm].astype(str) + " | " + plot_meta[tr].astype(str)
    group_order = list(pd.unique(plot_meta["group"]))

    # Group-level mean / SD / SEM / prevalence for every taxon.
    records = []
    long = rel.T.join(plot_meta, how="inner")
    taxa = list(rel.index)
    for group, sub in long.groupby("group", sort=False):
        abundance = sub[taxa]
        for taxon in taxa:
            vals = abundance[taxon].to_numpy(float)
            n = len(vals)
            sd = float(np.std(vals, ddof=1)) if n > 1 else np.nan
            records.append({
                "group": group,
                tr: str(sub[tr].iloc[0]),
                tm: str(sub[tm].iloc[0]),
                "n": n,
                "taxon": taxon,
                "mean_rarefied_relative_abundance": float(np.mean(vals)),
                "sd_rarefied_relative_abundance": sd,
                "sem_rarefied_relative_abundance": sd / math.sqrt(n) if n > 1 else np.nan,
                "prevalence": float(np.mean(vals > 0)),
            })
    group_summary = pd.DataFrame(records)
    if not group_summary.empty:
        group_summary.to_csv(outdir / "taxa_rarefied_relative_abundance_group_summary.tsv", sep="\t", index=False)
        unresolved_group = (group_summary[group_summary["taxon"].isin(unresolved["taxon"]) ]
                           .groupby(["group", tr, tm], as_index=False)
                           ["mean_rarefied_relative_abundance"].sum()
                           .rename(columns={"mean_rarefied_relative_abundance": "mean_unresolved_relative_abundance"}))
        unresolved_group.to_csv(outdir / "unresolved_abundance_by_group.tsv", sep="\t", index=False)

    # Headline top taxa are ALWAYS resolved biological taxa. Unresolved bins are
    # reported separately and never substituted into the top-taxa list.
    ranked = prevalence_table.loc[prevalence_table["taxon_status"] == "resolved"]
    top_taxa = ranked.head(max(1, top_n))["taxon"].tolist()
    pd.DataFrame({"rank": np.arange(1, len(top_taxa) + 1), "taxon": top_taxa}).to_csv(
        outdir / "top_taxa.tsv", sep="\t", index=False
    )

    if not group_summary.empty and top_taxa:
        means = group_summary[group_summary["taxon"].isin(top_taxa)].pivot(
            index="group", columns="taxon", values="mean_rarefied_relative_abundance"
        ).reindex(group_order)
        means = means.fillna(0.0)
        ordered_cols = [t for t in top_taxa if t in means.columns]
        means = means[ordered_cols]

        # Stacked barplot: named resolved top taxa + unresolved abundance + other resolved.
        unresolved_taxa = [t for t in rel.index if status.get(t, "resolved") == "unresolved"]
        unresolved_means = group_summary[group_summary["taxon"].isin(unresolved_taxa)].groupby("group")["mean_rarefied_relative_abundance"].sum()
        unresolved_means = unresolved_means.reindex(group_order).fillna(0.0)
        other_resolved = (1.0 - means.sum(axis=1) - unresolved_means).clip(lower=0.0)
        stack = means.copy()
        stack["Other resolved"] = other_resolved
        stack["Unresolved at selected rank"] = unresolved_means
        fig, ax = plt.subplots(figsize=(max(9, len(group_order) * 1.15), 7))
        bottom = np.zeros(len(stack))
        cmap = plt.get_cmap("tab20")
        for i, taxon in enumerate(ordered_cols):
            vals = stack[taxon].to_numpy()
            ax.bar(np.arange(len(stack)), vals, bottom=bottom, label=str(taxon),
                   color=cmap(i % 20), width=0.78)
            bottom += vals
        ax.bar(np.arange(len(stack)), stack["Other resolved"].to_numpy(), bottom=bottom,
               label="Other resolved", color="0.75", width=0.78)
        bottom += stack["Other resolved"].to_numpy()
        ax.bar(np.arange(len(stack)), stack["Unresolved at selected rank"].to_numpy(), bottom=bottom,
               label="Unresolved at selected rank", color="0.45", width=0.78)
        ax.set_xticks(np.arange(len(stack)))
        ax.set_xticklabels(stack.index, rotation=45, ha="right")
        ax.set_ylabel("Mean rarefied relative abundance")
        ax.set_xlabel(f"{tm} | {tr}")
        ax.set_title(f"Top {len(ordered_cols)} taxa by mean rarefied relative abundance")
        ax.legend(title="Taxon", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
        fig.tight_layout()
        fig.savefig(outdir / "taxa_top_abundance_stacked_bar.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

        # Heatmap of group means. Log10 scaling retains low-abundance taxa while
        # avoiding arbitrary z-scoring of compositional abundances.
        heat = np.log10(means.T + 1e-6)
        fig, ax = plt.subplots(figsize=(max(8, len(group_order) * 1.0), max(7, len(ordered_cols) * 0.35)))
        im = ax.imshow(heat.to_numpy(float), aspect="auto", interpolation="nearest", cmap="viridis")
        ax.set_yticks(np.arange(len(heat.index)))
        ax.set_yticklabels(heat.index, fontsize=8)
        ax.set_xticks(np.arange(len(heat.columns)))
        ax.set_xticklabels(heat.columns, rotation=45, ha="right", fontsize=8)
        ax.set_xlabel(f"{tm} | {tr}")
        ax.set_ylabel("Taxon")
        ax.set_title(f"Top {len(ordered_cols)} taxa — log10 mean rarefied relative abundance")
        fig.colorbar(im, ax=ax, label="log10(rarefied relative abundance + 1e-6)")
        fig.tight_layout()
        fig.savefig(outdir / "taxa_top_abundance_heatmap.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

        # Boxplots for the top taxa. This produces one figure with one panel per
        # taxon without requiring seaborn.
        n = len(ordered_cols)
        ncols = 2 if n > 1 else 1
        nrows = int(math.ceil(n / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(12, max(5, 3.2 * nrows)), squeeze=False)
        axes_flat = axes.ravel()
        box_labels = [f"{str(x)}" for x in group_order]
        group_arrays = [long.loc[long["group"] == g, taxon].to_numpy(float) for g in group_order for taxon in []]
        for i, taxon in enumerate(ordered_cols):
            ax = axes_flat[i]
            data = [long.loc[long["group"] == g, taxon].to_numpy(float) for g in group_order]
            ax.boxplot(data, tick_labels=box_labels, showfliers=True)
            ax.set_title(str(taxon), fontsize=9)
            ax.set_ylabel("Rarefied relative abundance")
            ax.tick_params(axis="x", rotation=55, labelsize=7)
        for ax in axes_flat[n:]:
            ax.axis("off")
        fig.suptitle(f"Top {len(ordered_cols)} taxa — sample-level rarefied relative abundance", y=1.0)
        fig.tight_layout()
        fig.savefig(outdir / "taxa_top_boxplots.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    return {
        "rarefied_relative_table": rel,
        "rarefied_presence_absence": presence_absence,
        "unrarefied_relative_table": raw_rel,
        "unrarefied_clr": unrarefied_clr,
        "prevalence_table": prevalence_table,
        "top_taxa": top_taxa,
        "group_summary": group_summary,
    }


def save_taxon_significance_plot(stats: pd.DataFrame, outdir: Path, top_n: int = 20) -> None:
    """Plot the strongest taxon-level FDR results for each model term."""
    if stats is None or stats.empty:
        return
    terms = [
        ("q_treatment", "Treatment"),
        ("q_time_or_phase", "Time/Phase"),
        ("q_treatment_x_time", "Treatment × Time/Phase"),
    ]
    for qcol, label in terms:
        if qcol not in stats.columns or "taxon" not in stats.columns:
            continue
        d = stats[["taxon", qcol]].dropna().copy()
        if d.empty:
            continue
        d[qcol] = pd.to_numeric(d[qcol], errors="coerce")
        d = d.dropna()
        if d.empty:
            continue
        d["score"] = -np.log10(np.clip(d[qcol].to_numpy(float), 1e-300, 1.0))
        d = d.sort_values("score", ascending=False).head(top_n).sort_values("score")
        fig, ax = plt.subplots(figsize=(8, max(5, 0.28 * len(d))))
        ax.barh(np.arange(len(d)), d["score"].to_numpy(float))
        ax.set_yticks(np.arange(len(d)))
        ax.set_yticklabels(d["taxon"].astype(str))
        ax.set_xlabel(r"$-\log_{10}(q)$")
        ax.set_ylabel("Taxon")
        ax.set_title(f"Taxon-level significance — {label}")
        ax.axvline(-math.log10(0.05), linestyle="--", linewidth=1.0, label="FDR 0.05")
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / f"taxa_significance_{_safe_name(label)}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


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
    p.add_argument("--normalization", default=None, help=argparse.SUPPRESS)  # deprecated, ignored
    p.add_argument("--rarefy-depth", default="auto",
                   help="Rarefy level-7 counts before alpha/beta analysis: 'auto' (= smallest non-zero depth, default), "
                        "an integer, or 'none'. Per-taxon CLR statistics always use the unrarefied counts.")
    p.add_argument("--beta-metrics", default="braycurtis,jaccard", help="Comma-separated beta metrics.")
    p.add_argument("--permutations", type=int, default=999, help="Number of PERMANOVA permutations (default: 999).")
    p.add_argument("--seed", type=int, default=42, help="Random seed for rarefaction and permutations.")
    p.add_argument(
        "--taxonomic-level",
        default="species",
        choices=TAXONOMIC_LEVELS + ["all"],
        help="Taxonomic resolution for analysis. Use 'all' to run every available rank (default: species).",
    )
    taxa_group = p.add_mutually_exclusive_group()
    taxa_group.add_argument("--taxa-stats", dest="taxa_stats", action="store_true",
                            help="Run per-taxon CLR mixed-model tests (default; retained for compatibility).")
    taxa_group.add_argument("--no-taxa-stats", dest="taxa_stats", action="store_false",
                            help="Skip per-taxon CLR mixed-model statistics; abundance tables and plots are still produced.")
    p.set_defaults(taxa_stats=True)
    p.add_argument("--top-n-taxa", type=int, default=20,
                   help="Number of top taxa used in abundance plots/tables (default: 20).")
    p.add_argument("--include-unclassified-taxa", action="store_true", help="Include unresolved taxonomy bins in the exploratory per-taxon CLR statistics. They are never included in the headline top-taxa list.")
    p.add_argument("--min-prevalence", type=float, default=0.10, help="Taxon prevalence threshold for --taxa-stats (default: 0.10).")
    p.add_argument("--min-mean-relative", type=float, default=1e-5, help="Minimum mean relative abundance for --taxa-stats.")
    p.add_argument("--no-connect-repeated", action="store_true", help="Do not connect repeated observations from the same Vessel/Pen/Animal on PCoA plots.")
    p.add_argument("--outdir", default="statistics_output", help="Output directory.")
    p.add_argument("--no-exploratory-factorial-permanova", action="store_true", help="Skip the composite-factor exploratory PERMANOVA.")
    args = p.parse_args()

    if args.normalization is not None:
        print("WARNING: --normalization is deprecated and ignored. The toolkit now uses explicit, analysis-specific "
              "representations: rarefied relative abundance, rarefied presence/absence, and unrarefied CLR.", file=sys.stderr)
    if args.permutations < 99:
        raise SystemExit("Use at least 99 permutations; 999 or 9999 is recommended for final results.")
    if args.top_n_taxa < 1:
        raise SystemExit("--top-n-taxa must be at least 1.")
    if not 0 <= args.min_prevalence <= 1:
        raise SystemExit("--min-prevalence must be between 0 and 1.")
    if args.min_mean_relative < 0:
        raise SystemExit("--min-mean-relative must be >= 0.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    table = read_table(Path(args.table))
    meta = read_metadata(Path(args.metadata))
    table = infer_orientation(table, meta.index)
    validate_count_table(table)
    meta.index = meta.index.astype(str).str.strip()

    samples = table.columns.intersection(meta.index)
    if len(samples) == 0:
        raise SystemExit("No samples remain after aligning abundance table and metadata.")
    table = table.loc[:, samples]
    meta = meta.loc[samples].copy()

    design = choose_columns(args, meta)
    meta = metadata_cleanup(meta, design)

    # ---- Sample exclusions (validated AFTER they are applied) ---------------
    excluded: list[tuple[str, str]] = []
    required = [design["treatment_col"], design["time_col"]] + ([design["random_col"]] if design["random_col"] else [])
    complete = meta[required].notna().all(axis=1)
    excluded += [(s, "missing design metadata") for s in meta.index[~complete]]
    meta = meta.loc[complete]
    table = table.loc[:, meta.index]
    if table.empty:
        raise SystemExit("No samples remain after design metadata filtering.")

    raw_depth = table.sum(axis=0)
    zero = raw_depth <= 0
    excluded += [(s, "zero total abundance") for s in table.columns[zero]]
    table = table.loc[:, ~zero]
    meta = meta.loc[table.columns]
    raw_depth = raw_depth[~zero]

    rar_arg = str(args.rarefy_depth).lower()
    rarefied_to = None
    if rar_arg != "none":
        rarefied_to = int(raw_depth.min()) if rar_arg == "auto" else int(rar_arg)
        low = raw_depth < rarefied_to
        excluded += [(s, f"depth {int(raw_depth[s])} < rarefaction depth {rarefied_to}") for s in table.columns[low]]
        table = table.loc[:, ~low]
        meta = meta.loc[table.columns]
        raw_depth = raw_depth[~low]
    if table.empty:
        raise SystemExit("No samples remain after depth filtering.")

    if excluded:
        pd.DataFrame(excluded, columns=["sample_id", "reason"]).to_csv(outdir / "excluded_samples.tsv", sep="\t", index=False)
        for s, why in excluded:
            print(f"EXCLUDED: {s}: {why}", file=sys.stderr)

    validation = validate_design(meta, design, excluded)
    validation.to_csv(outdir / "design_validation.tsv", sep="\t", index=False)

    # Explicit repeated-measures audit: useful when a sample is excluded for zero
    # abundance and leaves one vessel/pen with only a single retained observation.
    if design["random_col"]:
        rnd, tm, tr = design["random_col"], design["time_col"], design["treatment_col"]
        pairing_rows = []
        all_blocks = pd.Index(meta[rnd].dropna().astype(str).unique())
        for block in sorted(all_blocks):
            sub = meta[meta[rnd].astype(str) == block]
            times = list(pd.unique(sub[tm].astype(str)))
            pairing_rows.append({
                "experimental_unit": block,
                "treatment": "; ".join(pd.unique(sub[tr].astype(str))),
                "n_retained_samples": len(sub),
                "retained_time_or_phase": "; ".join(times),
                "paired_complete": bool(len(sub) == 2 and len(times) == 2),
                "status": "COMPLETE" if len(sub) == 2 and len(times) == 2 else "INCOMPLETE_AFTER_EXCLUSION",
            })
        if pairing_rows:
            pd.DataFrame(pairing_rows).to_csv(outdir / "repeated_measures_audit.tsv", sep="\t", index=False)

    for _, row in validation.iterrows():
        if row["status"] == "CHECK":
            print(f"DESIGN CHECK: {row['check']}: {row['detail']}", file=sys.stderr)

    table_raw = table
    if rarefied_to is not None:
        table_ab = rarefy(table_raw, rarefied_to, args.seed)
    else:
        table_ab = table_raw
    pd.DataFrame({"raw_depth": raw_depth, "analysis_depth": table_ab.sum(axis=0)}).to_csv(
        outdir / "sequencing_depth.tsv", sep="\t")

    available_levels = taxonomy_levels_available(table_raw)
    requested_levels = available_levels if args.taxonomic_level == "all" else [args.taxonomic_level]
    if args.taxonomic_level not in available_levels and args.taxonomic_level != "all":
        raise SystemExit(
            f"Requested taxonomic level '{args.taxonomic_level}' could not be recovered from the input taxonomy strings. "
            f"Available levels detected: {', '.join(available_levels) or 'none'}."
        )
    if not requested_levels:
        raise SystemExit("No taxonomic ranks could be recovered from the input table.")

    tr_col, tm_col, rnd_col = design["treatment_col"], design["time_col"], design["random_col"]
    run_summary = []

    for level in requested_levels:
        level_outdir = outdir / level if len(requested_levels) > 1 else outdir
        level_outdir.mkdir(parents=True, exist_ok=True)
        rank_map = taxonomy_rank_mapping(table_raw, level)
        level_raw = aggregate_taxonomic_level(table_raw, level)
        level_table = aggregate_taxonomic_level(table_ab, level)   # alpha / beta input
        status_by_taxon = (rank_map.groupby("aggregated_taxon")["taxon_status"]
                           .agg(lambda x: "unresolved" if (x == "unresolved").any() else "resolved"))
        rank_map.to_csv(level_outdir / "taxonomy_rank_mapping.tsv", sep="\t", index=False)
        pd.DataFrame({"taxon": status_by_taxon.index, "taxon_status": status_by_taxon.values}).to_csv(
            level_outdir / "taxonomy_rank_status.tsv", sep="\t", index=False
        )

        params = vars(args).copy()
        params["resolved_design"] = design
        params["taxonomic_level"] = level
        params["available_taxonomic_levels"] = available_levels
        params["rarefied_to"] = rarefied_to
        params["n_samples"] = int(level_table.shape[1])
        params["n_taxa"] = int(level_table.shape[0])
        params["normalisation_outputs"] = {
            "rarefied_relative_abundance": "rarefied_relative_abundance.tsv",
            "rarefied_presence_absence": "rarefied_presence_absence.tsv",
            "unrarefied_relative_abundance": "unrarefied_relative_abundance.tsv",
            "unrarefied_clr": "unrarefied_clr.tsv",
            "normalisation_manifest": "normalisation_manifest.tsv",
        }
        with open(level_outdir / "analysis_parameters.json", "w") as fh:
            json.dump(params, fh, indent=2, default=str)

        pd.DataFrame([{
            "study": design["study"], "taxonomic_level": level,
            "n_samples": level_table.shape[1], "n_taxa": level_table.shape[0],
            "treatments": "; ".join(sorted(meta[tr_col].dropna().astype(str).unique())),
            "time_or_phase": "; ".join(sorted(meta[tm_col].dropna().astype(str).unique())),
            "random_effect": rnd_col or "NONE (OLS/fixed-effect fallback)",
            "rarefied_to": rarefied_to if rarefied_to is not None else "none",
            "unresolved_taxa_handling": "retained in community composition; excluded from headline top-taxa and per-taxon CLR by default",
        }]).to_csv(level_outdir / "design_summary.tsv", sep="\t", index=False)
        meta.to_csv(level_outdir / "analysis_metadata.tsv", sep="\t")

        # Explicit community-level representations.
        rel = normalise(level_table, "relative")
        rel.to_csv(level_outdir / "rarefied_relative_abundance.tsv", sep="\t")
        presence_absence = (level_table > 0).astype(np.int8)
        presence_absence.to_csv(level_outdir / "rarefied_presence_absence.tsv", sep="\t")
        normalise(level_raw, "relative").to_csv(level_outdir / "unrarefied_relative_abundance.tsv", sep="\t")
        level_raw.to_csv(level_outdir / "unrarefied_counts.tsv", sep="\t")
        clr_samples(level_raw).to_csv(level_outdir / "unrarefied_clr.tsv", sep="\t")

        unresolved_labels = set(status_by_taxon.index[status_by_taxon == "unresolved"])
        alpha = alpha_metrics(level_table, unresolved_labels=unresolved_labels)
        alpha["raw_depth"] = raw_depth.reindex(alpha.index)
        alpha = alpha.join(meta, how="inner")
        alpha.to_csv(level_outdir / "alpha_diversity.tsv", sep="\t")
        save_alpha_plots(alpha[["shannon", "simpson", "observed_features", "observed_resolved_taxa", "read_count_or_total"]],
                         meta, tr_col, tm_col, level_outdir)

        metrics = ["shannon", "simpson", "observed_features", "observed_resolved_taxa"]
        omnibus = pd.concat([alpha_omnibus(alpha, m, tr_col, tm_col, rnd_col) for m in metrics], ignore_index=True)
        if not omnibus.empty:
            # BH across the three alpha metrics, separately for each model term.
            omnibus = _bh_within(omnibus, "term")
            omnibus.to_csv(level_outdir / "alpha_omnibus_statistics.tsv", sep="\t", index=False)
            fb = omnibus.loc[omnibus["note"].astype(str) != "", ["metric", "note"]].drop_duplicates()
            for _, r in fb.iterrows():
                print(f"MODEL NOTE ({level}/{r['metric']}): {r['note']}", file=sys.stderr)

        pairwise = pd.concat([alpha_pairwise(alpha, m, tr_col, tm_col, rnd_col) for m in metrics], ignore_index=True)
        if not pairwise.empty:
            pairwise.to_csv(level_outdir / "alpha_pairwise_contrasts.tsv", sep="\t", index=False)

        for metric in [m.strip() for m in args.beta_metrics.split(",") if m.strip()]:
            if metric.lower() == "jaccard":
                dist = beta_distance(presence_absence, metric)
            else:
                dist = beta_distance(rel, metric)
            dist.to_csv(level_outdir / f"beta_distance_{metric}.tsv", sep="\t")

            try:
                coords, variance_explained = pcoa(dist)
                coords = coords.join(meta[[tr_col, tm_col]], how="left")
                if rnd_col:
                    coords = coords.join(meta[[rnd_col]], how="left")
                coords.to_csv(level_outdir / f"beta_pcoa_coordinates_{metric}.tsv", sep="\t")
                variance_explained.to_csv(level_outdir / f"beta_pcoa_variance_explained_{metric}.tsv", sep="\t", header=True)
                if {"PC1", "PC2"}.issubset(coords.columns):
                    save_pcoa_plot(coords, meta, tr_col, tm_col, level_outdir, metric=metric,
                                   variance_explained=variance_explained,
                                   random_col=None if args.no_connect_repeated else rnd_col,
                                   taxonomic_level=level)
            except ValueError as exc:
                pd.DataFrame([{"metric": metric, "taxonomic_level": level, "status": "SKIPPED", "reason": str(exc)}]).to_csv(
                    level_outdir / f"beta_pcoa_status_{metric}.tsv", sep="\t", index=False)
                print(f"PCoA skipped for {level}/{metric}: {exc}", file=sys.stderr)

            permanova_time_stratified(dist, meta, tr_col, tm_col, args.permutations, args.seed).to_csv(
                level_outdir / f"beta_permanova_treatment_within_phase_{metric}.tsv", sep="\t", index=False)
            paired_time_permanova(dist, meta, tr_col, tm_col, rnd_col, args.permutations, args.seed, level_outdir, metric)
            paired_change_magnitude(dist, meta, tr_col, tm_col, rnd_col, args.permutations, args.seed, level_outdir, metric)
            if not args.no_exploratory_factorial_permanova:
                exploratory_full_factorial_permanova(dist, meta, tr_col, tm_col, args.permutations, args.seed).to_csv(
                    level_outdir / f"beta_permanova_exploratory_composite_{metric}.tsv", sep="\t", index=False)

        # Taxon-level abundance products are always created: a one-stop analysis
        # should leave both the inferential table and the normalised/compositional
        # summaries/figures needed to interpret it.
        taxon_outputs = taxon_abundance_outputs(
            level_raw, level_table, meta, design, level_outdir,
            top_n=args.top_n_taxa, taxonomy_status=status_by_taxon
        )

        taxon_stats_out = pd.DataFrame()
        if args.taxa_stats:
            # CLR is compositional: use the UNRAREFIED counts for the inferential
            # per-taxon model. The exact matrix is written as unrarefied_clr.tsv.
            # Relative-abundance plots use rarefied relative abundance.
            taxon_stats_out = taxon_statistics(
                level_raw, meta, design, level_outdir, args.min_prevalence,
                args.min_mean_relative, args.include_unclassified_taxa
            )
            save_taxon_significance_plot(taxon_stats_out, level_outdir, top_n=args.top_n_taxa)

        # A compact machine-readable summary for downstream reporting.
        summary_row = {
            "taxonomic_level": level,
            "n_samples": level_table.shape[1],
            "n_taxa": level_table.shape[0],
            "n_resolved_taxa": int((status_by_taxon == "resolved").sum()),
            "n_unresolved_taxa_bins": int((status_by_taxon == "unresolved").sum()),
            "n_top_taxa_plotted": len(taxon_outputs["top_taxa"]),
            "taxa_statistics_run": bool(args.taxa_stats),
            "n_taxa_statistics": int(len(taxon_stats_out)),
            "output_directory": str(level_outdir.resolve()),
        }
        run_summary.append(summary_row)

    pd.DataFrame(run_summary).to_csv(outdir / "taxonomic_level_run_summary.tsv", sep="\t", index=False)
    print("Analysis complete.")
    print(f"Study: {design['study']}")
    print(f"Samples analysed: {table.shape[1]}" + (f" (rarefied to {rarefied_to} reads)" if rarefied_to else ""))
    print(f"Treatment column: {tr_col}")
    print(f"Time/phase column: {tm_col}")
    print(f"Random/repeated-measures column: {rnd_col or 'NONE'}")
    print(f"Taxonomic levels analysed: {', '.join(requested_levels)}")
    print(f"Taxon abundance outputs: {'enabled' if args.top_n_taxa > 0 else 'disabled'}")
    print("Normalisation outputs: rarefied_relative_abundance, rarefied_presence_absence, unrarefied_relative_abundance, unrarefied_clr")
    print(f"Per-taxon CLR statistics: {'enabled' if args.taxa_stats else 'disabled'}")
    print(f"Results: {outdir.resolve()}")


if __name__ == "__main__":
    main()
