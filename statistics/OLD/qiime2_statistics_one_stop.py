#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from scipy.spatial.distance import pdist, squareform


def read_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=None, engine="python")
    df = df.dropna(how="all")
    first_col = df.columns[0]
    df = df.set_index(first_col)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    return df


def orient_table(table: pd.DataFrame, sample_ids: pd.Index) -> pd.DataFrame:
    row_hits = table.index.astype(str).intersection(sample_ids.astype(str))
    col_hits = table.columns.astype(str).intersection(sample_ids.astype(str))
    if len(col_hits) >= len(row_hits) and len(col_hits) > 0:
        oriented = table.copy()
        oriented.columns = oriented.columns.astype(str).str.strip()
        return oriented.loc[:, oriented.columns.intersection(sample_ids.astype(str))]
    if len(row_hits) > 0:
        oriented = table.copy()
        oriented.index = oriented.index.astype(str).str.strip()
        return oriented.loc[oriented.index.intersection(sample_ids.astype(str))]
    raise SystemExit("No overlapping samples were found between the table and metadata.")


def read_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path, sep=None, engine="python")
    meta = meta.dropna(how="all")
    sample_col = meta.columns[0]
    meta = meta.rename(columns={sample_col: "sample_id"})
    return meta.set_index("sample_id")


def normalize(table: pd.DataFrame, method: str) -> pd.DataFrame:
    method = method.lower()
    if method == "none":
        return table.copy()
    if method == "relative":
        return table.div(table.sum(axis=0), axis=1).fillna(0)
    if method == "clr":
        pseudo = 1.0
        x = table + pseudo
        gm = np.exp(np.log(x).mean(axis=0))
        return np.log(x.div(gm, axis=1))
    if method == "zscore":
        return table.sub(table.mean(axis=1), axis=0).div(table.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    raise ValueError(f"Unknown normalization method: {method}")


def shannon(col: pd.Series) -> float:
    p = col[col > 0]
    p = p / p.sum()
    return float(-(p * np.log(p)).sum())


def simpson(col: pd.Series) -> float:
    p = col[col > 0]
    p = p / p.sum()
    return float(1 - np.sum(p ** 2))


def alpha_metrics(table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    samples = table.columns if table.shape[1] <= table.shape[0] else table.index
    for sample in samples:
        col = table[sample] if sample in table.columns else table.loc[sample]
        rows.append(
            {
                "sample_id": sample,
                "observed_features": int((col > 0).sum()),
                "shannon": shannon(col),
                "simpson": simpson(col),
                "richness": float(col.sum()),
            }
        )
    return pd.DataFrame.from_records(rows).set_index("sample_id")


def beta_distance(table: pd.DataFrame, metric: str) -> pd.DataFrame:
    metric = metric.lower()
    x = table.T.values if table.shape[1] <= table.shape[0] else table.values
    labels = table.columns if table.shape[1] <= table.shape[0] else table.index
    if metric == "braycurtis":
        dist = pdist(x, metric="braycurtis")
    elif metric == "jaccard":
        dist = pdist(x > 0, metric="jaccard")
    elif metric == "euclidean":
        dist = pdist(x, metric="euclidean")
    else:
        raise ValueError(f"Unknown beta metric: {metric}")
    return pd.DataFrame(squareform(dist), index=labels, columns=labels)


def save_alpha_plots(alpha: pd.DataFrame, meta: pd.DataFrame, group_col: str, outdir: Path) -> None:
    merged = alpha.join(meta[[group_col]], how="inner")
    for metric in ["shannon", "simpson", "observed_features"]:
        plt.figure(figsize=(7, 5))
        merged.boxplot(column=metric, by=group_col)
        plt.title(f"{metric.title()} by {group_col}")
        plt.suptitle("")
        plt.tight_layout()
        plt.savefig(outdir / f"alpha_{metric}_by_{group_col}.png", dpi=300)
        plt.close()


def save_alpha_timecourse(alpha: pd.DataFrame, meta: pd.DataFrame, group_col: str, time_col: str, outdir: Path) -> None:
    # time_col can be numeric or categorical representing ordered timepoints
    merged = alpha.join(meta[[group_col, time_col]], how="inner").reset_index()
    # try to coerce time to numeric if possible for ordering
    try:
        merged[time_col] = pd.to_numeric(merged[time_col])
        sort_by_time = True
    except Exception:
        # leave as categorical
        merged[time_col] = merged[time_col].astype(str)
        sort_by_time = False
    groups = merged[group_col].unique()
    plt.figure(figsize=(8, 6))
    for g in groups:
        sub = merged[merged[group_col] == g]
        if sort_by_time:
            agg = sub.groupby(time_col)["shannon"].agg(["mean", "sem"]).reset_index()
            x = agg[time_col]
            y = agg["mean"]
            yerr = agg["sem"]
        else:
            # treat time as categorical ordered as encountered
            agg = sub.groupby(time_col)["shannon"].agg(["mean", "sem"]).reindex().reset_index()
            x = agg[time_col]
            y = agg["mean"]
            yerr = agg["sem"]
        plt.errorbar(x, y, yerr=yerr, label=str(g), marker='o')
    plt.xlabel(time_col)
    plt.ylabel("Shannon index")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / f"alpha_timecourse_shannon_by_{group_col}_vs_{time_col}.png", dpi=300)
    plt.close()


def save_pcoa(dist: pd.DataFrame, meta: pd.DataFrame, group_col: str, outdir: Path) -> None:
    coords = MDS(n_components=2, dissimilarity="precomputed", random_state=42).fit_transform(dist.values)
    pcoa = pd.DataFrame(coords, index=dist.index, columns=["PC1", "PC2"]).join(meta[[group_col]], how="left")
    plt.figure(figsize=(7, 6))
    for group, sub in pcoa.groupby(group_col):
        plt.scatter(sub["PC1"], sub["PC2"], label=str(group), s=60)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / f"beta_pcoa_{group_col}.png", dpi=300)
    plt.close()
    pcoa.to_csv(outdir / "beta_pcoa_coordinates.tsv", sep="\t")


def save_pca(table: pd.DataFrame, meta: pd.DataFrame, group_col: str, outdir: Path) -> None:
    x = table.T.values if table.shape[1] <= table.shape[0] else table.values
    labels = table.columns if table.shape[1] <= table.shape[0] else table.index
    coords = PCA(n_components=2, random_state=42).fit_transform(x)
    pca = pd.DataFrame(coords, index=labels, columns=["PC1", "PC2"]).join(meta[[group_col]], how="left")
    plt.figure(figsize=(7, 6))
    for group, sub in pca.groupby(group_col):
        plt.scatter(sub["PC1"], sub["PC2"], label=str(group), s=60)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / f"pca_{group_col}.png", dpi=300)
    plt.close()
    pca.to_csv(outdir / "pca_coordinates.tsv", sep="\t")


def main() -> None:
    p = argparse.ArgumentParser(description="One-stop QIIME2 statistics script")
    p.add_argument("--table", required=True, help="QIIME2 level-7 CSV/TSV feature table")
    p.add_argument("--metadata", required=True, help="Sample metadata TSV/CSV with sample IDs in first column")
    p.add_argument("--group-col", help="Metadata column that labels control/trial groups (kept for compatibility)")
    p.add_argument("--group-cols", help="Comma-separated metadata columns to combine as a single grouping (e.g. Treatment,Replicate)", default=None)
    p.add_argument("--time-col", help="Optional metadata column for timepoints (will produce time-course alpha plots)", default=None)
    p.add_argument("--normalization", default="relative", choices=["none", "relative", "clr", "zscore"])
    p.add_argument("--beta-metric", default="braycurtis", choices=["braycurtis", "jaccard", "euclidean"])
    p.add_argument("--outdir", default="statistics_output", help="Output directory")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    table = read_table(Path(args.table))
    meta = read_metadata(Path(args.metadata))
    table.columns = table.columns.astype(str).str.strip()
    meta.index = meta.index.astype(str).str.strip()
    table = orient_table(table, meta.index)
    if table.index.dtype == object and len(table.columns) and table.columns.isin(meta.index).any():
        table = table.loc[:, table.columns.intersection(meta.index)]
        meta = meta.loc[table.columns]
    else:
        meta = meta.loc[table.index]
    # determine grouping column(s)
    if args.group_cols:
        cols = [c.strip() for c in args.group_cols.split(",") if c.strip()]
        missing = [c for c in cols if c not in meta.columns]
        if missing:
            raise SystemExit(f"Missing metadata columns: {', '.join(missing)}")
        combined_name = "group_combined"
        meta[combined_name] = meta[cols].astype(str).agg("|".join, axis=1)
        group_col_used = combined_name
    elif args.group_col:
        if args.group_col not in meta.columns:
            raise SystemExit(f"Missing metadata column: {args.group_col}")
        group_col_used = args.group_col
    else:
        # try to auto-detect common treatment/group column names
        candidates = [c for c in ["Treatment", "treatment", "TreatmentID", "treat", "Group", "group", "Condition", "condition"] if c in meta.columns]
        if candidates:
            group_col_used = candidates[0]
            print(f"Auto-detected grouping column: {group_col_used}")
        else:
            raise SystemExit("No grouping column provided or detected. Use --group-col or --group-cols.")

    # detect or validate time column for longitudinal plots
    if args.time_col:
        if args.time_col not in meta.columns:
            raise SystemExit(f"Missing time column: {args.time_col}")
        time_col = args.time_col
    else:
        time_candidates = [c for c in ["Time", "time", "Timepoint", "timepoint", "Day", "day", "TP"] if c in meta.columns]
        time_col = time_candidates[0] if time_candidates else None
        if time_col:
            print(f"Auto-detected time column: {time_col}")

    if table.empty:
        raise SystemExit("No overlapping samples were found between the table and metadata.")

    norm = normalize(table, args.normalization)
    norm.to_csv(outdir / f"normalized_table_{args.normalization}.tsv", sep="\t")

    alpha = alpha_metrics(norm)
    alpha_joined = alpha.join(meta[[group_col_used]])
    alpha_joined.to_csv(outdir / "alpha_diversity.tsv", sep="\t")
    save_alpha_plots(alpha, meta, group_col_used, outdir)
    if time_col:
        save_alpha_timecourse(alpha, meta, group_col_used, time_col, outdir)

    dist = beta_distance(norm, args.beta_metric)
    dist.to_csv(outdir / f"beta_distance_{args.beta_metric}.tsv", sep="\t")
    save_pcoa(dist, meta, group_col_used, outdir)
    save_pca(norm, meta, group_col_used, outdir)


if __name__ == "__main__":
    main()
