"""End-to-end orchestration."""

from __future__ import annotations

import json
import platform
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__, alpha as alpha_mod, beta as beta_mod, first_pass as first_pass_mod, plots, taxa as taxa_mod
from .common import Notes, safe_name, sha256
from .design import DesignSpec, FixedDesign, analyse_rank, build_frame, validate_design
from .io import InputError, read_feature_table, read_metadata, validate_counts
from .rarefaction import choose_depth, rarefaction_curves, rarefy_matrix
from .taxonomy import LEVELS, RANKS, available_levels, build_level_map, lineage_frame


def _w(df: pd.DataFrame, path: Path, index: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=index)


def _split(v) -> list:
    if v is None:
        return []
    if isinstance(v, (list, tuple)):
        return [str(x).strip() for x in v if str(x).strip()]
    return [x.strip() for x in str(v).split(",") if x.strip()]


def resolve_levels(arg: str, lineages, notes) -> list[str]:
    avail = available_levels(lineages) if lineages is not None else []
    a = (arg or "auto").strip().lower()
    if a == "auto":
        if "genus" in avail:
            return ["genus"]
        if avail:
            return [avail[-1]]
        return ["feature"]
    if a == "all":
        return list(dict.fromkeys(["feature"] + avail))
    out = []
    for lv in _split(a):
        if lv not in LEVELS:
            raise InputError(f"Unknown taxonomic level '{lv}'. Choose from {LEVELS}, 'all' or 'auto'.")
        if lv != "feature" and lv not in avail:
            notes.warn("levels", f"Level '{lv}' not available in the taxonomy (available: {avail}); skipped.")
            continue
        out.append(lv)
    if not out:
        raise InputError("No requested taxonomic level is available.")
    return out


def _primary_report_level(request: str, levels: list[str]) -> str:
    requested = str(request).strip().lower()
    if requested != "auto":
        if requested not in levels:
            raise InputError(
                f"Primary report level '{requested}' must be among requested levels: {levels}."
            )
        return requested
    if "genus" in levels:
        return "genus"
    taxonomic_levels = [level for level in levels if level != "feature"]
    if taxonomic_levels:
        return max(taxonomic_levels, key=LEVELS.index)
    return "feature"


def _alpha_sensitivity_comparison(primary: pd.DataFrame, sensitivity: pd.DataFrame) -> pd.DataFrame:
    columns = ["metric", "term", "p_value", "p_adj_across_metrics"]
    left = primary[[c for c in columns if c in primary]].rename(
        columns={"p_value": "unrarefied_p", "p_adj_across_metrics": "unrarefied_adjusted_p"})
    right = sensitivity[[c for c in columns + ["p_min_over_draws", "p_max_over_draws"] if c in sensitivity]].rename(
        columns={"p_value": "rarefied_p", "p_adj_across_metrics": "rarefied_adjusted_p",
                 "p_min_over_draws": "rarefied_p_min_across_draws",
                 "p_max_over_draws": "rarefied_p_max_across_draws"})
    if left.empty or right.empty:
        return pd.DataFrame()
    comparison = left.merge(right, on=["metric", "term"], how="inner")
    base_sig = comparison["unrarefied_adjusted_p"] < 0.05
    rare_sig = comparison["rarefied_adjusted_p"] < 0.05
    missing = comparison["unrarefied_adjusted_p"].isna() | comparison["rarefied_adjusted_p"].isna()
    comparison["threshold_classification"] = np.select(
        [missing, base_sig == rare_sig],
        ["not estimable in one analysis", "same at 0.05 threshold"],
        default="differs at 0.05 threshold",
    )
    return comparison


def _beta_sensitivity_comparison(primary: pd.DataFrame, sensitivity: pd.DataFrame) -> pd.DataFrame:
    if (primary.empty or sensitivity.empty or "term" not in primary or "term" not in sensitivity
            or "p_value" not in primary or "p_value" not in sensitivity):
        return pd.DataFrame()
    optional = [c for c in ("min_attainable_p", "exact") if c in primary and c in sensitivity]
    left = primary[["term", "p_value", *optional]].rename(
        columns={"p_value": "unrarefied_p", "min_attainable_p": "unrarefied_min_attainable_p",
                 "exact": "unrarefied_exact"})
    right = sensitivity[["term", "p_value", *optional]].rename(
        columns={"p_value": "rarefied_p", "min_attainable_p": "rarefied_min_attainable_p",
                 "exact": "rarefied_exact"})
    comparison = left.merge(right, on="term", how="inner")
    base_sig = comparison["unrarefied_p"] < 0.05
    rare_sig = comparison["rarefied_p"] < 0.05
    missing = comparison["unrarefied_p"].isna() | comparison["rarefied_p"].isna()
    comparison["threshold_classification"] = np.select(
        [missing, base_sig == rare_sig],
        ["not estimable in one analysis", "same at 0.05 threshold"],
        default="differs at 0.05 threshold",
    )
    return comparison


def _rarefied_taxon_prevalence(aggregated: np.ndarray, ids, level_map, frame, draw: int) -> pd.DataFrame:
    counts = pd.DataFrame(aggregated, index=ids, columns=level_map.labels)
    detected = counts.gt(0)
    meta = frame.df.loc[ids]
    unit = meta["upper"].where(meta["upper"].notna(), pd.Series(ids, index=ids)).astype(str)
    rows = []
    for time in (frame.time_levels if frame.has_time else [None]):
        for group in frame.group_levels:
            mask = meta["group"].eq(group)
            if time is not None:
                mask &= meta["time"].eq(time)
            sample_ids = meta.index[mask]
            unit_detection = detected.loc[sample_ids].groupby(unit.loc[sample_ids], sort=False).max()
            if unit_detection.empty:
                continue
            stratum = str(group) if time is None else f"{group} / {time}"
            for taxon, prevalence in unit_detection.mean(axis=0).items():
                rows.append({
                    "draw": draw,
                    "taxon": taxon,
                    "group": group,
                    "time": time if time is not None else "",
                    "stratum": stratum,
                    "mean_prevalence": float(prevalence),
                    "n_units": len(unit_detection),
                })
    return pd.DataFrame(rows)


def _run_rarefaction_sensitivity(counts: pd.DataFrame, level_map, level: str, metrics_a: list[str],
                                 metrics_b: list[str], frame, fixed, rank, args, notes, out: Path,
                                 primary_alpha: dict, primary_beta: dict) -> dict:
    ids = list(frame.df.index)
    depth = int(counts[ids].sum(0).min())
    iterations = int(args.rarefaction_sensitivity_iterations)
    rng = np.random.default_rng(np.random.SeedSequence([int(args.seed), 8341]))
    C = counts[ids].T.to_numpy()
    raw = pd.DataFrame(np.asarray(level_map.matrix @ C.T).T, index=ids, columns=level_map.labels)
    sens_out = out / "sensitivity" / "rarefaction"
    _w(pd.DataFrame([{"method": "without-replacement multinomial hypergeometric draws",
                      "depth": depth, "iterations": iterations, "samples": len(ids),
                      "taxonomic_level": level}]), sens_out / "settings.tsv")
    notes.info("rarefaction_sensitivity",
               f"Compared the unrarefied primary analysis with {iterations} without-replacement draws at "
               f"{depth} reads per sample; primary results were not replaced.")
    curves = rarefaction_curves(C, ids, rng, max_depth=depth)
    _w(curves, out / "qc" / "rarefaction_sensitivity_curves.tsv")
    plots.rarefaction_plot(curves, depth, out / "qc" / "rarefaction_sensitivity_curves.png")
    if depth < 0.1 * float(counts[ids].sum(0).median()):
        notes.warn("rarefaction_sensitivity",
                   f"The minimum-library sensitivity depth ({depth}) is below 10% of the median library size; "
                   "this comparison may be strongly affected by the shallowest sample.")

    alpha_comparison = pd.DataFrame()
    if not args.skip_alpha:
        raw_alpha = alpha_mod.compute_alpha(raw.to_numpy(), ids, metrics_a, level == "feature")
        alpha_draws = []
        for _ in range(iterations):
            rarefied = rarefy_matrix(C, depth, rng)
            aggregated = np.asarray(level_map.matrix @ rarefied.T).T
            draw = alpha_mod.compute_alpha(aggregated, ids, metrics_a, level == "feature")
            for metric in ("chao1", "goods_coverage"):
                if metric in draw and metric in raw_alpha:
                    draw[metric] = raw_alpha[metric]
            alpha_draws.append(draw)
        alpha_result = alpha_mod.run_alpha(
            alpha_draws, frame, fixed, rank, args.alpha_transform, args.p_adjust, notes,
            f"{level}/rarefaction_sensitivity/alpha",
        )
        alpha_comparison = _alpha_sensitivity_comparison(primary_alpha["omnibus"], alpha_result["omnibus"])
        _w(alpha_result["omnibus"], sens_out / "rarefied_alpha_omnibus.tsv")
        _w(alpha_comparison, sens_out / "alpha_comparison.tsv")
        alpha_fig_dir = sens_out / "figures"
        for metric in raw_alpha.columns:
            plots.rarefaction_alpha_comparison(
                raw_alpha, alpha_draws, frame, metric, alpha_fig_dir / f"alpha_{safe_name(metric)}.png"
            )

    beta_comparisons = {}
    beta_metrics = [m for m in metrics_b if m in ("braycurtis", "jaccard")]
    prevalence_draws = []
    if not args.skip_beta:
        for metric in beta_metrics:
            distance_sum = np.zeros((len(ids), len(ids)), dtype=float)
            for _ in range(iterations):
                rarefied = rarefy_matrix(C, depth, rng)
                aggregated = np.asarray(level_map.matrix @ rarefied.T).T
                table = pd.DataFrame(aggregated.T, columns=ids)
                source = (table.div(table.sum(0).replace(0, np.nan), axis=1).fillna(0)
                          if metric == "braycurtis" else table)
                distance_sum += beta_mod.distance_from_table(source, metric).to_numpy()
            mean_distance = pd.DataFrame(distance_sum / iterations, index=ids, columns=ids)
            _w(mean_distance, sens_out / f"rarefied_{safe_name(metric)}_mean_distance.tsv", index=True)
            try:
                coordinates, variance, _ = beta_mod.pcoa(mean_distance)
                _w(coordinates, sens_out / f"rarefied_{safe_name(metric)}_pcoa_coordinates.tsv", index=True)
                plots.pcoa_plot(
                    coordinates, variance, frame,
                    sens_out / "figures" / f"pcoa_{safe_name(metric)}.png",
                    f"{level}: rarefied {metric} PCoA", not args.no_connect_repeated,
                )
            except ValueError as e:
                notes.warn(f"{level}/rarefaction_sensitivity/beta/{metric}", f"PCoA skipped: {e}")
            beta_result = beta_mod.run_beta(
                mean_distance, frame, fixed, int(args.permutations), rng, args.p_adjust, notes,
                f"{level}/rarefaction_sensitivity/beta/{metric}", args.change_magnitude,
            )
            comparison = _beta_sensitivity_comparison(
                primary_beta.get(metric, {}).get("terms", pd.DataFrame()),
                beta_result.get("terms", pd.DataFrame()),
            )
            _w(beta_result.get("terms", pd.DataFrame()), sens_out / f"rarefied_beta_{safe_name(metric)}_terms.tsv")
            _w(comparison, sens_out / f"beta_{safe_name(metric)}_comparison.tsv")
            beta_comparisons[metric] = comparison

    for draw_index in range(iterations):
        rarefied = rarefy_matrix(C, depth, rng)
        aggregated = np.asarray(level_map.matrix @ rarefied.T).T
        prevalence_draws.append(_rarefied_taxon_prevalence(aggregated, ids, level_map, frame, draw_index + 1))
    prevalence_by_draw = pd.concat(prevalence_draws, ignore_index=True) if prevalence_draws else pd.DataFrame()
    if len(prevalence_by_draw):
        prevalence_summary = (prevalence_by_draw.groupby(["taxon", "group", "time", "stratum"], sort=False)
                              .agg(mean_prevalence=("mean_prevalence", "mean"),
                                   lower_5pct=("mean_prevalence", lambda x: float(x.quantile(0.05))),
                                   upper_95pct=("mean_prevalence", lambda x: float(x.quantile(0.95))),
                                   n_units=("n_units", "first"))
                              .reset_index())
    else:
        prevalence_summary = pd.DataFrame()
    if len(prevalence_summary):
        raw_prevalence = _rarefied_taxon_prevalence(raw.to_numpy(), ids, level_map, frame, 0)
        raw_prevalence = raw_prevalence[["taxon", "group", "time", "mean_prevalence"]].rename(
            columns={"mean_prevalence": "unrarefied_prevalence"}
        )
        prevalence_summary = prevalence_summary.merge(
            raw_prevalence, on=["taxon", "group", "time"], how="left"
        )
    _w(prevalence_by_draw, sens_out / "rarefied_taxon_prevalence_by_draw.tsv")
    _w(prevalence_summary, sens_out / "rarefied_taxon_prevalence_summary.tsv")
    if len(prevalence_summary):
        plots.rarefaction_taxa_prevalence(
            prevalence_summary, sens_out / "figures" / "taxon_detection_prevalence.png"
        )

    return dict(depth=depth, iterations=iterations, level=level, alpha=alpha_comparison,
                beta=beta_comparisons, metrics=beta_metrics, alpha_metrics=list(raw_alpha.columns) if not args.skip_alpha else [],
                beta_figures=beta_metrics, taxon_prevalence=prevalence_summary)


def run(args) -> int:
    notes = Notes()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    plots.set_dpi(args.dpi)
    nperm = int(args.permutations)

    # ---------------- load and align -----------------
    meta = read_metadata(Path(args.metadata), args.sample_id_col)
    counts, lineages = read_feature_table(Path(args.table), meta.index, Path(args.taxonomy) if args.taxonomy else None, notes)
    validate_counts(counts)
    counts = counts.astype(np.int64)
    excluded: list[tuple] = []
    for s in counts.columns.difference(meta.index):
        excluded.append((s, "table sample missing from metadata", "alignment"))
    for s in meta.index.difference(counts.columns):
        excluded.append((s, "metadata sample missing from table", "alignment"))
    shared = [s for s in meta.index if s in set(counts.columns)]
    drop = set(_split(args.exclude_samples))
    if args.exclude_regex:
        rx = re.compile(args.exclude_regex)
        drop |= {s for s in shared if rx.search(s)}
    for s in shared:
        if s in drop:
            excluded.append((s, "excluded by user", "user"))
    shared = [s for s in shared if s not in drop]
    meta = meta.loc[shared]
    counts = counts[shared]

    if args.exclude_taxa_regex:
        rx = re.compile(args.exclude_taxa_regex, re.I)
        text = lineages.reindex(counts.index).fillna("").astype(str) if lineages is not None else counts.index.astype(str)
        m = pd.Series([bool(rx.search(a) or rx.search(str(b))) for a, b in zip(text, counts.index)], index=counts.index)
        notes.info("filter", f"--exclude-taxa-regex removed {int(m.sum())} features (e.g. host/chloroplast/mitochondria).")
        counts = counts.loc[~m]

    # ---------------- design -----------------
    spec = DesignSpec(group_col=args.group_col, time_col=args.time_col, subject_col=args.subject_col,
                      cluster_col=args.cluster_col, covariates=tuple(_split(args.covariates)),
                      group_order=_split(args.group_order) or None, time_order=_split(args.time_order) or None,
                      reference_group=args.reference_group)
    frame, ex_design = build_frame(meta, spec, notes)
    excluded += [(s, r, "design metadata") for s, r in ex_design]
    depths = counts[list(frame.df.index)].sum(0)
    depth, kept, why = choose_depth(depths, args.rarefy_depth, args.min_depth, notes)
    excluded += [(s, r, "depth") for s, r in why.items()]
    if len(kept) < len(frame.df):
        frame = frame.subset(kept)
    if len(frame.df) < 6:
        raise InputError(f"Only {len(frame.df)} samples remain; too few for statistical analysis.")
    ids = list(frame.df.index)
    counts = counts[ids]
    counts = counts.loc[counts.sum(1) > 0]
    if lineages is not None:
        lineages = lineages.reindex(counts.index)
    fixed = FixedDesign(frame)
    rank = analyse_rank(fixed)
    checks, cells, audit = validate_design(frame, fixed, rank, [(s, r) for s, r, _ in excluded], notes)
    exdf = pd.DataFrame(excluded, columns=["sample_id", "reason", "stage"])
    qc = out / "qc"
    _w(exdf, qc / "excluded_samples.tsv"); _w(checks, qc / "design_validation.tsv")
    _w(cells.reset_index() if isinstance(cells, pd.DataFrame) else cells, qc / "design_cells.tsv")
    if audit is not None and len(audit):
        _w(audit, qc / "repeated_measures_audit.tsv")
    sample_info = frame.df.copy()
    sample_info["reads"] = depths.reindex(ids)
    _w(sample_info, qc / "analysis_samples.tsv", index=True)
    plots.depth_plot(depths.reindex(ids), depth, qc / "library_sizes.png")
    if depth:
        curves = rarefaction_curves(counts.T.to_numpy(), ids, rng, max_depth=depth)
        _w(curves, qc / "rarefaction_curves.tsv")
        plots.rarefaction_plot(curves, depth, qc / "rarefaction_curves.png")
    if depth is None:
        notes.info("rarefaction", "No rarefaction (default). Observed richness and Jaccard presence/absence are depth-sensitive; "
                                  "Bray-Curtis uses within-sample proportions and Aitchison/CLR use log-ratio representations.")

    levels = resolve_levels(args.taxonomic_level, lineages, notes)
    primary_level = _primary_report_level(args.primary_report_level, levels)
    metrics_a = _split(args.alpha_metrics) or alpha_mod.ALL_METRICS
    metrics_b = _split(args.beta_metrics)
    ext_d = {}
    for item in args.external_distance or []:
        name, _, p = item.partition("=")
        if not p:
            raise InputError("--external-distance expects name=path")
        d = beta_mod.read_distance_matrix(Path(p))
        miss = [i for i in ids if i not in d.index]
        if miss:
            raise InputError(f"External distance '{name}' lacks samples: {miss[:5]}")
        ext_d[name] = d.loc[ids, ids]
    ext_a = None
    if args.external_alpha:
        ext_a = pd.read_csv(args.external_alpha, sep=None, engine="python", index_col=0, dtype={0: str})
        ext_a.index = ext_a.index.astype(str)
        ext_a = ext_a.reindex(ids).apply(pd.to_numeric, errors="coerce")

    # ---------------- rarefaction draws -----------------
    n_draws = max(1, int(args.rarefaction_iterations)) if depth else 1
    C = counts.T.to_numpy()  # samples x features
    maps = {lv: build_level_map(counts.index, lineages, lv) for lv in levels}
    alpha_draws = {lv: [] for lv in levels}
    dist_sum = {lv: {m: 0 for m in metrics_b if m in ("braycurtis", "jaccard")} for lv in levels}
    aggregated_raw = {}
    for lv in levels:
        aggregated_raw[lv] = pd.DataFrame(np.asarray(maps[lv].matrix @ C.T).T, index=ids, columns=maps[lv].labels)
    alpha_raw = {
        lv: alpha_mod.compute_alpha(aggregated_raw[lv].to_numpy(), ids, metrics_a, lv == "feature")
        for lv in levels
    } if not args.skip_alpha else {}
    for d in range(n_draws):
        R = rarefy_matrix(C, depth, rng) if depth else C
        for lv in levels:
            Ag = np.asarray(maps[lv].matrix @ R.T).T
            if not args.skip_alpha:
                adraw = alpha_mod.compute_alpha(Ag, ids, metrics_a, lv == "feature")
                for metric in ("chao1", "goods_coverage"):
                    if metric in adraw and metric in alpha_raw[lv]:
                        adraw[metric] = alpha_raw[lv][metric]
                alpha_draws[lv].append(adraw)
            if not args.skip_beta:
                tab = pd.DataFrame(Ag.T, columns=ids)
                for m in dist_sum[lv]:
                    src = tab.div(tab.sum(0).replace(0, np.nan), axis=1).fillna(0) if m == "braycurtis" else tab
                    dist_sum[lv][m] = dist_sum[lv][m] + beta_mod.distance_from_table(src, m)
    results = dict(levels={}, frame=frame, depth=depth, n_draws=n_draws, checks=checks, cells=cells, audit=audit,
                   excluded=exdf, notes=notes, args=vars(args), version=__version__, levels_run=levels,
                   primary_level=primary_level, rank=rank, fixed=fixed)

    first = True
    for lv in levels:
        L = results["levels"][lv] = dict(first_pass=None, alpha=None, beta={}, taxa=None, labels=list(maps[lv].labels))
        base = out / safe_name(lv)
        _w(maps[lv].table.reset_index().rename(columns={"index": "feature"}), base / "tables" / "feature_to_taxon_map.tsv")
        raw = aggregated_raw[lv]
        _w(raw.T, base / "tables" / f"{lv}_counts.tsv", index=True)
        scope = lv
        # The descriptive first pass is deliberately unrarefied. CPM is a
        # within-sample library-size scaling for summaries/plots, not the
        # inferential model used by the downstream analyses.
        if args.first_pass:
            fpdir = base / "first_pass"
            cpm = first_pass_mod.counts_per_million(raw)
            means = first_pass_mod.grouped_means(cpm, frame)
            _w(cpm.T, fpdir / "cpm.tsv", index=True)
            _w(means, fpdir / "mean_cpm_by_group.tsv")
            group_keys = ["group"] + (["time"] if frame.has_time else [])
            plots.mean_cpm_barplot(means, group_keys, fpdir / "mean_cpm_top_taxa.png", args.top_n_taxa)

            pca_coords, pca_loadings = first_pass_mod.pca(cpm)
            if len(pca_coords.columns) >= 2:
                _w(pca_coords, fpdir / "pca_coordinates.tsv", index=True)
                _w(pca_loadings, fpdir / "pca_loadings.tsv")
                plots.first_pass_ordination(
                    pca_coords, frame, fpdir / "pca.png", f"{lv}: PCA of log2(CPM + 1)",
                    pca_coords.attrs.get("variance_explained"))

            lda_coords, lda_loadings = first_pass_mod.lda(cpm, frame.df["group"])
            if len(lda_coords.columns):
                _w(lda_coords, fpdir / "lda_coordinates.tsv", index=True)
                _w(lda_loadings, fpdir / "lda_loadings.tsv")
                plots.first_pass_ordination(lda_coords, frame, fpdir / "lda.png", f"{lv}: descriptive LDA projection")

            lefse = first_pass_mod.lefse_style(cpm, frame, args.min_prevalence, args.alpha_level)
            _w(lefse, fpdir / "lefse_style_screen.tsv")
            plots.stacked_bar(cpm.div(1_000_000), frame, fpdir / "relative_abundance_top_taxa.png",
                              args.top_n_taxa, "group")
            L["first_pass"] = dict(cpm=cpm, means=means, pca=pca_coords, pca_loadings=pca_loadings,
                                   lda=lda_coords, lda_loadings=lda_loadings, lefse=lefse)
        # ---- alpha ----
        if not args.skip_alpha:
            a = alpha_mod.run_alpha(alpha_draws[lv], frame, fixed, rank, args.alpha_transform, args.p_adjust, notes,
                                    f"{scope}/alpha", ext_a if first else None)
            L["alpha"] = a
            ad = base / "alpha"
            _w(a["mean_values"], ad / "alpha_values.tsv", index=True)
            for k in ("omnibus", "contrasts", "variance", "diagnostics", "descriptives"):
                if len(a[k]):
                    _w(a[k], ad / f"alpha_{k}.tsv")
            for m in a["mean_values"].columns:
                plots.alpha_plot(a["mean_values"][m], frame, m, ad / "figures" / f"alpha_{safe_name(m)}.png", not args.no_connect_repeated)
                if m in a["residuals"]:
                    plots.residual_plot(*a["residuals"][m], m, ad / "figures" / f"diagnostics_{safe_name(m)}.png")
        # ---- beta ----
        if not args.skip_beta:
            dists = {m: s / n_draws for m, s in dist_sum[lv].items()}
            if "aitchison" in metrics_b:
                clr = taxa_mod.clr_transform(raw, args.clr_pseudocount, args.clr_zero_method)
                dists["aitchison"] = beta_mod.distance_from_table(clr.T, "aitchison")
            if first:
                dists.update(ext_d)
            for m, D in dists.items():
                bd = base / "beta" / safe_name(m)
                _w(D, bd / "distance_matrix.tsv", index=True)
                corrected, correction, neg_fraction = beta_mod.lingoes_correct(D.to_numpy(float))
                _w(pd.DataFrame(corrected, index=D.index, columns=D.columns),
                   bd / "inference_distance_matrix.tsv", index=True)
                _w(pd.DataFrame([{"method": "Lingoes" if correction else "none",
                                  "constant": correction,
                                  "original_negative_eigenvalue_fraction": neg_fraction}]),
                   bd / "distance_correction.tsv")
                try:
                    coords, var, negf = beta_mod.pcoa(D)
                    _w(coords, bd / "pcoa_coordinates.tsv", index=True)
                    plots.pcoa_plot(coords, var, frame, bd / "pcoa.png", f"{lv}: {m}", not args.no_connect_repeated)
                except ValueError as e:
                    notes.warn(f"{scope}/beta/{m}", f"PCoA skipped: {e}")
                bres = beta_mod.run_beta(D, frame, fixed, nperm, rng, args.p_adjust, notes,
                                         f"{scope}/beta/{m}", args.change_magnitude)
                L["beta"][m] = bres
                for k, df in bres.items():
                    if isinstance(df, pd.DataFrame) and len(df):
                        _w(df, bd / f"{k}.tsv")
        # ---- taxa ----
        if args.taxa_stats and (lv != "feature" or len(raw.columns) <= 5000):
            t = taxa_mod.run_taxa(raw, pd.Series(maps[lv].unresolved_mask, index=maps[lv].labels), frame, fixed, rank,
                                  args.clr_pseudocount, args.min_prevalence, args.min_mean_relative,
                                  args.include_unclassified_taxa, notes, f"{scope}/taxa", args.clr_zero_method,
                                  args.min_unit_prevalence)
            L["taxa"] = t
            td = base / "taxa"
            _w(t["summary"], td / "taxon_summary.tsv")
            for k in ("omnibus", "contrasts", "failures"):
                if len(t[k]):
                    _w(t[k], td / f"taxa_{k}.tsv")
            _w(t["clr"].T, td / "clr_values.tsv", index=True)
            fg = td / "figures"
            plots.stacked_bar(t["rel"], frame, fg / "top_taxa_by_group.png", args.top_n_taxa, "group")
            plots.stacked_bar(t["rel"], frame, fg / "top_taxa_by_sample.png", args.top_n_taxa, "sample")
            ctd = t["contrasts"]
            if len(ctd):
                keys = ["family", "level_a", "level_b", "time", "group"]
                for i, (k, d) in enumerate(ctd.groupby(keys, dropna=False)):
                    if i >= 24:
                        break
                    lab = " ".join(str(x) for x in k if str(x))
                    plots.volcano_plot(d, f"{lv}: {lab}", fg / f"volcano_{safe_name(lab)}.png", args.alpha_level)
            om = t["omnibus"]
            if len(om):
                top = om.sort_values("p_value").drop_duplicates("taxon")["taxon"].head(30).tolist()
                plots.heatmap(t["clr"], frame, top, fg / "heatmap_top_taxa.png", f"{lv}: top taxa by smallest omnibus p (CLR z-scores)")
        first = False

    sensitivity = None
    if not args.no_rarefaction_sensitivity and depth is None:
        primary = results["levels"][primary_level]
        sensitivity = _run_rarefaction_sensitivity(
            counts, maps[primary_level], primary_level, metrics_a, metrics_b, frame, fixed, rank,
            args, notes, out,
            primary.get("alpha") or {"omnibus": pd.DataFrame()},
            primary.get("beta", {}),
        )
    results["rarefaction_sensitivity"] = sensitivity

    # ---------------- cross-level summary -----------------
    summary = build_summary(results)
    _w(summary, out / "results_summary.tsv")
    results["summary"] = summary
    _w(notes.frame(), out / "warnings.tsv")
    manifest = dict(toolkit_version=__version__, python=platform.python_version(),
                    packages={k: __import__(k).__version__ for k in ("numpy", "pandas", "scipy", "statsmodels", "matplotlib")},
                    command=" ".join(sys.argv), arguments={k: v for k, v in vars(args).items()},
                    inputs={k: {"path": str(v), "sha256": sha256(Path(v))} for k, v in
                            dict(table=args.table, metadata=args.metadata, **({"taxonomy": args.taxonomy} if args.taxonomy else {})).items()},
                    rarefaction_depth=depth, rarefaction_draws=n_draws, samples_analysed=len(ids),
                    levels=levels, primary_report_level=primary_level,
                    rarefaction_sensitivity=bool(sensitivity),
                    rarefaction_sensitivity_depth=sensitivity["depth"] if sensitivity else None,
                    rarefaction_sensitivity_iterations=sensitivity["iterations"] if sensitivity else 0)
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2, default=str))
    from . import report
    report.write_reports(results, out)
    print(f"Done. Results in {out}/ (open report.html)")
    return 0


def build_summary(results) -> pd.DataFrame:
    rows = []
    for lv, L in results["levels"].items():
        a = L["alpha"]
        if a is not None and len(a["omnibus"]):
            for _, r in a["omnibus"].iterrows():
                rows.append(dict(level=lv, analysis="alpha", metric=r["metric"], term=r["term"], statistic=r["statistic"],
                                 p_value=r["p_value"], p_adj_within_analysis=r.get("p_adj_across_metrics", np.nan),
                                 permutations=np.nan, exact=np.nan, min_attainable_p=np.nan, note=r["note"]))
        for m, b in L["beta"].items():
            t = b.get("terms")
            if t is not None and len(t) and "p_value" in t:
                for _, r in t.iterrows():
                    rows.append(dict(level=lv, analysis="beta", metric=m, term=r["term"], statistic=r.get("pseudo_F", np.nan),
                                     p_value=r.get("p_value", np.nan), p_adj_within_analysis=np.nan,
                                     permutations=r.get("permutations", np.nan), exact=r.get("exact", False),
                                     min_attainable_p=r.get("min_attainable_p", np.nan),
                                     note=f"R2={r.get('R2', np.nan):.3f}; min attainable p={r.get('min_attainable_p', np.nan):.3g}; "
                                          + str(r.get("note", ""))))
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    from .contrasts import adjust_p
    df["p_adj_across_levels"] = np.nan
    for _, ix in df.groupby(["analysis", "metric", "term"]).groups.items():
        if len(ix) > 1:
            df.loc[ix, "p_adj_across_levels"] = adjust_p(df.loc[ix, "p_value"].to_numpy(float), "fdr_bh")
        else:
            df.loc[ix, "p_adj_across_levels"] = df.loc[ix, "p_value"]
    return df
