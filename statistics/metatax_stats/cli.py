"""Command-line interface."""

from __future__ import annotations

import argparse
import sys

from . import __version__
from .io import InputError
from .pipeline import run

PRESETS = {
    "rusitec": dict(group_col="Treatment", time_col="Timepoint", subject_col="Vessel"),
    "piglet": dict(group_col="Treatment", time_col=None, cluster_col="Pen", subject_col="Animal"),
}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="metataxonomic_stats", description="Design-aware statistics for metataxonomic data "
                                "(16S/18S/ITS; QIIME2 or any count table).")
    p.add_argument("--version", action="version", version=__version__)
    g = p.add_argument_group("inputs")
    g.add_argument("--table", required=True, help="Count table (CSV/TSV, features x samples or samples x features; raw integer counts)")
    g.add_argument("--metadata", required=True, help="Sample metadata (CSV/TSV)")
    g.add_argument("--taxonomy", help="Optional taxonomy file (feature id -> lineage) for ASV/OTU tables")
    g.add_argument("--sample-id-col", help="Metadata sample-ID column (default: first column)")
    d = p.add_argument_group("study design")
    d.add_argument("--study", choices=sorted(PRESETS), help="Preset supplying default column names only (can be overridden)")
    d.add_argument("--group-col", "--treatment-col", dest="group_col", help="Main between-unit factor (e.g. treatment)")
    d.add_argument("--time-col", help="Repeated-measures factor (time point); omit for cross-sectional data")
    d.add_argument("--subject-col", "--random-col", dest="subject_col", help="Repeatedly sampled unit (e.g. vessel, animal)")
    d.add_argument("--cluster-col", help="Higher-level unit that received the treatment (e.g. pen); subject is nested in it")
    d.add_argument("--covariates", help="Comma-separated covariates (numeric or categorical)")
    d.add_argument("--group-order"); d.add_argument("--time-order"); d.add_argument("--reference-group")
    t = p.add_argument_group("taxonomy, depth and filtering")
    t.add_argument("--taxonomic-level", default="auto", help="Comma list of kingdom..species/feature, 'all', or 'auto' (genus, else deepest)")
    t.add_argument("--primary-report-level", default="auto",
                   help="Taxonomic level for core report ('auto': genus if available, otherwise deepest requested level; or an explicit requested level)")
    t.add_argument("--rarefy-depth", default="none", help="'none' (default), 'auto' (minimum depth), or an integer")
    t.add_argument("--min-depth", type=int, help="Drop samples below this depth before choosing the rarefaction depth")
    t.add_argument("--rarefaction-iterations", type=int, default=10)
    t.add_argument("--rarefaction-sensitivity-iterations", type=int, default=20,
                   help="Draws for the default unrarefied-versus-rarefied sensitivity analysis (default: 20)")
    t.add_argument("--no-rarefaction-sensitivity", action="store_true",
                   help="Disable the default rarefaction sensitivity analysis")
    t.add_argument("--exclude-samples", help="Comma-separated sample IDs"); t.add_argument("--exclude-regex")
    t.add_argument("--exclude-taxa-regex", help="e.g. 'Chloroplast|Mitochondria' (matched against lineage and feature id)")
    a = p.add_argument_group("analysis")
    a.add_argument("--alpha-metrics", help="Default: observed_features,shannon,simpson,pielou_evenness,chao1,goods_coverage")
    a.add_argument("--alpha-transform", choices=["auto", "none", "log"], default="auto")
    a.add_argument("--external-alpha", help="Table sample x metric (e.g. Faith PD)")
    a.add_argument("--beta-metrics", default="braycurtis,aitchison",
                   help="Default: braycurtis,aitchison; Jaccard is optional and depth-sensitive without rarefaction")
    a.add_argument("--external-distance", action="append", metavar="NAME=PATH", help="Precomputed distance matrix (e.g. UniFrac); repeatable")
    a.add_argument("--permutations", type=int, default=999)
    a.add_argument("--p-adjust", default="holm", choices=["holm", "fdr_bh", "bonferroni"])
    a.add_argument("--alpha-level", type=float, default=0.05)
    a.add_argument("--taxa-stats", dest="taxa_stats", action="store_true", default=True)
    a.add_argument("--no-taxa-stats", dest="taxa_stats", action="store_false")
    a.add_argument("--first-pass", dest="first_pass", action="store_true", default=True,
                   help="Run descriptive CPM means, PCA, LDA projection, LEfSe-style screen and taxa bars (default)")
    a.add_argument("--no-first-pass", dest="first_pass", action="store_false")
    a.add_argument("--clr-pseudocount", type=float, default=1.0)
    a.add_argument("--clr-zero-method", choices=["multiplicative", "pseudocount"], default="multiplicative",
                   help="CLR zero handling; multiplicative replacement is default, pseudocount mode is for sensitivity")
    a.add_argument("--min-prevalence", type=float, default=0.1, help="Min fraction of samples with the taxon (for taxon tests)")
    a.add_argument("--min-mean-relative", type=float, default=0.0)
    a.add_argument("--min-unit-prevalence", type=float, default=0.0,
                   help="Optional minimum fraction of independent assignment units with the taxon in at least one group; 0 disables this filter")
    a.add_argument("--include-unclassified-taxa", action="store_true")
    a.add_argument("--top-n-taxa", type=int, default=15)
    a.add_argument("--skip-alpha", action="store_true"); a.add_argument("--skip-beta", action="store_true")
    a.add_argument("--change-magnitude", action="store_true",
                   help="Add exploratory magnitude-only repeated-time contrast (not a substitute for group:time)")
    o = p.add_argument_group("output")
    o.add_argument("--outdir", default="metatax_results"); o.add_argument("--seed", type=int, default=12345)
    o.add_argument("--dpi", type=int, default=200); o.add_argument("--no-connect-repeated", action="store_true")
    return p


def main(argv=None) -> int:
    p = build_parser()
    args = p.parse_args(argv)
    for name in ("min_prevalence", "min_unit_prevalence"):
        value = getattr(args, name)
        if not 0 <= value <= 1:
            p.error(f"--{name.replace('_', '-')} must be between 0 and 1")
    if args.rarefaction_sensitivity_iterations < 1:
        p.error("--rarefaction-sensitivity-iterations must be at least 1")
    if args.study:
        for k, v in PRESETS[args.study].items():
            if getattr(args, k) is None:
                setattr(args, k, v)
    if not args.group_col:
        p.error("--group-col is required (or use --study)")
    try:
        return run(args)
    except InputError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2
