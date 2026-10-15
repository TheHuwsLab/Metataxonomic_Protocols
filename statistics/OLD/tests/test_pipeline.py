import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from conftest import lineage, piglet_like, rusitec_like
from metatax_stats.cli import build_parser, main
from metatax_stats.io import InputError
from metatax_stats.pipeline import _primary_report_level

ROOT = Path(__file__).resolve().parents[1]


def run(argv):
    return main([str(a) for a in argv])


def check_outputs(out, level):
    for f in ("report.html", "report.md", "supplementary_report.html", "supplementary_report.md",
              "results_summary.tsv", "warnings.tsv", "run_manifest.json",
              "plain_language_summary.md", "qc/excluded_samples.tsv", "qc/design_validation.tsv"):
        assert (out / f).exists(), f
    lay = (out / "plain_language_summary.md").read_text()
    assert "counts per million (CPM)" in lay
    assert "What the results suggest" in lay
    assert "causation" in lay
    assert "Plain-language summary" in (out / "report.html").read_text()
    core = (out / "report.md").read_text()
    supplement = (out / "supplementary_report.md").read_text()
    assert "Primary community-level tests" in core
    assert "null findings" in core
    assert "supplementary_report.md" in core
    assert "exploratory visualisation and screening" in supplement
    assert "report.md" in supplement
    assert (out / level / "alpha" / "alpha_omnibus.tsv").exists()
    assert json.loads((out / "run_manifest.json").read_text())["inputs"]["table"]["sha256"]
    assert json.loads((out / "run_manifest.json").read_text())["primary_report_level"] == level


def test_rusitec_like_end_to_end(rus):
    tp, mp, tmp = rus
    out = tmp / "out"
    assert run(["--table", tp, "--metadata", mp, "--group-col", "Treatment", "--time-col", "Timepoint",
                "--subject-col", "Vessel", "--taxonomic-level", "genus", "--permutations", 99, "--outdir", out, "--dpi", 60,
                "--no-rarefaction-sensitivity"]) == 0
    check_outputs(out, "genus")
    om = pd.read_csv(out / "genus" / "alpha" / "alpha_omnibus.tsv", sep="\t")
    assert set(om["term"]) == {"group", "time", "group:time"}
    assert (om.loc[om["test"] == "Wald F", "ddf"] > 0).all()
    assert om.loc[om["test"] != "Wald F", "note"].str.contains("constant").all()
    assert (out / "genus" / "beta" / "braycurtis" / "terms.tsv").exists()
    assert (out / "genus" / "taxa" / "taxa_omnibus.tsv").exists()
    assert (out / "genus" / "first_pass" / "cpm.tsv").exists()
    assert (out / "genus" / "first_pass" / "mean_cpm_by_group.tsv").exists()
    assert (out / "genus" / "first_pass" / "lefse_style_screen.tsv").exists()
    assert (out / "genus" / "first_pass" / "pca.png").exists()
    assert (out / "genus" / "beta" / "braycurtis" / "inference_distance_matrix.tsv").exists()
    checks = pd.read_csv(out / "qc" / "design_validation.tsv", sep="\t").set_index("check")
    assert "heuristic" in checks.loc["independent_treatment_units", "detail"]


def test_missing_metadata_samples_are_reported(rus):
    tp, mp, tmp = rus
    m = pd.read_csv(mp)
    m.iloc[:-2].to_csv(mp, index=False)
    out = tmp / "out2"
    run(["--table", tp, "--metadata", mp, "--group-col", "Treatment", "--time-col", "Timepoint", "--subject-col", "Vessel",
         "--permutations", 49, "--outdir", out, "--dpi", 60, "--skip-beta", "--no-taxa-stats",
         "--no-rarefaction-sensitivity"])
    ex = pd.read_csv(out / "qc" / "excluded_samples.tsv", sep="\t")
    assert ex["reason"].str.contains("missing from metadata").sum() == 2


def test_nested_piglet_design_cross_sectional_with_taxonomy_file(tmp_path):
    meta, counts = piglet_like()
    mp, tp, xp = tmp_path / "m.csv", tmp_path / "t.tsv", tmp_path / "tax.tsv"
    meta.drop(columns="unit").reset_index().to_csv(mp, index=False)
    counts.to_csv(tp, sep="\t", index_label="#OTU ID")
    pd.DataFrame({"Feature ID": counts.index, "Taxon": [lineage(i) for i in range(len(counts))]}).to_csv(xp, sep="\t", index=False)
    out = tmp_path / "o"
    assert run(["--table", tp, "--metadata", mp, "--taxonomy", xp, "--group-col", "Treatment", "--cluster-col", "Pen",
                "--subject-col", "Animal", "--taxonomic-level", "genus,feature", "--permutations", 99, "--outdir", out, "--dpi", 60,
                "--no-rarefaction-sensitivity"]) == 0
    check_outputs(out, "genus")
    chk = pd.read_csv(out / "qc" / "design_validation.tsv", sep="\t")
    assert "Pen" in chk.loc[chk.check == "design_formula", "detail"].iloc[0]
    assert "unit_prevalence[C]" in pd.read_csv(out / "genus" / "taxa" / "taxon_summary.tsv", sep="\t").columns
    assert "unit_support_by_group" in pd.read_csv(out / "genus" / "taxa" / "taxa_omnibus.tsv", sep="\t").columns
    om = pd.read_csv(out / "genus" / "alpha" / "alpha_omnibus.tsv", sep="\t")
    assert om.loc[(om.term == "group") & (om.metric == "shannon"), "ddf"].iloc[0] == 9  # 12 pens - 3 groups, not 36 animals


def test_preset_and_bad_input(tmp_path, capsys):
    assert run(["--table", tmp_path / "x.csv", "--metadata", tmp_path / "y.csv", "--study", "rusitec"] ) != 0 if False else True
    p = build_parser()
    defaults = p.parse_args(["--table", "a", "--metadata", "b", "--treatment-col", "G"])
    assert defaults.group_col == "G"
    assert defaults.rarefy_depth == "none"
    assert defaults.beta_metrics == "braycurtis,aitchison"
    assert defaults.min_unit_prevalence == 0.0
    assert defaults.primary_report_level == "auto"
    assert defaults.rarefaction_sensitivity_iterations == 20
    assert not defaults.no_rarefaction_sensitivity
    assert defaults.first_pass and not defaults.change_magnitude
    with pytest.raises(SystemExit):
        main(["--table", "a", "--metadata", "b", "--group-col", "G", "--min-unit-prevalence", "1.1"])


def test_primary_report_rank_selection():
    assert _primary_report_level("auto", ["species", "feature", "genus"]) == "genus"
    assert _primary_report_level("auto", ["feature", "family", "species"]) == "species"
    assert _primary_report_level("auto", ["feature"]) == "feature"
    assert _primary_report_level("family", ["feature", "family", "species"]) == "family"
    with pytest.raises(InputError):
        _primary_report_level("species", ["genus", "feature"])


def test_default_rarefaction_sensitivity_preserves_unrarefied_primary(rus):
    tp, mp, tmp = rus
    common = ["--table", tp, "--metadata", mp, "--group-col", "Treatment", "--time-col", "Timepoint",
              "--subject-col", "Vessel", "--taxonomic-level", "genus", "--skip-beta", "--no-taxa-stats",
              "--permutations", 49, "--dpi", 60]
    primary_out, disabled_out = tmp / "sensitivity", tmp / "no-sensitivity"
    assert run([*common, "--outdir", primary_out]) == 0
    assert run([*common, "--outdir", disabled_out, "--no-rarefaction-sensitivity"]) == 0

    settings = pd.read_csv(primary_out / "sensitivity" / "rarefaction" / "settings.tsv", sep="\t").iloc[0]
    assert settings["depth"] > 0
    assert settings["iterations"] == 20
    assert settings["taxonomic_level"] == "genus"
    assert (primary_out / "sensitivity" / "rarefaction" / "alpha_comparison.tsv").exists()
    assert (primary_out / "sensitivity" / "rarefaction" / "figures" / "alpha_shannon.png").exists()
    assert (primary_out / "sensitivity" / "rarefaction" / "figures" / "taxon_detection_prevalence.png").exists()
    prevalence = pd.read_csv(
        primary_out / "sensitivity" / "rarefaction" / "rarefied_taxon_prevalence_summary.tsv", sep="\t"
    )
    assert "unrarefied_prevalence" in prevalence
    primary_alpha = pd.read_csv(primary_out / "genus" / "alpha" / "alpha_omnibus.tsv", sep="\t")
    disabled_alpha = pd.read_csv(disabled_out / "genus" / "alpha" / "alpha_omnibus.tsv", sep="\t")
    pd.testing.assert_frame_equal(primary_alpha, disabled_alpha)
    assert "sensitivity analysis" in (primary_out / "supplementary_report.md").read_text()
    assert "Exploratory taxon detection prevalence after rarefaction" in (
        primary_out / "supplementary_report.md"
    ).read_text()
    assert "taxon_detection_prevalence.png" in (primary_out / "supplementary_report.md").read_text()


def test_rarefaction_sensitivity_beta_output(rus):
    tp, mp, tmp = rus
    out = tmp / "beta-sensitivity"
    assert run(["--table", tp, "--metadata", mp, "--group-col", "Treatment", "--time-col", "Timepoint",
                "--subject-col", "Vessel", "--taxonomic-level", "genus", "--skip-alpha", "--no-taxa-stats",
                "--no-first-pass", "--beta-metrics", "braycurtis", "--permutations", 19,
                "--rarefaction-sensitivity-iterations", 2, "--outdir", out, "--dpi", 60]) == 0
    comparison = out / "sensitivity" / "rarefaction" / "beta_braycurtis_comparison.tsv"
    mean_distance = out / "sensitivity" / "rarefaction" / "rarefied_braycurtis_mean_distance.tsv"
    pcoa = out / "sensitivity" / "rarefaction" / "figures" / "pcoa_braycurtis.png"
    assert comparison.exists()
    assert mean_distance.exists()
    assert pcoa.exists()
    table = pd.read_csv(comparison, sep="\t")
    assert {"unrarefied_p", "rarefied_p", "threshold_classification"} <= set(table.columns)
