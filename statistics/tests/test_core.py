import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import pdist, squareform

from conftest import lineage, make_counts, piglet_like, rusitec_like
from metatax_stats import beta as B
from metatax_stats.common import Notes
from metatax_stats.design import DesignSpec, FixedDesign, analyse_rank, build_frame, containment_ddf
from metatax_stats.first_pass import counts_per_million, lefse_style
from metatax_stats.models import contrast, fit_model, wald_term
from metatax_stats.report import _key_findings, _report_cell
from metatax_stats.taxa import run_taxa
from metatax_stats.taxonomy import build_level_map, parse_lineage


def test_parse_lineage_styles():
    g = parse_lineage("d__Bacteria;p__Firmicutes;c__;o__X;f__;g__Gen;s__Gen sp1")
    assert g["class"] is None and g["order"] == "X" and g["species"] == "Gen sp1"
    s = parse_lineage("D_0__Bacteria;D_1__Firmicutes;D_2__Bacilli")
    assert s["class"] == "Bacilli" and s["order"] is None
    assert parse_lineage("Bacteria;Firmicutes;Bacilli")["class"] == "Bacilli"


def test_homonym_genera_stay_separate_and_unresolved_labelled():
    lin = pd.Series({"a": "d__B;p__P1;c__C;o__O;f__F;g__Same;s__x", "b": "d__B;p__P2;c__C;o__O;f__F;g__Same;s__y",
                     "c": "d__B;p__P1;c__C;o__O;f__F;g__;s__"})
    lm = build_level_map(pd.Index(["a", "b", "c"]), lin, "genus")
    assert len(lm.labels) == 3
    assert "F_unclassified" in set(lm.labels)
    assert lm.unresolved_mask.sum() == 1


def _design(seed=0, nv=16):
    meta, _ = rusitec_like(seed)
    spec = DesignSpec("Treatment", "Timepoint", "Vessel", group_order=["CTL", "BCP", "TAN", "SAP"])
    fr, _ = build_frame(meta, spec, Notes(verbose=False))
    fx = FixedDesign(fr)
    return fr, fx, analyse_rank(fx)


def test_reml_matches_statsmodels_and_ddf():
    import statsmodels.formula.api as smf
    fr, fx, rk = _design()
    rng = np.random.default_rng(3)
    ue = dict(zip(fr.df["upper"].unique(), rng.normal(size=16)))
    y = fr.df["upper"].map(ue).to_numpy() + rng.normal(size=len(fr.df)) * 0.5 + (fr.df["time"] == "4h") * 1.0
    fit = fit_model(y.to_numpy(float), fx, rk)
    d = fr.df.assign(y=y.to_numpy())
    m = smf.mixedlm("y ~ C(group, Sum) * C(time, Sum)", d, groups=d["upper"]).fit(reml=True)
    assert np.allclose(fx.X.to_numpy() @ fit.beta, np.asarray(m.model.exog) @ np.asarray(m.fe_params), atol=1e-4)
    assert abs(fit.variance["upper"] - float(np.asarray(m.cov_re)[0, 0])) < 1e-3
    assert wald_term(fit, "group")["ddf"] == 12


def test_reml_matches_statsmodels_for_nested_random_intercepts():
    import statsmodels.api as sm

    meta, _ = piglet_like(seed=11)
    repeated = []
    for sample_id, row in meta.iterrows():
        for time in ("T0", "T1"):
            repeated.append({**row.to_dict(), "sample_id": f"{sample_id}_{time}", "Time": time})
    meta = pd.DataFrame(repeated).set_index("sample_id")
    frame, _ = build_frame(meta, DesignSpec("Treatment", "Time", "Animal", "Pen"),
                           Notes(verbose=False))
    fixed = FixedDesign(frame)
    rank = analyse_rank(fixed)
    rng = np.random.default_rng(17)
    pen_effect = {p: rng.normal(0, 0.8) for p in frame.df["upper"].unique()}
    animal_effect = {s: rng.normal(0, 0.4) for s in frame.df["lower"].unique()}
    y = np.array([
        pen_effect[r.upper] + animal_effect[r.lower] +
        (r.time == "T1") * 0.5 + rng.normal(0, 0.5)
        for r in frame.df.itertuples()
    ])
    fit = fit_model(y, fixed, rank)
    d = frame.df.copy()
    for j in range(fixed.X.shape[1]):
        d[f"x{j}"] = fixed.X.iloc[:, j].to_numpy()
    d["y"] = y
    formula = "y ~ 0 + " + " + ".join(f"x{j}" for j in range(fixed.X.shape[1]))
    reference = sm.MixedLM.from_formula(
        formula, d, groups=d["upper"], re_formula="1",
        vc_formula={"lower": "0 + C(lower)"}
    ).fit(reml=True, method="lbfgs", disp=False)
    assert fit.ok and np.allclose(fit.beta, reference.fe_params.to_numpy(), atol=2e-3)
    assert np.isclose(fit.variance["upper"], float(reference.cov_re.iloc[0, 0]), atol=2e-2)
    assert np.isclose(fit.variance["lower"], float(reference.vcomp[0]), atol=2e-2)
    assert np.isclose(fit.sigma2, float(reference.scale), atol=2e-2)


def test_unit_aware_taxon_support_and_optional_filter():
    meta, _ = piglet_like(seed=13)
    frame, _ = build_frame(meta, DesignSpec("Treatment", None, "Animal", "Pen"),
                           Notes(verbose=False))
    fixed = FixedDesign(frame)
    counts = pd.DataFrame({"common": 50, "rare": 0}, index=frame.df.index)
    counts.loc[frame.df["upper"] == "Pen0", "rare"] = 5
    result = run_taxa(
        counts, pd.Series(False, index=counts.columns), frame, fixed, analyse_rank(fixed),
        1.0, 0.0, 0.0, False, Notes(verbose=False), "taxa",
        min_unit_prevalence=0.5
    )
    summary = result["summary"].set_index("taxon")
    assert summary.loc["common", "units_present"] == 12
    assert summary.loc["common", "tested"]
    assert summary.loc["rare", "units_present"] == 1
    assert not summary.loc["rare", "tested"]
    assert summary.loc["rare", "filter_reason"] == "unit prevalence"
    assert "C=1/4" in summary.loc["rare", "unit_support_by_group"]


def test_single_observation_nested_subject_variance_is_not_separately_fitted():
    meta, _ = piglet_like(seed=3)
    frame, _ = build_frame(meta, DesignSpec("Treatment", None, "Animal", "Pen"),
                           Notes(verbose=False))
    fixed = FixedDesign(frame)
    rank = analyse_rank(fixed)
    fit = fit_model(np.random.default_rng(5).normal(size=len(frame.df)), fixed, rank)
    assert fit.ok
    assert "lower" not in fit.variance
    assert "not separately identifiable" in fit.note


def test_report_marks_monte_carlo_p_value_resolution_limit():
    row = pd.Series({"exact": False, "min_attainable_p": 0.001, "permutations": 999})
    assert _report_cell("p_value", 0.001, row) == "≤ 0.001"
    row["exact"] = True
    assert _report_cell("p_value", 0.001, row) == "0.001"
    findings = _key_findings({
        "args": {"alpha_level": 0.05},
        "summary": pd.DataFrame([{
            "analysis": "beta", "level": "genus", "metric": "braycurtis", "term": "time",
            "p_value": 0.001, "p_adj_across_levels": 0.001, "exact": False,
            "min_attainable_p": 0.001,
        }]),
        "levels": {},
    })
    assert "p ≤ 0.001 (permutation resolution limit)" in findings[0]


def test_ddf_split_plot_16_vessels():
    fr, fx, rk = _design()
    ddf = containment_ddf(fx, rk.kept)
    g_cols = [fx.names.index(c) for c in fx.terms["group"]]
    assert all(ddf[j] == 12 for j in g_cols)


def test_nonestimable_term_is_flagged_not_hidden():
    meta, _ = rusitec_like()
    meta = meta[~((meta.Treatment == "SAP") & (meta.Timepoint == "4h"))]
    fr, _ = build_frame(meta, DesignSpec("Treatment", "Timepoint", "Vessel"), Notes(verbose=False))
    fx = FixedDesign(fr)
    rk = analyse_rank(fx)
    assert rk.dropped_terms


def _dist(fr, k=6, seed=0, effect_group=0.0, time_effect=0.0, inter=0.0):
    rng = np.random.default_rng(seed)
    gi = {g: i for i, g in enumerate(fr.group_levels)}
    up = fr.df["upper"].fillna(pd.Series(fr.df.index, index=fr.df.index))
    ue = {u: rng.normal(size=k) for u in up.unique()}
    X = []
    for sid, r in fr.df.iterrows():
        x = ue[up[sid]] + rng.normal(size=k) * 0.5
        x[0] += effect_group * gi[r["group"]]
        if r["time"] == fr.time_levels[-1]:
            x[1] += time_effect
            x[2] += inter * gi[r["group"]]
        X.append(x)
    X = np.array(X)
    return pd.DataFrame(squareform(pdist(X)), index=fr.df.index, columns=fr.df.index)


def test_splitplot_beta_is_calibrated_and_powerful():
    fr, fx, _ = _design()
    rng = np.random.default_rng(1)
    ps = {"group": [], "group:time": []}
    for i in range(40):
        D = _dist(fr, seed=100 + i)
        t = B.run_beta(D, fr, fx, 99, rng, "holm", Notes(verbose=False), "t")["terms"].set_index("term")
        for k in ps:
            ps[k].append(t.loc[k, "p_value"])
    for k, v in ps.items():
        assert np.mean(np.array(v) <= 0.05) < 0.2, k
    D = _dist(fr, seed=5, effect_group=2.5, inter=2.5)
    t = B.run_beta(D, fr, fx, 499, rng, "holm", Notes(verbose=False), "t")["terms"].set_index("term")
    assert t.loc["group", "p_value"] < 0.05 and t.loc["group:time", "p_value"] < 0.05


def test_interaction_is_directional():
    """Same-magnitude changes in opposite directions between groups must be detected by group:time."""
    fr, fx, _ = _design()
    rng = np.random.default_rng(2)
    D_rows = []
    gi = {g: i for i, g in enumerate(fr.group_levels)}
    r2 = np.random.default_rng(9)
    ue = {u: r2.normal(size=4) for u in fr.df["upper"].unique()}
    X = []
    for sid, r in fr.df.iterrows():
        x = ue[r["upper"]] + r2.normal(size=4) * 0.3
        if r["time"] == fr.time_levels[-1]:
            x[0] += 3 if gi[r["group"]] % 2 == 0 else -3
        X.append(x)
    X = np.array(X)
    D = pd.DataFrame(squareform(pdist(X)), index=fr.df.index, columns=fr.df.index)
    out = B.run_beta(D, fr, fx, 499, rng, "holm", Notes(verbose=False), "t", include_change_magnitude=True)
    assert out["terms"].set_index("term").loc["group:time", "p_value"] < 0.05
    assert out["change_magnitude"]["p_value"].iloc[0] > 0.05  # magnitude identical


def test_interaction_uses_treatment_assignment_units_for_denominator():
    meta, _ = piglet_like()
    # Give each pen two repeated observations and retain three animals per pen.
    repeated = []
    for sid, row in meta.iterrows():
        for time in ("d0", "d7"):
            repeated.append({**row.to_dict(), "sample_id": f"{sid}_{time}", "Time": time})
    m = pd.DataFrame(repeated).set_index("sample_id")
    fr, _ = build_frame(m, DesignSpec("Treatment", "Time", "Animal", "Pen"), Notes(verbose=False))
    fx = FixedDesign(fr)
    rng = np.random.default_rng(4)
    X = rng.normal(size=(len(fr.df), 4))
    d = pd.DataFrame(squareform(pdist(X)), index=fr.df.index, columns=fr.df.index)
    out = B.run_beta(d, fr, fx, 49, rng, "holm", Notes(verbose=False), "test")
    interaction = out["terms"].set_index("term").loc["group:time"]
    assert interaction["n_units"] == 12
    assert interaction["df_den"] == 9


def test_lingoes_correction_preserves_all_corrected_distances():
    binary = np.array([[0, 0, 1, 0, 0], [0, 1, 0, 0, 0], [0, 1, 0, 0, 1],
                       [0, 0, 0, 1, 0], [0, 1, 0, 1, 0], [0, 0, 1, 0, 1]], dtype=float)
    d = squareform(pdist(binary, metric="jaccard"))
    corrected, constant, negative = B.lingoes_correct(d)
    coords, info = B.gower_coords(corrected)
    assert negative > 0
    assert constant > 0
    assert np.allclose(squareform(pdist(coords)), corrected, atol=1e-8)
    assert info["negative_fraction"] < 1e-8


def test_permdisp_recomputes_centroids_under_permutation():
    rng = np.random.default_rng(91)
    small = rng.normal(0, 0.2, size=(12, 3))
    wide = rng.normal(0, 2.0, size=(12, 3))
    coords = np.vstack([small, wide])
    labels = np.array(["A"] * 12 + ["B"] * 12)
    result = B.permdisp(coords, labels, 499, np.random.default_rng(5))
    assert result is not None
    assert result["dispersion_p"] < 0.05


def test_permdisp_null_is_not_systematically_anti_conservative():
    rng = np.random.default_rng(121)
    labels = np.array(["A"] * 10 + ["B"] * 10)
    p_values = []
    for _ in range(30):
        coords = rng.normal(size=(20, 5))
        result = B.permdisp(coords, labels, 99, rng)
        p_values.append(result["dispersion_p"])
    assert np.mean(np.asarray(p_values) <= 0.05) < 0.25


def test_cpm_and_lefse_style_screen_use_library_scaled_values():
    meta, counts = piglet_like()
    # Apply a strong, reproducible group-associated change in a few features.
    counts = counts.T.copy()
    group = meta["Treatment"]
    counts.loc[group == "T1", counts.columns[:4]] *= 5
    cpm = counts_per_million(counts)
    assert np.allclose(cpm.sum(axis=1), 1_000_000)
    fr, _ = build_frame(meta, DesignSpec("Treatment", None, None, "Pen"), Notes(verbose=False))
    result = lefse_style(cpm, fr, min_prevalence=0.0)
    assert not result.empty
    assert {"kw_q", "lda_like_separation", "selected"} <= set(result.columns)


def test_independent_design_freedman_lane():
    meta, _ = rusitec_like()
    meta = meta.drop(columns="Vessel")
    meta["Timepoint"] = "x"
    fr, _ = build_frame(meta, DesignSpec("Treatment"), Notes(verbose=False))
    fx = FixedDesign(fr)
    D = _dist(fr, effect_group=1.5)
    out = B.run_beta(D, fr, fx, 199, np.random.default_rng(0), "holm", Notes(verbose=False), "t")
    assert out["terms"].set_index("term").loc["group", "p_value"] < 0.05
