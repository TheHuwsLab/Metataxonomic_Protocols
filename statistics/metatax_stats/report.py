"""Markdown and self-contained HTML reports."""

from __future__ import annotations

import base64
import html
from pathlib import Path

import numpy as np
import pandas as pd

from .common import fmt_p, safe_name


def _fmt(v):
    if isinstance(v, (float, np.floating)):
        if not np.isfinite(v):
            return ""
        return f"{v:.4g}"
    return "" if v is None else str(v)


def _md_table(df: pd.DataFrame, max_rows: int = 40) -> str:
    if df is None or df.empty:
        return "_none_\n"
    d = df.head(max_rows)
    head = "| " + " | ".join(map(str, d.columns)) + " |\n|" + "---|" * len(d.columns) + "\n"
    body = "".join(
        "| " + " | ".join(_report_cell(c, r[c], r).replace("|", "/") for c in d.columns) + " |\n"
        for _, r in d.iterrows()
    )
    more = f"\n_({len(df) - max_rows} more rows in the TSV files)_\n" if len(df) > max_rows else ""
    return head + body + more


def _html_table(df: pd.DataFrame, max_rows: int = 40) -> str:
    if df is None or df.empty:
        return "<p><em>none</em></p>"
    d = df.head(max_rows)
    h = "<table><thead><tr>" + "".join(f"<th>{html.escape(str(c))}</th>" for c in d.columns) + "</tr></thead><tbody>"
    for _, r in d.iterrows():
        cells = []
        for c in d.columns:
            t = _report_cell(c, r[c], r)
            cls = ""
            cells.append(f"<td{cls}>{html.escape(t)}</td>")
        h += "<tr>" + "".join(cells) + "</tr>"
    h += "</tbody></table>"
    if len(df) > max_rows:
        h += f"<p><em>{len(df) - max_rows} more rows in the TSV files.</em></p>"
    return h


def _report_cell(column, value, row) -> str:
    if (column in {"p_value", "dispersion_p"} or column.endswith("_p")) and pd.notna(value):
        p = float(value)
        prefix = column.rsplit("_p", 1)[0] if column.endswith("_p") and column not in {"p_value", "dispersion_p"} else ""
        exact_value = row.get(f"{prefix}_exact", row.get("exact", False)) if prefix else row.get("exact", False)
        exact = bool(exact_value) if pd.notna(exact_value) else False
        floor = row.get(f"{prefix}_min_attainable_p", row.get("min_attainable_p", np.nan)) if prefix else row.get("min_attainable_p", np.nan)
        if not exact and pd.notna(floor) and p <= float(floor) + 1e-12:
            return f"≤ {float(floor):.4g}"
    return _fmt(value)


def _img(path: Path, out: Path) -> str:
    return str(path.relative_to(out))


def _methods(R) -> str:
    a = R["args"]
    fr = R["frame"]
    spec = fr.spec
    parts = [f"Sequence count tables were analysed with metatax_stats v{R['version']} (Python; numpy, pandas, scipy, statsmodels)."]
    if R["depth"]:
        parts.append(f"Samples were rarefied without replacement to {R['depth']} reads per sample; alpha diversity and "
                     f"Bray-Curtis/Jaccard distances were averaged over {R['n_draws']} independent rarefactions.")
    else:
        parts.append("Counts were not rarefied (the default); alpha richness and unrarefied presence/absence are "
                     "depth-sensitive, while other analyses use their metric-specific transformations.")
    sensitivity = R.get("rarefaction_sensitivity")
    if sensitivity:
        parts.append(f"A separate rarefaction sensitivity analysis compared the unrarefied primary results with "
                     f"{sensitivity['iterations']} draws at {sensitivity['depth']} reads per sample at "
                     f"{sensitivity['level']}; it did not replace the primary results.")
    if a.get("first_pass", True):
        parts.append("The descriptive first pass reports counts per million (CPM = counts/library total x 1,000,000), "
                     "group/time mean CPM, PCA of centered log2(CPM+1), regularized descriptive LDA, and a "
                     "LEfSe-inspired Kruskal-Wallis/rank-sum marker screen on experimental-unit-level CPM. "
                     "These are screening/descriptive outputs, not substitutes for the adjusted inferential models.")
    parts.append("Taxa were aggregated by full lineage at the chosen level(s) (not by name alone); unresolved taxa were "
                 "labelled by their closest classified ancestor.")
    parts.append("Alpha diversity (observed features, Shannon [natural log], Simpson [1-D], Pielou evenness, and, at feature level, "
                 "Chao1 and Good's coverage) was modelled with a Gaussian linear model with sum-to-zero coded "
                 f"{spec.group_col}{' x ' + spec.time_col if fr.has_time else ''}"
                 f"{' + ' + ', '.join(spec.covariates) if spec.covariates else ''} fixed effects")
    if fr.has_units:
        re = " and ".join(([f"{spec.cluster_col}"] if spec.cluster_col else []) + ([f"{spec.subject_col}"] if spec.subject_col else []))
        parts[-1] += (f" and random intercept(s) for {re}" + (" (subject nested in cluster)" if spec.cluster_col and spec.subject_col else "") +
                      ", fitted by restricted maximum likelihood (REML). Fixed-effect tests are Wald F tests and contrasts are "
                      "t tests using containment-based approximate denominator degrees of freedom (the smallest of the "
                      "involved coefficients; no Kenward-Roger or Satterthwaite correction). A nested subject random "
                      "intercept is included only when subjects are observed more than once; with one observation per "
                      "subject its variance is confounded with residual variance and is omitted.")
    else:
        parts[-1] += "; samples were treated as independent (ordinary least squares)."
    parts.append(f"Pairwise contrasts were adjusted by {a['p_adjust']} within each metric and contrast family; omnibus p-values "
                 "were adjusted across metrics by Benjamini-Hochberg.")
    parts.append(f"Beta-diversity ({a['beta_metrics']}) was tested with permutation methods ({a['permutations']} permutations) on "
                 "all axes of Euclidean Lingoes-corrected distances when correction was needed: " + (
                     "Freedman-Lane residual permutation (type II) for independent samples" if not fr.has_units else
                     "restricted permutation split-plot tests (whole-unit permutation for between-unit effects, within-unit permutation for time, and "
                         "whole-unit permutation of time-adjusted data for a custom group x time interaction test; nested samples are averaged "
                         "at assignment-unit x time)") +
                     ". PERMDISP recomputes centroids under each permutation; significant PERMANOVA and PERMDISP together leave "
                     "location and dispersion effects unresolved. Non-Euclidean distances are Lingoes-corrected before inference, so tests use "
                     "the corrected Euclidean representation; the correction magnitude is reported.")
    if a["taxa_stats"]:
        zero_method = a.get("clr_zero_method", "multiplicative")
        zero_text = ("multiplicative replacement of zero proportions (delta is sample-specific, capped to leave positive mass)"
                     if zero_method == "multiplicative" else
                     f"adding a count-scale pseudocount of {a['clr_pseudocount']}")
        parts.append(f"Taxa were analysed after a centred log-ratio transform of unrarefied counts using {zero_text}, with the "
                     "same model as alpha diversity per taxon; q-values are Benjamini-Hochberg across taxa within each term/contrast. "
                     "The CLR denominator uses every aggregated taxon; prevalence filtering is applied to the tested responses afterward. "
                     f"Sample prevalence filtering uses --min-prevalence={a['min_prevalence']}; independent-unit support is reported by group. "
                     f"The optional --min-unit-prevalence={a['min_unit_prevalence']} filter is applied when positive. "
                     "Effect sizes are log2 differences in CLR values, not absolute fold changes.")
    parts.append("Repeated-measures models currently use random intercepts only; they do not model subject-specific time slopes or "
                 "an explicit residual autocorrelation structure.")
    return " ".join(parts)


def _key_findings(R) -> list[str]:
    a = R["args"]["alpha_level"]
    out = []
    S = R["summary"]
    if S is not None and len(S):
        for analysis in ("alpha", "beta"):
            p_column = "p_adj_across_levels"
            p_values = S[p_column].where(S[p_column].notna(), S["p_value"])
            s = S[(S["analysis"] == analysis) & (p_values < a)]
            for _, r in s.iterrows():
                p_label = f"p = {fmt_p(r['p_value'])}"
                exact = r.get("exact", False)
                if analysis == "beta" and not (pd.notna(exact) and bool(exact)):
                    floor = r.get("min_attainable_p", np.nan)
                    if pd.notna(floor) and r["p_value"] <= float(floor) + 1e-12:
                        p_label = f"p ≤ {float(floor):.4g} (permutation resolution limit)"
                adjusted = r.get(p_column)
                extra = f", across-level adjusted p = {fmt_p(adjusted)}" if pd.notna(adjusted) else ""
                out.append(f"{analysis} / {r['level']} / {r['metric']} / {r['term']}: {p_label}{extra}")
        if not out:
            out.append(f"No community-level (alpha/beta) term reached p < {a}.")
    for lv, L in R["levels"].items():
        t = L["taxa"]
        if t and len(t["omnibus"]):
            om = t["omnibus"]
            for term, d in om.groupby("term"):
                out.append(f"taxa / {lv} / {term}: {int((d['q_value'] < a).sum())} of {len(d)} tested taxa with q < {a}")
        if t and len(t["contrasts"]):
            for k, d in t["contrasts"].groupby(["family", "level_a", "level_b", "time", "group"], dropna=False):
                n = int((d["q_value"] < a).sum())
                if n:
                    lab = " ".join(str(x) for x in k if str(x))
                    out.append(f"taxa / {lv} / {lab}: {n} taxa with q < {a}")
    return out


def _lay_summary(R) -> str:
    """Plain-language methods and cautious interpretation for non-specialist readers."""
    frame, args = R["frame"], R["args"]
    n = len(frame.df)
    excluded = len(R["excluded"])
    lines = [
        "# Plain-language summary",
        "",
        "## What was done",
        "",
        f"This analysis compared the microbial DNA sequence profiles in {n} samples. "
        f"{excluded} {'sample was' if excluded == 1 else 'samples were'} excluded before analysis.",
    ]
    if frame.has_units:
        unit_name = frame.spec.cluster_col or frame.spec.subject_col or "independent unit"
        n_units = frame.df["upper"].nunique()
        per_group = frame.df.groupby("group")["upper"].nunique().reindex(frame.group_levels)
        lines.append(
            f"The study design identified {n_units} independent upper-level unit(s) ({unit_name}) "
            f"across the groups: " + ", ".join(f"{g}: {int(v)}" for g, v in per_group.items()) + "."
        )
    else:
        lines.append(
            "No subject or cluster identifier was supplied, so each sample was treated as independent. "
            "If samples came from the same animal, pen, vessel or other shared unit, this design information "
            "needs to be added before interpreting group comparisons."
        )
    lines.extend([
        "",
        "The input remained a table of raw sequence counts; normalized values were calculated for the relevant "
        "analyses rather than written back over the input. For descriptive tables and plots, counts were converted "
        "to counts per million (CPM), not TPM: each sample's counts were divided by its total and multiplied by "
        "1,000,000. This puts sample profiles on a common total scale for display. PCA and the descriptive LDA "
        "projection used log2(CPM + 1). PCA is an unsupervised way to display similarities among samples; LDA "
        "is a supervised group-separation display, not a validated prediction test. These plots and the "
        "LEfSe-inspired marker screen are exploratory, not proof that a taxon causes a change.",
        "",
        "For statistical tests, different calculations were used for different questions: diversity measures "
        "were calculated from the counts; Bray-Curtis compares the proportions of taxa between samples; "
        "Aitchison compares the relative balance of taxa after replacing zero proportions; and taxon tests "
        "compare each taxon's relative log-ratio with the overall composition. Shannon and Simpson are summary "
        "scores describing the number and balance of taxa, not whether microbes are beneficial or harmful. "
        "There is no one normalization that is suitable for every question.",
    ])
    sensitivity = R.get("rarefaction_sensitivity")
    if sensitivity:
        lines.append(
            f"A separate sensitivity check compared the unrarefied primary analysis with {sensitivity['iterations']} "
            f"random subsamples of {sensitivity['depth']} reads per sample at the {sensitivity['level']} level. "
            "It is reported in the supplementary report and does not replace the main results."
        )
    if R["depth"] is None:
        lines.append(
            "Rarefaction (randomly reducing each sample to the same read count) was not used. "
            "Observed richness and presence/absence can be affected by sequencing depth."
        )
    else:
        lines.append(
            f"Samples were rarefied to {R['depth']} reads and {R['n_draws']} draws were averaged for the "
            "analyses that use rarefaction. Chao1 and Good's coverage used the original unrarefied counts."
        )
    lines.append(
        f"Taxon tables report how many independent units contained each taxon, separately by group. "
        f"The optional unit-support filter was set to {args['min_unit_prevalence']}; zero means it was not used."
    )
    if frame.has_time:
        lines.append(
            f"The design compared {frame.spec.group_col} groups across {frame.spec.time_col} time/stage levels."
        )
    else:
        lines.append(f"The design compared groups defined by {frame.spec.group_col}.")
    lines.extend(["", "## What the results suggest", ""])

    summary = R.get("summary", pd.DataFrame())
    primary_levels = [R.get("primary_level", "genus")]
    for level in primary_levels:
        rows = summary[summary["level"] == level] if len(summary) else pd.DataFrame()
        if rows.empty:
            continue
        alpha_rows = rows[(rows["analysis"] == "alpha") & rows["p_adj_across_levels"].notna()]
        beta_rows = rows[(rows["analysis"] == "beta") & rows["p_adj_across_levels"].notna()]
        def signal(terms):
            if not terms:
                return pd.DataFrame()
            return pd.concat([alpha_rows[alpha_rows["term"].isin(terms)],
                              beta_rows[beta_rows["term"].isin(terms)]], ignore_index=True)

        time_signal = signal(["time"])
        group_signal = signal(["group"])
        interaction_signal = signal(["group:time"])
        level_data = R["levels"][level]

        def any_adjusted(df):
            return bool(len(df) and (df["p_adj_across_levels"] < args["alpha_level"]).any())

        if frame.has_time:
            if any_adjusted(time_signal):
                sig = time_signal[time_signal["p_adj_across_levels"] < args["alpha_level"]]
                beta_metrics = sorted(set(sig.loc[sig["analysis"] == "beta", "metric"].astype(str)))
                beta_labels = {"braycurtis": "Bray-Curtis", "aitchison": "Aitchison", "jaccard": "Jaccard"}
                alpha_labels = {
                    "observed_features": "number of detected taxa",
                    "shannon": "Shannon diversity",
                    "simpson": "Simpson diversity",
                    "pielou_evenness": "evenness",
                    "chao1": "estimated richness",
                    "goods_coverage": "estimated sequencing coverage",
                }
                alpha_metrics = sorted(set(sig.loc[sig["analysis"] == "alpha", "metric"].astype(str)))
                findings = []
                if beta_metrics:
                    names = " and ".join(beta_labels.get(m, m) for m in beta_metrics)
                    noun = "comparison" if len(beta_metrics) == 1 else "comparisons"
                    findings.append(f"overall community composition (compared using {names})")
                if alpha_metrics:
                    alpha_names = ", ".join(alpha_labels.get(m, m.replace("_", " ")) for m in alpha_metrics)
                    findings.append(f"diversity measures ({alpha_names})")
                lines.append(
                    f"At the {level} level, the analysis found evidence of a time/stage difference in "
                    f"{' and '.join(findings)}. This describes a difference in the measured samples; it does not "
                    "by itself explain why that difference occurred."
                )
                alpha = level_data.get("alpha")
                shannon_signal = ((time_signal["analysis"] == "alpha") &
                                  (time_signal["metric"] == "shannon") &
                                  (time_signal["p_adj_across_levels"] < args["alpha_level"])).any()
                if (alpha is not None and shannon_signal and len(frame.time_levels) == 2
                        and "shannon" in alpha["mean_values"]):
                    values = alpha["mean_values"]["shannon"].reindex(frame.df.index)
                    time_means = values.groupby(frame.df["time"]).mean().reindex(frame.time_levels)
                    if time_means.notna().all():
                        direction = "higher" if time_means.iloc[1] > time_means.iloc[0] else "lower"
                        lines.append(
                            f"For context, the average Shannon diversity was {time_means.iloc[0]:.3f} at "
                            f"{frame.time_levels[0]} and {time_means.iloc[1]:.3f} at {frame.time_levels[1]} "
                            f"({direction} at the second point). This is a descriptive average across groups."
                        )
            else:
                lines.append(
                    f"At the {level} level, the analysis did not find clear adjusted evidence of an overall "
                    "time/stage difference."
                )
            lines.append(
                f"At the {level} level, the analysis "
                + ("found some adjusted evidence of group differences." if any_adjusted(group_signal)
                   else "did not find clear adjusted evidence of an overall group difference.")
            )
            lines.append(
                f"At the {level} level, the analysis "
                + ("found some adjusted evidence that time/stage changes differed among groups."
                   if any_adjusted(interaction_signal)
                   else "did not find clear adjusted evidence that time/stage changes differed among groups.")
            )
        else:
            lines.append(
                f"At the {level} level, the analysis "
                + ("found some adjusted evidence of group differences." if any_adjusted(group_signal)
                   else "did not find clear adjusted evidence of an overall group difference.")
            )

        taxa = level_data.get("taxa")
        if taxa is not None and len(taxa.get("omnibus", pd.DataFrame())):
            tested = taxa["omnibus"]
            for term in (["time", "group", "group:time"] if frame.has_time else ["group"]):
                d = tested[tested["term"] == term]
                if d.empty:
                    continue
                n_sig = int((d["q_value"] < args["alpha_level"]).sum())
                label = {"time": "time/stage", "group": "group", "group:time": "group-by-time"}[term]
                plural = "taxon" if n_sig == 1 else "taxa"
                lines.append(
                    f"The relative taxon analysis at {level} found {n_sig} {plural} that met the adjusted "
                    f"statistical threshold (q < {args['alpha_level']}) for the {label} test, out of "
                    f"{len(d)} tested. The adjustment accounts for testing many taxa. "
                    "These are relative log-ratio signals, not measurements of absolute microbial load."
                )
                if n_sig:
                    top = d[d["q_value"] < args["alpha_level"]].sort_values("q_value").head(5)
                    names = ", ".join(top["taxon"].astype(str))
                    support_summary = taxa["summary"].set_index("taxon")
                    support_examples = []
                    for taxon in top["taxon"].astype(str).head(2):
                        if taxon in support_summary.index:
                            support_examples.append(
                                f"{taxon} ({support_summary.loc[taxon, 'unit_support_by_group']})"
                            )
                    lines.append(
                        f"Examples with the strongest adjusted statistical evidence at this level: {names}. "
                        "They are not necessarily the most biologically important changes. Their independent-unit "
                        "detection counts by group are shown in the taxon tables; sparse support should temper "
                        "interpretation rather than being treated as equivalent to detection across many units."
                    )
                    if support_examples:
                        lines.append("Example independent-unit detection counts: " + " | ".join(support_examples) + ".")
        fp = level_data.get("first_pass")
        if fp is not None and len(fp.get("lefse", pd.DataFrame())):
            selected = int(fp["lefse"]["selected"].sum())
            lines.append(
                f"The exploratory LEfSe-inspired screen selected {selected} candidate marker(s) at {level}. "
                "This screen is not canonical LEfSe and should not be treated as confirmatory evidence."
            )

    design_checks = R.get("checks", pd.DataFrame())
    if len(design_checks):
        small = design_checks[design_checks["check"] == "independent_treatment_units"]
        if len(small) and (small["status"] == "WARN").any():
            lines.extend([
                "",
                "## Important cautions",
                "",
                small.iloc[0]["detail"],
            ])
    if not frame.has_units:
        lines.extend([
            "",
            "## Important caution",
            "",
            "Because no subject/cluster identifiers were supplied, the analysis cannot account for repeated or "
            "nested samples. If samples came from the same animal, pen or other shared unit, the group-level "
            "statistical evidence may be too strong; confirm the design and rerun with the correct identifiers."
        ])
    lines.extend([
        "",
        "The statistical results are evidence about differences in the measured sequence profiles under the "
        "specified model. They do not establish causation, clinical importance, or absolute microbial abundance. "
        "Interpret the effect estimates and confidence intervals alongside the design, sample sizes, model "
        "diagnostics and full result tables. Results across taxonomic levels are related, not independent replications.",
        "",
        "This plain-language summary is automatically generated to help readers navigate the results. "
        "It is not a substitute for scientific review or a final thesis/manuscript interpretation.",
    ])
    return "\n".join(lines) + "\n"


def _supplementary_blocks(R, out: Path) -> list:
    fr, a = R["frame"], R["args"]
    B = [("h1", "Supplementary statistical report"), ("p", f"metatax_stats v{R['version']}"),
         ("p", "Detailed and exploratory analyses are provided here; the concise core report is available here:"),
         ("link", ("Core report (Markdown)", "report.md")),
         ("link", ("Core report (HTML)", "report.html"))]
    B += [("h2", "1. Design"), ("pre", R["fixed"].frame.spec.describe()), ("table", R["checks"])]
    B += [("h3", "Samples per group (and time)"), ("table", R["cells"].reset_index() if hasattr(R["cells"], "reset_index") else R["cells"])]
    B += [("h2", "2. Sample flow and sequencing depth"),
          ("p", f"{len(fr.df)} samples analysed. Rarefaction depth: {R['depth'] or 'none'}"
                + (f"; {R['n_draws']} rarefactions averaged." if R['depth'] else ".")),
          ("table", R["excluded"]), ("img", out / "qc" / "library_sizes.png"), ("img", out / "qc" / "rarefaction_curves.png")]
    kf = _key_findings(R)
    B += [("h2", "3. Inferential findings"), ("list", kf),
          ("p", "Exploratory PCA/LDA projections, taxa bar plots and the LEfSe-inspired screen are not inferential "
                "tests; they are presented separately within each taxonomic-level section below.")]
    warn = R["notes"].frame()
    B += [("h2", "4. Warnings and caveats"), ("table", warn[warn["severity"] != "INFO"], 200)]
    for lv, L in R["levels"].items():
        B.append(("h2", f"5. Results at level: {lv}"))
        fp = L.get("first_pass")
        if fp is not None:
            B += [("h3", f"{lv}: exploratory visualisation and screening"),
                  ("p", "CPM is a descriptive library-size scaling. Mean CPM and LEfSe-inspired results are exploratory; "
                        "the screen uses experimental-unit means and does not replace the repeated-measures models."),
                  ("table", fp["means"].iloc[:, :min(12, fp["means"].shape[1])], 30),
                  ("table", fp["lefse"].query("selected").sort_values(["kw_q", "lda_like_separation"],
                                                                       ascending=[True, False]) if len(fp["lefse"])
                   else fp["lefse"], 30)]
            for name in ("mean_cpm_top_taxa.png", "relative_abundance_top_taxa.png", "pca.png", "lda.png"):
                B.append(("img", out / lv / "first_pass" / name))
        al = L["alpha"]
        if al is not None:
            B += [("h3", f"{lv}: alpha diversity - omnibus tests"),
                  ("table", al["omnibus"][["metric", "term", "statistic", "df", "ddf", "p_value", "p_adj_across_metrics",
                                           "p_min_over_draws", "p_max_over_draws", "note"]] if len(al["omnibus"]) else al["omnibus"], 60),
                  ("h3", f"{lv}: alpha diversity - contrasts (Holm within metric x family)"),
                  ("table", al["contrasts"][["metric", "family", "time", "group", "level_a", "level_b", "estimate", "ci_low", "ci_high",
                                             "df", "p_value", "p_adj"]] if len(al["contrasts"]) else al["contrasts"], 40),
                  ("h3", f"{lv}: model diagnostics"), ("table", al["diagnostics"])]
            for m in al["mean_values"].columns[:3]:
                B.append(("img", out / lv / "alpha" / "figures" / f"alpha_{m}.png"))
        for m, b in L["beta"].items():
            B.append(("h3", f"{lv}: beta diversity - {m}"))
            B.append(("table", b.get("terms")))
            for k, title in (("group_within_stratum", "Group effect within each time (PERMANOVA + PERMDISP)"),
                             ("pairwise_group", "Pairwise group tests (Holm within stratum)"),
                             ("time_within_group", "Time effect within each group (paired restricted permutation)"),
                             ("change_magnitude", "Magnitude of change (not direction)")):
                if k in b and len(b[k]):
                    B += [("h3", title), ("table", b[k])]
            B.append(("img", out / lv / "beta" / m / "pcoa.png"))
        t = L["taxa"]
        if t is not None and len(t["omnibus"]):
            B.append(("h3", f"{lv}: inferential taxon analysis (CLR + mixed model)"))
            B.append(("table", t["omnibus"].sort_values("p_value")[["taxon", "term", "statistic", "ddf", "p_value", "q_value",
                                                                     "prevalence", "mean_relative_abundance", "units_present",
                                                                     "units_total", "unit_support_by_group"]], 25))
            support_cols = ["taxon", "units_present", "units_total", "unit_prevalence",
                            "max_group_unit_prevalence", "unit_support_by_group", "tested", "filter_reason"]
            support = t["summary"][[c for c in support_cols if c in t["summary"]]]
            B += [("h4", f"{lv}: independent-unit support and taxon filters"), ("table", support, 40)]
            if len(t["contrasts"]):
                B.append(("table", t["contrasts"].sort_values("p_value")[["taxon", "family", "time", "group", "level_a", "level_b",
                                                                          "log2_estimate", "log2_ci_low", "log2_ci_high", "p_value",
                                                                          "q_value", "prevalence", "units_present", "units_total",
                                                                          "unit_support_by_group"]], 25))
            B.append(("img", out / lv / "taxa" / "figures" / "top_taxa_by_group.png"))
            B.append(("img", out / lv / "taxa" / "figures" / "heatmap_top_taxa.png"))
    B += [("h2", "6. Methods (auto-generated)"), ("p", _methods(R))]
    B += [("h2", "7. Cross-level summary of community tests"), ("table", R["summary"], 100)]
    B += [("h2", "8. Files"), ("pre",
          "plain_language_summary.md, results_summary.tsv, warnings.tsv, run_manifest.json\n"
          "qc/: sample flow, design validation, depth, rarefaction\n"
          "<level>/tables, alpha, beta/<metric>, taxa: all tests as TSV with figures in figures/ subfolders")]
    sensitivity = R.get("rarefaction_sensitivity")
    if sensitivity:
        B += [
            ("h2", "9. Rarefaction sensitivity comparison"),
            ("p", f"The primary analyses were unrarefied. This separate sensitivity analysis used {sensitivity['iterations']} "
                  f"without-replacement draws at the minimum retained library size ({sensitivity['depth']} reads per sample), "
                  f"at the primary reporting level ({sensitivity['level']}). It compares threshold classifications at 0.05; "
                  "a difference is a sensitivity signal, not a second confirmatory finding. Alpha figures show unrarefied "
                  "group means beside the spread across rarefaction draws (5th–95th percentile); this spread is not a "
                  "biological confidence interval. Chao1 and Good's coverage remain based on unrarefied counts. "
                  "Aitchison distances, taxon-level tests, descriptive CPM outputs and supplied external metrics are not "
                  "rarefied and are not part of this comparison. Alpha classifications use across-metric adjusted p-values; "
                  "beta classifications use raw permutation p-values and are descriptive only, without a new multiplicity correction."),
            ("h3", "Alpha-diversity comparison (adjusted p-values)"),
            ("table", sensitivity["alpha"], 100),
            ("img", out / "qc" / "rarefaction_sensitivity_curves.png"),
        ]
        for metric in sensitivity.get("alpha_metrics", []):
            B.append(("img", out / "sensitivity" / "rarefaction" / "figures" /
                      f"alpha_{safe_name(metric)}.png"))
        for metric, comparison in sensitivity["beta"].items():
            B += [("h3", f"Beta-diversity comparison: {metric} (raw permutation p-values)"),
                  ("table", comparison, 100),
                  ("p", "PCoA is calculated from the mean pairwise-distance matrix across rarefaction draws; "
                        "it visualizes the sensitivity distance, not a separate inferential test."),
                  ("img", out / "sensitivity" / "rarefaction" / "figures" /
                   f"pcoa_{safe_name(metric)}.png")]
        if sensitivity.get("taxon_prevalence") is not None and len(sensitivity["taxon_prevalence"]):
            B += [
                ("h3", "Exploratory taxon detection prevalence after rarefaction"),
                ("p", "This plot shows taxa whose detection varies most across rarefaction draws, summarized "
                      "within groups and time points at the independent-unit level (or sample level when no unit ID "
                      "was supplied). Teal points and bars show the rarefied mean and 5th–95th percentile; orange "
                      "diamonds show the unrarefied prevalence. This is descriptive presence support only: it is not a differential-abundance test, and taxon-level "
                      "CLR inference remains based on the unrarefied counts."),
                ("img", out / "sensitivity" / "rarefaction" / "figures" / "taxon_detection_prevalence.png"),
            ]
        B += [("p", f"Machine-readable sensitivity tables are in sensitivity/rarefaction/. Beta sensitivity tests use "
                    f"{R['args']['permutations']} permutations on the mean distance matrix across draws.")]
    elif R["depth"] is not None:
        B += [("h2", "9. Rarefaction sensitivity"),
              ("p", "No separate rarefaction sensitivity comparison was run because the primary analysis itself used "
                    "rarefaction. Use an unrarefied run to compare the default rarefaction sensitivity against the "
                    "unrarefied primary analysis.")]
    elif R["args"].get("no_rarefaction_sensitivity"):
        B += [("h2", "9. Rarefaction sensitivity"),
              ("p", "The default rarefaction sensitivity comparison was disabled with --no-rarefaction-sensitivity.")]
    return B


def _significant(df: pd.DataFrame, threshold: float, columns=("p_adj", "q_value", "p_adj_across_metrics")) -> pd.DataFrame:
    if df is None or df.empty:
        return df
    for column in columns:
        if column in df:
            return df[df[column].notna() & (df[column] < threshold)]
    return df.iloc[0:0]


def _core_blocks(R, out: Path) -> list:
    level = R.get("primary_level", "genus")
    L = R["levels"].get(level, {})
    threshold = R["args"]["alpha_level"]
    B = [
        ("h1", "Metataxonomic statistics report: core results"),
        ("p", f"metatax_stats v{R['version']}. Primary taxonomic level: {level}. "
              "The core report retains planned community-level tests, including null findings, and selected adjusted taxon results. "
              "Exploratory and exhaustive outputs are in supplementary_report.md / supplementary_report.html."),
        ("h2", "1. Study design and data quality"),
        ("pre", R["fixed"].frame.spec.describe()),
        ("table", R["checks"]),
        ("h3", "Samples per group (and time)"),
        ("table", R["cells"].reset_index() if hasattr(R["cells"], "reset_index") else R["cells"]),
        ("p", f"{len(R['frame'].df)} samples analysed; {len(R['excluded'])} excluded. "
              f"Rarefaction depth: {R['depth'] or 'none'}"
              + (f"; {R['n_draws']} rarefactions averaged." if R["depth"] else ".")),
    ]
    warn = R["notes"].frame()
    B += [("h3", "Important warnings"), ("table", warn[warn["severity"] != "INFO"], 100)]

    summary = R.get("summary", pd.DataFrame())
    primary = summary[summary["level"] == level] if len(summary) else pd.DataFrame()
    B += [("h2", "2. Primary community-level tests"),
          ("p", "These are the complete alpha- and beta-diversity omnibus results at the selected primary level. "
                "Null results are retained to avoid significance-based selection. Across-level adjusted p-values, where shown, "
                "account for tests at the other requested ranks as well."),
          ("table", primary, 100)]
    if L.get("alpha") is not None:
        alpha = L["alpha"]
        alpha_contrasts = _significant(alpha["contrasts"], threshold, ("p_adj",))
        if len(alpha_contrasts):
            alpha_contrasts = alpha_contrasts.sort_values("p_adj")
        B += [("h3", "Adjusted alpha-diversity contrasts"),
              ("p", f"Showing contrasts with adjusted p < {threshold}; all contrasts and diagnostics are in the supplement and TSV outputs."),
              ("table", alpha_contrasts)]
        for metric in alpha["mean_values"].columns[:3]:
            B.append(("img", out / level / "alpha" / "figures" / f"alpha_{metric}.png"))
    for metric, beta in L.get("beta", {}).items():
        terms = beta.get("terms")
        B += [("h3", f"Beta diversity: {metric}"), ("table", terms)]
        for key, title in (
            ("group_within_stratum", "Adjusted group effects within strata"),
            ("pairwise_group", "Adjusted pairwise group tests"),
            ("time_within_group", "Adjusted time effects within groups"),
        ):
            selected = _significant(beta.get(key, pd.DataFrame()), threshold,
                                    ("p_adj", "q_value", "p_adj_across_metrics"))
            if selected is not None and len(selected):
                pcol = next((c for c in ("p_adj", "q_value", "p_adj_across_metrics") if c in selected), None)
                if pcol:
                    selected = selected.sort_values(pcol)
                B += [("h4", title), ("table", selected)]
        B.append(("img", out / level / "beta" / metric / "pcoa.png"))

    taxa = L.get("taxa")
    if taxa is not None and len(taxa.get("omnibus", pd.DataFrame())):
        omnibus = taxa["omnibus"]
        selected = _significant(omnibus, threshold, ("q_value",))
        support = [c for c in ("taxon", "term", "q_value", "prevalence", "mean_relative_abundance",
                               "units_present", "units_total", "unit_support_by_group") if c in selected]
        B += [("h2", "3. Adjusted taxon-level candidates"),
              ("p", f"{len(selected)} of {len(omnibus)} tested taxon-term results met q < {threshold} at {level}. "
                    "These are candidates for interpretation, not automatic publication claims; they represent relative "
                    "log-ratio differences and should be considered with independent-unit support."),
              ("table", selected.sort_values("q_value")[support], 40)]
        contrasts = _significant(taxa.get("contrasts", pd.DataFrame()), threshold, ("q_value",))
        if len(contrasts):
            cols = [c for c in ("taxon", "family", "time", "group", "level_a", "level_b", "log2_estimate",
                                "log2_ci_low", "log2_ci_high", "q_value", "units_present", "units_total",
                                "unit_support_by_group") if c in contrasts]
            B += [("h3", "Adjusted taxon contrasts"), ("table", contrasts.sort_values("q_value")[cols], 40)]
    B += [
        ("h2", "4. Interpretation and methods"),
        ("p", "The plain-language summary is included above. The primary-level omnibus tests use the study-design-aware "
              "models described in the supplementary methods; multiplicity adjustment is applied within each stated test "
              "family. Repeated or clustered samples are handled according to the supplied identifiers. Results describe "
              "relative sequence profiles and do not establish causation, clinical importance or absolute microbial load. "
              "See the supplementary report for full transformations, permutation restrictions, model diagnostics and limitations."),
        ("h2", "5. Supplementary material and machine-readable outputs"),
        ("p", "See supplementary_report.md / supplementary_report.html for secondary taxonomic ranks, exploratory "
              "PCA/LDA and LEfSe-inspired screens, complete contrasts, diagnostics, detailed taxon tables and additional plots. "
              "The TSV files retain the complete machine-readable results."),
        ("link", ("Supplementary report (Markdown)", "supplementary_report.md")),
        ("link", ("Supplementary report (HTML)", "supplementary_report.html")),
        ("pre", "plain_language_summary.md, results_summary.tsv, warnings.tsv, run_manifest.json\n"
                "supplementary_report.md / supplementary_report.html\n"
                "qc/: sample flow, design validation, depth, rarefaction\n"
                "<level>/tables, alpha, beta/<metric>, taxa: complete tests as TSV with figures in figures/ subfolders"),
    ]
    return B


def _render_report(B: list, out: Path, stem: str, intro: str = "") -> None:
    md, hh = [], []
    if intro:
        md.append(intro)
        hh.append(f"<h2>{html.escape(intro.splitlines()[0].lstrip('# '))}</h2>")
        lay_html = []
        for block in intro.split("\n\n")[1:]:
            text = block.strip()
            if not text:
                continue
            if text.startswith("## "):
                lay_html.append(f"<h3>{html.escape(text[3:])}</h3>")
            else:
                lay_html.append(f"<p>{html.escape(text.replace(chr(10), ' '))}</p>")
        hh.append("<div class='plain-language'>" + "".join(lay_html) + "</div>")
    for blk in B:
        kind, content, *rest = blk
        n = rest[0] if rest else 40
        if kind in ("h1", "h2", "h3", "h4"):
            md.append("#" * int(kind[1]) + f" {content}\n")
            hh.append(f"<{kind}>{html.escape(content)}</{kind}>")
        elif kind == "p":
            md.append(content + "\n"); hh.append(f"<p>{html.escape(content)}</p>")
        elif kind == "pre":
            md.append(f"```\n{content}\n```\n"); hh.append(f"<pre>{html.escape(content)}</pre>")
        elif kind == "list":
            md.append("\n".join(f"- {x}" for x in content) + "\n")
            hh.append("<ul>" + "".join(f"<li>{html.escape(x)}</li>" for x in content) + "</ul>")
        elif kind == "link":
            label, target = content
            md.append(f"[{label}]({target})\n")
            hh.append(f'<p><a href="{html.escape(target, quote=True)}">{html.escape(label)}</a></p>')
        elif kind == "table":
            md.append(_md_table(content, n)); hh.append(_html_table(content, n))
        elif kind == "img":
            if content.exists():
                md.append(f"![{content.name}]({_img(content, out)})\n")
                b64 = base64.b64encode(content.read_bytes()).decode()
                hh.append(f'<figure><img src="data:image/png;base64,{b64}" style="max-width:100%"><figcaption>{content.name}</figcaption></figure>')
    (out / f"{stem}.md").write_text("\n".join(md))
    css = ("body{font-family:sans-serif;max-width:1100px;margin:2em auto;padding:0 1em}table{border-collapse:collapse;font-size:12px;margin:1em 0}"
           "th,td{border:1px solid #ccc;padding:2px 6px}th{background:#eee}pre{background:#f6f6f6;padding:8px;white-space:pre-wrap}"
           "div.t{overflow-x:auto}")
    title = "Supplementary statistical report" if stem.startswith("supplementary") else "Metataxonomic statistics: core results"
    (out / f"{stem}.html").write_text(f"<!doctype html><html><head><meta charset='utf-8'><title>{title}</title><style>{css}</style></head><body>"
                                      + "\n".join(hh) + "</body></html>")


def write_reports(R, out: Path) -> None:
    lay = _lay_summary(R)
    (out / "plain_language_summary.md").write_text(lay)
    _render_report(_core_blocks(R, out), out, "report", lay)
    _render_report(_supplementary_blocks(R, out), out, "supplementary_report")
