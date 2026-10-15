# metatax_stats (v6.2): metataxonomic analysis workflow

This toolkit analyses 16S/18S/ITS and similar amplicon count tables. It runs a **simple descriptive first pass** followed by design-aware alpha, beta and taxon-level inference. The analyses use different representations because they target different quantities; there is no universal normalization that is correct for every analysis.

Run from the `statistics/` directory:

```bash
python -m pip install -r requirements.txt
python metataxonomic_stats.py \
  --table counts.csv --metadata metadata.csv \
  --group-col Treatment --time-col Timepoint --subject-col Vessel \
  --outdir results
```

Open `results/report.html` for the concise core report. Detailed contrasts, diagnostics, secondary taxonomic ranks, and exploratory outputs are in `results/supplementary_report.html`. For all CLI options, use `python metataxonomic_stats.py --help`.

## Input files

**Count table:** CSV/TSV with non-negative integer raw counts; features may be rows or columns. The toolkit detects orientation by matching sample IDs. The first column contains row/feature identifiers. Taxonomy may be encoded as lineage strings in IDs/table columns or supplied using `--taxonomy`.

**Metadata:** one row per sample, with a sample ID column (first column unless `--sample-id-col` is specified) and all design columns named in the command.

Do not provide relative abundances, CPM, TPM, or previously normalized values as the count table. Table/metadata mismatches, explicit exclusions, missing design metadata and depth-based exclusions are recorded in `qc/excluded_samples.tsv`.

## Define the experimental design

| Study design | Arguments |
|---|---|
| Independent cross-sectional groups | `--group-col Treatment` |
| Repeated observations from one vessel/animal | `--group-col Treatment --time-col Timepoint --subject-col Vessel` |
| Treatment assigned to pens; animals sampled within pens | `--group-col Treatment --cluster-col Pen --subject-col Animal` |
| Animals nested in pens, measured over time | `--group-col Treatment --cluster-col Pen --subject-col Animal --time-col Day` |

Treatment effects are replicated at the unit that received the treatment. If treatment was allocated to pens, pens—not animals within pens—are the independent treatment units. `--study rusitec` and `--study piglet` only supply default metadata column names; they do not infer missing design information.

The QC report separates samples from upper-level independent treatment units. Its small-unit warning is triggered below five units in any group as a screening heuristic only—not a statistical validity or adequacy cutoff. Interpret the actual allocation, unit counts, effect sizes and permutation resolution; neither four nor five units automatically determines whether a design is adequate.

Use `--group-order`, `--time-order`, and `--reference-group` to control factor order and reference contrasts. `--covariates Age,Weight` adds covariates to alpha/taxon models. Current repeated-unit beta tests do not adjust for covariates.

## Shared processing

1. Read and validate raw counts, read metadata and align samples by sample ID.
2. Apply requested sample/taxon exclusions. `--exclude-taxa-regex` matches feature IDs and available lineage strings.
3. Exclude zero-read samples and samples missing required design metadata. Rarefaction is disabled by default; if requested, depth-based sample exclusions are also recorded.
4. Check design cells, fixed-effect estimability, group allocation and repeated-unit completeness; outputs are in `qc/`.
5. Aggregate features to requested ranks using the full lineage. Homonymous terminal names with different parents remain distinct. Unresolved taxa are assigned to the closest resolved ancestor with an `_unclassified` suffix. See each level's `tables/feature_to_taxon_map.tsv`.
6. Run the CPM descriptive first pass, then the inferential analyses. CPM, rarefaction, relative abundance, presence/absence and CLR are distinct processing choices, not one interchangeable normalization.

## Default first pass: CPM, means, PCA, LDA and marker screen

The first pass runs at each requested taxonomic level and is saved to `<level>/first_pass/`.

* **CPM:** counts are divided by each sample's total and multiplied by 1,000,000. This is relative abundance on a CPM scale—not TPM and not an inferential count model.
* **Means:** mean CPM is reported for each group and, when applicable, group × time. These are descriptive sample means; the independent experimental units determine inferential replication.
* **Taxa bar plots:** top taxa are selected by mean CPM. The mean CPM plot displays group means; the stacked relative-abundance plot displays within-sample proportions.
* **PCA:** centered `log2(CPM + 1)` values, without unit-variance scaling of taxa, are used for standard PCA. PCA is an unsupervised visualization and has no attached p-value.
* **LDA:** a regularized Fisher linear-discriminant projection of `log2(CPM + 1)` is output as sample coordinates and taxon loadings. It is supervised and descriptive, not a validated prediction model; no cross-validation accuracy or inferential p-value is claimed.
* **LEfSe-inspired marker screen:** CPM values are averaged to one value per upper experimental unit within each time stratum. It applies Kruskal-Wallis tests, BH adjustment across tested taxa within each stratum, and a pairwise rank-sum consistency gate with Holm adjustment. Candidate features are ranked by a standardized between-group separation score.

This is **not the canonical LEfSe implementation**: its effect score is not the original LEfSe LDA effect-size score, it has no subclass-consistency option, and it does not adjust for covariates. It is a discovery screen, not the primary inferential result. Confirm candidates using the model-based analyses and independent biological evidence. Disable the first pass with `--no-first-pass`.

## Rarefaction

Rarefaction is **off by default** (`--rarefy-depth none`) to retain observed reads. To request it, use `--rarefy-depth auto` (minimum remaining library size) or specify an integer depth. `--min-depth` removes shallower libraries before the depth is selected. Subsampling is without replacement by a multivariate hypergeometric draw; `--rarefaction-iterations` sets the number of draws and `--seed` controls reproducibility.

When enabled, the workflow averages draw-specific alpha metrics other than Chao1/Good's coverage, and Bray-Curtis/Jaccard distance matrices. Chao1 and Good's coverage always use full unrarefied feature counts. Alpha draw-wise minimum and maximum p-values, when present, are stochastic sensitivity summaries—not confidence intervals or the primary p-value.

By default, a **separate, small rarefaction sensitivity analysis** also runs alongside the unrarefied primary analysis: 20 without-replacement draws at the minimum retained sample depth. It reports primary-level alpha-diversity comparisons and, when selected, Bray-Curtis/Jaccard comparisons using the average distance matrix across draws. Supplementary plots compare unrarefied alpha means with rarefied draw variability, show PCoA of the mean rarefied beta-distance matrix, and summarize taxa whose detection prevalence varies most across draws at the independent-unit level (or sample level where unit IDs are unavailable), against their unrarefied prevalence. These are descriptive robustness displays; the prevalence plot is not a differential-abundance test. The sensitivity analysis does not replace the primary unrarefied results or constitute a second confirmatory analysis. Chao1 and Good's coverage remain unrarefied. Aitchison, inferential taxon-level tests, CPM descriptions and external metrics are not part of the sensitivity comparison. Results and figures are in `supplementary_report` and `sensitivity/rarefaction/`. Alpha/beta ordination and prevalence figure intervals summarize rarefaction draw variability, not biological confidence intervals. Adjust the draw count with `--rarefaction-sensitivity-iterations` or disable it with `--no-rarefaction-sensitivity`. If the primary run itself uses `--rarefy-depth`, this additional comparison is skipped.

Rarefaction discards reads and does not equalize precision. The sensitivity depth is set by the shallowest retained library and may therefore be very low; the report warns when it is below 10% of the median. Observed richness and unrarefied Jaccard are depth-sensitive. Review `qc/library_sizes.png`, rarefaction curves when generated, and the exclusion table.

## Normalization and computation by analysis

### Alpha diversity

Metrics are calculated after taxonomic aggregation. By default, they use unrarefied counts. With explicit rarefaction, observed richness, Shannon, Simpson and Pielou are calculated for each rarefied draw; Chao1 and Good's coverage are calculated on full, unrarefied feature counts.

| Metric | Computation |
|---|---|
| Observed features | Number of taxa/features with count greater than zero. Depth-sensitive. |
| Shannon | `-sum(p * ln(p))`, where `p = count / sample total`; natural logarithm. |
| Simpson | `1 - sum(p^2)` using within-sample proportions. |
| Pielou evenness | Shannon divided by `ln(Sobs)`; undefined if `Sobs <= 1`. |
| Chao1 | Bias-corrected `Sobs + F1*(F1-1)/(2*(F2+1))`; feature level, unrarefied. |
| Good's coverage | `1 - F1/N`; feature level, unrarefied. |

Observed richness and Chao1 are log-transformed by default (`--alpha-transform auto`); use `none` or `log` to change this. A Gaussian REML model fits sum-coded group, time, group × time and covariates, with random intercepts for specified units. A nested subject intercept is included only when subjects are observed more than once; with one observation per subject its variance is not separately identifiable from residual variance and is omitted with a QC note. It reports Wald F tests, model-based contrasts and confidence intervals, residual diagnostics, variance components and ICC. Bounded Simpson/Pielou responses use a Gaussian approximation; assess diagnostics and response ranges. The repeated-measures model uses random intercepts only; it does not fit subject-specific time slopes or residual autocorrelation.

Contrasts are adjusted within metric and contrast family using `--p-adjust`; omnibus metric tests are BH-adjusted within term and taxonomic level. Denominator degrees of freedom use containment approximation, not Kenward-Roger/Satterthwaite.

### Beta diversity, PCA and PCoA

Default beta metrics are Bray-Curtis and Aitchison. Jaccard is optional (`--beta-metrics braycurtis,aitchison,jaccard`) because unrarefied detection is depth-sensitive.

| Distance | Data representation |
|---|---|
| Bray-Curtis | Rarefied counts if requested, otherwise raw counts converted to within-sample relative abundance (`count / sample total`). Matrices are averaged across rarefaction draws if used. |
| Jaccard | Presence/absence (`count > 0`) from each rarefaction draw or from unrarefied counts; draw-specific matrices are averaged. |
| Aitchison | Unrarefied counts converted to relative abundance, zeros multiplicatively replaced by default, CLR-transformed, then Euclidean distance is computed. |
| External distance | Matrix supplied by the user (e.g. UniFrac) is used as provided; provenance and normalization are the user's responsibility. It is analyzed at the first requested taxonomic level only. |

The separate first-pass PCA uses centered log2 CPM. Beta PCoA is based on the selected beta distance. For non-Euclidean distances, Lingoes correction is applied before PCoA and permutation inference; both the original and corrected distance matrices and correction information are saved.

Independent designs use Freedman-Lane residual permutations. For repeated designs, whole treatment-assignment units are permuted for between-unit group effects, time labels are permuted within complete repeated units, and the directional group × time test compares time profiles after removing the pooled time profile and permutes labels at the treatment-assignment-unit level. This is a custom split-plot permutation procedure; it is not a drop-in implementation of a named package's PERMANOVA. Its operating characteristics have simulation tooling, but applicability beyond simulated design structures should not be assumed. Its denominator df are based on upper-level independent units, not nested subjects.

In repeated/clustered beta tests, multiple rows in an upper assignment-unit × time cell are averaged before inference. With animals nested in pens, inference is therefore at the pen-level trajectory; animal-level variation within pens is not retained in that beta test. Incomplete upper units are excluded from repeated beta tests and reported. Averaging is appropriate only when observations within the cell are subsamples/replicates representing the same assignment unit, not genuine independent treatment replicates.

PERMDISP accompanies group comparisons and recomputes group centroids/distances for every permuted allocation. A significant PERMANOVA with significant PERMDISP does not distinguish centroid shifts from dispersion effects; PERMDISP diagnoses this ambiguity but does not correct it. Outputs include pseudo-F, R-squared, permutation p-value and minimum attainable p-value.

For Monte Carlo tests, the minimum possible reported p-value is `1 / (permutations + 1)`; report tables and key findings display this boundary as `p ≤ ...` when the observed value reaches it. Exact restricted tests are labelled separately. More permutations improve Monte Carlo resolution but cannot overcome a small number of independent treatment units or too few distinct label arrangements.

### Taxon-level compositional differential analysis

Taxa are aggregated from unrarefied counts. The default zero method is multiplicative replacement of zero proportions: within each sample, zeros receive a small sample-specific replacement based on the minimum positive proportion, capped to preserve positive mass. `--clr-zero-method pseudocount` selects the alternative legacy count-scale pseudocount approach (`--clr-pseudocount`, default 1) for sensitivity analysis.

The CLR is:

`CLR_ij = ln(replaced proportion_ij) - mean_j[ln(replaced proportion_ij)]`

All taxa at the level contribute to the CLR geometric-mean denominator. Prevalence/abundance filtering is applied **after** CLR calculation and controls which taxon responses are modeled, not the denominator. Taxa with prevalence below `--min-prevalence` (default 0.10) or mean relative abundance below `--min-mean-relative` (default 0) are not tested. Unresolved taxa are excluded unless requested. Filters are summarized in `taxon_summary.tsv`.

The taxon summary also reports the number and fraction of **independent assignment units** with each taxon detected, both overall and separately by group. Repeated samples within a unit count once; if treatment was assigned to pens, pen is the support unit. This support information helps identify sparse findings but does not itself establish adequate replication. By default it does not impose a minimum unit-support filter. To require detection in a specified fraction of assignment units in at least one group, use `--min-unit-prevalence 0.25` (valid range 0–1); the default `0` disables this additional filter. The report includes the exact per-group support for leading q-value candidates rather than treating rare and widely observed taxa as equally supported.

Each taxon's CLR is fitted with the alpha REML model. This is a **compositional log-ratio analysis**, not a test of absolute abundance. Effects are converted from natural-log CLR differences to log2 units by dividing by `ln(2)`; they are not exact fold changes. Sparse-taxon results can be sensitive to zero handling. BH q-values are adjusted across taxa within each term/contrast family separately at each taxonomic level. Results across levels are not one multiplicity-controlled family.

Relative abundance for bar plots is `count / sample total`, not the taxon-model response. CPM plots are exactly that proportion multiplied by one million.

### External alpha metrics

`--external-alpha file.csv` accepts supplied sample-level numeric metrics (e.g. Faith's PD calculated elsewhere). Values are not recalculated by this toolkit for each taxonomy rank and are not rarefied or scaled by the input step. The optional alpha transform may log-transform them. They are stored under the first requested level for organization only.

## Multiplicity and interpretation

`--p-adjust` selects Holm (default), BH or Bonferroni for relevant within-family adjustments. Alpha omnibus tests are BH-adjusted across metrics within term and level. Taxon tests are BH-adjusted across taxa within term/contrast family and level. Beta summary results also report across-level BH adjustment. Testing many ranks increases the opportunity for chance findings; the first-pass marker screen is separately exploratory. Interpret adjusted results alongside effects, confidence intervals, experimental-unit counts, design checks and diagnostics.

## Outputs

| Path | Contents |
|---|---|
| `report.html`, `report.md` | Concise core report: design/QC, all primary-level alpha/beta omnibus results (including null results), selected adjusted taxon candidates with independent-unit support, and methods/limitations. |
| `supplementary_report.html`, `supplementary_report.md` | Full detailed report, including secondary ranks, exploratory PCA/LDA and LEfSe-inspired screens, complete contrasts, diagnostics and detailed taxon outputs. |
| `plain_language_summary.md` | Standalone plain-language account of the processing, main findings and design cautions. |

The primary reporting rank defaults to genus when genus is among the requested levels; otherwise it selects the deepest requested taxonomic rank (or feature level if no taxonomic rank is available). Set `--primary-report-level` to an explicitly requested level to override this. Statistical significance alone does not determine whether planned primary-level community tests are included: null results remain in the core report, while exhaustive and exploratory material is routed to the supplementary report. The core report is a navigation and reporting aid, not an automated claim that a result is publication-ready.
| `results_summary.tsv` | Cross-level alpha/beta summary. |
| `warnings.tsv` | Warnings and run notes. |
| `run_manifest.json` | Arguments, versions, input hashes, seed and sample/level summary. |
| `qc/` | Exclusions, design checks/cells, repeated-measures audit, library sizes and rarefaction curves. |
| `sensitivity/rarefaction/` | Default rarefaction robustness check: settings, alpha comparisons and figures, mean rarefied beta distances/PCoA, beta comparisons, and descriptive taxon detection prevalence tables/plot. |
| `<level>/first_pass/` | CPM, group means, PCA/LDA coordinates and loadings, LEfSe-inspired screen and taxa plots. |
| `<level>/tables/` | Aggregated counts and feature-to-taxon map. |
| `<level>/alpha/` | Alpha values, tests, contrasts, diagnostics and plots. |
| `<level>/beta/<metric>/` | Raw/corrected distances, PCoA and permutation test tables. |
| `<level>/taxa/` | CLR values, filters, mixed-model tests, contrasts and figures. |

## Reading the results

| Question | Where to look | How to read it |
|---|---|---|
| Is the data usable? | `qc/`, `warnings.tsv` | Check excluded samples, library sizes, unit counts per group and the small-unit warning first. |
| Do overall communities differ? | `results_summary.tsv` (beta rows), `<level>/beta/<metric>/terms.tsv` | `group` = treatment effect averaged over time; `time` = change over time averaged over groups; `group:time` = whether time trajectories differ between groups. R² is the share of distance variation attributed to the term. |
| Does within-sample diversity differ? | `<level>/alpha/alpha_omnibus.tsv`, `alpha_contrasts.tsv` | Use `p_adj_across_metrics` for omnibus tests and `p_adj` for contrasts. Richness contrasts are on the log scale; `ratio` is the multiplicative difference (first level / second level). |
| Which taxa differ? | `<level>/taxa/taxa_omnibus.tsv`, `taxa_contrasts.tsv`, `taxon_summary.tsv` | Use `q_value` (FDR). Effects are CLR differences in log2 units, i.e. changes relative to the rest of the community, not absolute abundance. Always check `unit_support_by_group`. |
| Do conclusions depend on depth/rarefaction? | `sensitivity/rarefaction/` | `threshold_classification` shows whether each alpha/beta call changed at 0.05. |

Practical guidance:

* **Interaction first.** Because models use sum coding, a `group` main effect is an average over time points and `time` an average over groups. If `group:time` is credible, interpret contrasts within time or within group instead of main effects.
* **Omnibus before contrasts.** Treat pairwise/contrast tables as follow-up to a meaningful omnibus result, or label them exploratory. Contrast adjustment is within a family only; it does not control error across all families and taxa.
* **Absence of evidence is not evidence of absence.** With few independent units, non-significant results mainly reflect limited power. Report effect sizes and confidence intervals alongside p-values.
* **Presence/absence-driven taxa.** A taxon detected in only one group (e.g. 0/4 vs 4/4 units) yields a very large CLR difference that depends on the zero-replacement value. Treat it as a detection difference, not a quantified fold change.
* **Different time points are not different replicates.** Samples from the same vessel/animal/pen are correlated; only upper-level independent units count as replicates.
* **Taxonomic levels are not independent tests.** The same data are analysed at every rank; choose the primary level in advance (default genus) and use ranks above/below it as supporting context.
* **Raw p-values from permutation tests are bounded** by `1/(permutations+1)`; report them as `p ≤ 0.001`, not `p = 0.001`.

## Reporting in a paper

Suggested methods wording (edit to match the options actually used; `run_manifest.json` records versions, arguments and seed):

> Count tables were analysed with metatax_stats v6.2 (Python; numpy, pandas, scipy, statsmodels). Samples were not rarefied for primary analyses; a sensitivity analysis with 20 rarefied draws at the minimum library size was run. Alpha diversity (observed features, Shannon, Simpson, Pielou, Chao1, Good's coverage) was computed from raw counts after aggregation to the stated taxonomic rank and analysed with Gaussian REML models (sum-coded treatment, time and their interaction; random intercept for the repeatedly sampled unit; richness metrics log-transformed; Wald F tests with containment degrees of freedom). Beta diversity (Bray-Curtis on relative abundances; Aitchison on CLR-transformed data after multiplicative zero replacement) was tested by permutation (999 permutations; Freedman-Lane for independent samples; restricted split-plot permutations for repeated measures), with PERMDISP as a dispersion check. Taxon-level differences used CLR-transformed data in the same mixed models, with Benjamini-Hochberg correction across taxa within each term. Taxa with prevalence < 10% of samples were not tested. Results at multiple taxonomic ranks are descriptive of the same data and were not treated as independent replication.

Report at minimum: number of independent units per group, the primary taxonomic level, the adjustment method, effect sizes with confidence intervals, and exact/bounded p-values. Describe the PCA/LDA/LEfSe-inspired outputs as exploratory.

## Limitations

* Small numbers of independent clusters make treatment inference unstable; nested sample counts do not replace cluster replication.
* The fewer-than-five-unit QC message is a heuristic warning, not a universal adequacy threshold. Permutation resolution depends on distinct unit-level label arrangements.
* Containment df are approximate. Random-intercept-only models may not represent complex temporal covariance. Use the included benchmark script to compare REML estimates with statsmodels MixedLM for simulated structures; agreement in those scenarios does not establish correctness for every design.
* Repeated beta tests use a custom restricted split-plot procedure; simulation validation should be run for the designs of interest. They require complete units for time/interaction, average samples within assignment-unit × time cells, and do not adjust for covariates.
* The first-pass LDA/LEfSe-inspired screen is exploratory, not canonical LEfSe and not a cross-validated classifier.
* CLR taxon tests depend on zero replacement and compositional denominator choices; they do not estimate absolute abundance. Unit-level taxon support is reported; no universal support threshold is imposed by default.
* The LEfSe-inspired screen requires every pairwise rank-sum comparison to pass a Holm-adjusted threshold. With the exact two-sided Mann-Whitney test, 4 vs 4 units has a minimum p of 0.029 (0.17 after Holm over the 6 pairs of four groups), so **no taxon can be selected with four groups of four units**; with 5 per group the minimum Holm-adjusted p is 0.048, so selection needs complete separation between every pair. A result of "0 candidates" in small studies is therefore expected and uninformative.
* With few units, exact restricted permutation tests have coarse resolution (e.g. 4 vessels give a minimum exact p of 0.125 for a within-group time test). Use pooled or model-based tests for such questions; check `min_attainable_p`.
* Alpha omnibus/contrast Wald tests assume approximately normal, homoscedastic residuals; check `alpha_diagnostics.tsv` and the diagnostics plots.
* If time points come from the *same* animals/pens/vessels, supply `--subject-col`/`--cluster-col`; without them every sample is treated as independent and p-values may be anti-conservative (the report warns about this).
* External metrics are trusted as supplied. Chao1 and Good's coverage are only computed at feature level.

## Method-validation utilities

The validation scripts are reproducible checks, not a guarantee of general validity. They write scenario-level data for interpretation and should be expanded to match the experimental design and effect sizes being studied.

```bash
python validation/validate_beta.py --replicates 200 --permutations 999
python validation/validate_models.py --replicates 100
```

The beta script estimates type-I error and power for group, time and group × time tests in independent, repeated and nested clustered simulations, with Wilson intervals. Independent designs use Freedman-Lane term tests; repeated and nested designs use restricted unit/time permutations. It writes a fit-level TSV and a Markdown summary. The REML script compares fixed-effect estimates, standard errors and variance components with `statsmodels.MixedLM` for balanced/unbalanced, nested and boundary-variance scenarios. It tries multiple reference optimizers and chooses the converged fit with the highest REML likelihood, then writes a TSV and Markdown summary. Containment-based denominator df are an additional approximate inference step and are not expected to match statsmodels' asymptotic tests. Review convergence/failure rows and simulation uncertainty; do not interpret a small run as proof of calibration.
