#!/usr/bin/env python3
"""
Compare normalised top 10 taxa across multiple QIIME2 taxonomy CSVs.

Uses the `load_and_normalize()` function from the user's existing script.

Expects CSVs that contain:
    feature-id   abundance_column
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# ----------------------------------------------------------------------
# IMPORT EXISTING FUNCTIONS DIRECTLY FROM PROVIDED SCRIPT
# (Copied here verbatim so the script is self-contained)
# ----------------------------------------------------------------------

def tmm_normalize(counts):
    counts = counts.copy()
    counts = pd.to_numeric(counts, errors='coerce')
    counts = counts.fillna(0)
    counts_calc = counts.replace(0, np.nan)
    lib_size = counts.sum()
    if lib_size == 0:
        return counts
    norm_factor = lib_size / 1e6
    return counts / norm_factor


def load_and_normalize(filepath, feature_col='feature-id'):
    import pandas as pd
    import numpy as np

    df = pd.read_csv(filepath, index_col=0)

    first_row_numeric = (
        pd.to_numeric(df.iloc[0], errors='coerce').notna().sum()
        / len(df.columns) > 0.5
    )

    if first_row_numeric:
        df = df.T
        print(f"  Detected transposed format in {filepath}, transposing...")

    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    df = df.fillna(0)

    for col in df.columns:
        df[col] = tmm_normalize(df[col])

    df = df.reset_index()
    df = df.rename(columns={'index': feature_col})

    return df

# ----------------------------------------------------------------------
# HELPER FUNCTIONS
# ----------------------------------------------------------------------

def extract_taxa(df, feature_col, abundance_col):
    """
    Convert a loaded dataframe into a standardised:
        taxon | abundance
    table.

    `feature_col` → becomes 'taxon'
    `abundance_col` → becomes 'abundance'
    """
    if feature_col not in df.columns:
        raise ValueError(f"Missing feature column: {feature_col}")

    if abundance_col not in df.columns:
        raise ValueError(f"Missing abundance column: {abundance_col}")

    out = df[[feature_col, abundance_col]].copy()
    out = out.rename(columns={feature_col: "taxon", abundance_col: "abundance"})

    return out


def get_top10(df):
    """Extract top-10 taxa by abundance."""
    df2 = df.sort_values("abundance", ascending=False).reset_index(drop=True)
    return df2.head(10)


# ----------------------------------------------------------------------
# PLOTTING
# ----------------------------------------------------------------------

def create_plots(gold_df, test_data, output_prefix):
    pdf_path = f"{output_prefix}_top10_taxa_plots.pdf"

    sns.set_style("whitegrid")
    plt.rcParams['figure.figsize'] = (14, 8)

    # Top-10 from gold standard defines union set of taxa
    gold_top10 = get_top10(gold_df)['taxon'].tolist()

    # Build combined table
    combined = pd.DataFrame({'taxon': gold_top10})

    for fn, df in test_data.items():
        sub = df[df['taxon'].isin(gold_top10)][['taxon', 'abundance']]
        combined[fn] = (
            sub.set_index('taxon')['abundance']
            .reindex(gold_top10)
            .fillna(0)
            .values
        )

    # Ranks
    ranks = combined.set_index('taxon').rank(ascending=False, method='min')

    # Correlations
    corr = combined.set_index('taxon').corr()

    # Bray–Curtis
    all_samples = {'gold': gold_df}
    all_samples.update(test_data)
    sample_names = list(all_samples.keys())

    bc = pd.DataFrame(
        np.zeros((len(sample_names), len(sample_names))),
        index=sample_names,
        columns=sample_names
    )

    for a in sample_names:
        for b in sample_names:
            v1 = all_samples[a].set_index('taxon')['abundance'].reindex(gold_top10, fill_value=0)
            v2 = all_samples[b].set_index('taxon')['abundance'].reindex(gold_top10, fill_value=0)
            bc.loc[a, b] = np.sum(np.abs(v1 - v2)) / np.sum(v1 + v2)

    # Write PDF
    with PdfPages(pdf_path) as pdf:

        # Page 1 – stacked barplot
        fig, ax = plt.subplots(figsize=(14, 8))
        melted = combined.melt(id_vars='taxon', var_name='sample', value_name='abundance')
        sns.barplot(data=melted, x='sample', y='abundance', hue='taxon', ax=ax)
        plt.xticks(rotation=90)
        ax.set_title("Top 10 Taxa – Normalised Abundances")
        ax.legend(title="Taxon", bbox_to_anchor=(1.02, 1), loc='upper left')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Page 2 – abundance heatmap
        fig, ax = plt.subplots(figsize=(12, 8))
        sns.heatmap(combined.set_index('taxon'), cmap='viridis')
        ax.set_title("Top 10 Taxa – Heatmap of Abundance")
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Page 3 – rank-shift
        fig, ax = plt.subplots(figsize=(14, 8))
        long_rank = ranks.reset_index().melt(id_vars='taxon',
                                             var_name='sample',
                                             value_name='rank')
        sns.lineplot(data=long_rank, x='sample', y='rank', hue='taxon', marker='o', ax=ax)
        plt.gca().invert_yaxis()
        plt.xticks(rotation=90)
        ax.set_title("Rank Shift of Top 10 Taxa Across Parameter Sets")
        ax.set_ylabel("Rank (1 = highest abundance)")
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Page 4 – correlation heatmap
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(corr, annot=True, cmap='coolwarm', vmin=-1, vmax=1)
        ax.set_title("Correlation of Top-10 Taxon Profiles")
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Page 5 – Bray–Curtis
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(bc, annot=True, cmap='magma', vmin=0, vmax=1)
        ax.set_title("Bray–Curtis Dissimilarity Between Samples")
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    print(f"\nSaved: {pdf_path}")


# ----------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compare top 10 taxa using the load_and_normalize() function."
    )

    parser.add_argument("-g", "--gold-standard", required=True,
                        help="Gold standard CSV for normalised taxa.")
    parser.add_argument("-d", "--directory", required=True,
                        help="Directory of test CSVs.")
    parser.add_argument("-p", "--pattern", default="*.csv",
                        help="File pattern (default: *.csv).")
    parser.add_argument("--feature-col", default="feature-id",
                        help="Column containing feature/taxon names.")
    parser.add_argument("--abundance-col", default=None,
                        help="Column containing abundance/count values. "
                             "If not given, the script uses the first numeric column.")

    parser.add_argument("-o", "--output", default="top10_comparison",
                        help="Output prefix")

    args = parser.parse_args()

    # Load gold standard
    print(f"Loading gold standard: {args.gold_standard}")
    gold_raw = load_and_normalize(args.gold_standard, args.feature_col)

    # Determine abundance column
    if args.abundance_col is None:
        possible = [c for c in gold_raw.columns if c != args.feature_col]
        if not possible:
            print("ERROR: No numeric columns available.", file=sys.stderr)
            sys.exit(1)
        abundance_col = possible[0]
        print(f"  Auto-selected abundance column: {abundance_col}")
    else:
        abundance_col = args.abundance_col

    gold = extract_taxa(gold_raw, args.feature_col, abundance_col)

    # Load test files
    d = Path(args.directory)
    files = list(d.glob(args.pattern))

    if not files:
        print("No matching CSV files found.")
        sys.exit(1)

    test_data = {}
    for f in files:
        print(f"Loading: {f.name}")
        test_raw = load_and_normalize(f, args.feature_col)
        test_df = extract_taxa(test_raw, args.feature_col, abundance_col)
        test_data[f.name] = test_df

    # Plot
    create_plots(gold, test_data, args.output)


if __name__ == "__main__":
    main()
