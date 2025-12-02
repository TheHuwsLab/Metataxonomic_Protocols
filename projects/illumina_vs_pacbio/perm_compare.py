#!/usr/bin/env python3
"""
Compare top 10 most abundant taxa across multiple QIIME2 parameter sets.

This script loads multiple level-6.csv files (different QIIME2 parameters),
normalizes the data, identifies the top 10 most abundant taxa in each,
and visualizes whether these top taxa are consistent across parameter sets.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import Counter


def load_and_normalize(filepath):
    """
    Load QIIME2 level-6 CSV and normalize to relative abundance.

    Args:
        filepath: Path to level-6.csv file

    Returns:
        DataFrame with normalized abundances (samples as columns, taxa as rows)
    """
    # Read CSV with first column as index (sample IDs)
    df = pd.read_csv(filepath, index_col=0)

    # Transpose so taxa are rows and samples are columns
    df = df.T

    # Ensure all values are numeric
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df = df.fillna(0)

    # Normalize each sample to relative abundance (sum to 1)
    df_norm = df.div(df.sum(axis=0), axis=1)

    return df_norm


def get_top_n_taxa(df, n=10, aggregate='mean'):
    """
    Get top N most abundant taxa across all samples.

    Args:
        df: DataFrame with taxa as rows, samples as columns
        n: Number of top taxa to return
        aggregate: How to aggregate across samples ('mean', 'median', 'max')

    Returns:
        List of top N taxa names
    """
    if aggregate == 'mean':
        abundance = df.mean(axis=1)
    elif aggregate == 'median':
        abundance = df.median(axis=1)
    elif aggregate == 'max':
        abundance = df.max(axis=1)
    else:
        raise ValueError(f"Unknown aggregate method: {aggregate}")

    top_taxa = abundance.nlargest(n).index.tolist()
    return top_taxa


def shorten_taxa_name(taxa_name, max_length=60):
    """
    Shorten long taxonomic names for better visualization.

    Args:
        taxa_name: Full taxonomic string
        max_length: Maximum length for display

    Returns:
        Shortened taxa name
    """
    if len(taxa_name) <= max_length:
        return taxa_name

    # Try to extract genus level (last level before __)
    parts = taxa_name.split(';')
    for part in reversed(parts):
        if part and not part.endswith('__'):
            genus = part.split('__')[-1] if '__' in part else part
            return f"...;{genus}"

    return taxa_name[:max_length - 3] + "..."


def create_comparison_plots(data_dict, output_prefix, top_n=10, aggregate='mean'):
    """
    Create comprehensive visualization of top taxa across parameter sets.

    Args:
        data_dict: Dictionary mapping filename to normalized DataFrame
        output_prefix: Prefix for output files
        top_n: Number of top taxa to analyze
        aggregate: How to aggregate abundances across samples
    """
    # Collect top taxa from each parameter set
    all_top_taxa = {}
    taxa_counter = Counter()

    print(f"\nTop {top_n} taxa for each parameter set:")
    print("=" * 80)

    for filename, df in data_dict.items():
        top_taxa = get_top_n_taxa(df, n=top_n, aggregate=aggregate)
        all_top_taxa[filename] = top_taxa
        taxa_counter.update(top_taxa)

        print(f"\n{filename}:")
        for i, taxa in enumerate(top_taxa, 1):
            avg_abundance = df.loc[taxa].mean()
            print(f"  {i:2d}. {shorten_taxa_name(taxa)}")
            print(f"      Mean abundance: {avg_abundance:.4f}")

    # Get union of all top taxa across all parameter sets
    all_unique_taxa = list(taxa_counter.keys())

    print(f"\n\nTotal unique taxa in top {top_n} across all parameter sets: {len(all_unique_taxa)}")
    print(
        f"Taxa appearing in all {len(data_dict)} parameter sets: {sum(1 for count in taxa_counter.values() if count == len(data_dict))}")

    # Set style
    sns.set_style("whitegrid")

    # Figure 1: Heatmap of taxa presence in top 10
    fig, ax = plt.subplots(figsize=(12, max(8, len(all_unique_taxa) * 0.3)))

    # Create presence/absence matrix
    presence_matrix = pd.DataFrame(0,
                                   index=all_unique_taxa,
                                   columns=list(data_dict.keys()))

    for filename, top_taxa in all_top_taxa.items():
        for rank, taxa in enumerate(top_taxa, 1):
            presence_matrix.loc[taxa, filename] = rank

    # Sort by frequency of appearance
    taxa_frequency = (presence_matrix > 0).sum(axis=1).sort_values(ascending=False)
    presence_matrix = presence_matrix.loc[taxa_frequency.index]

    # Shorten taxa names for display
    display_names = [shorten_taxa_name(taxa) for taxa in presence_matrix.index]

    # Create heatmap
    sns.heatmap(presence_matrix.values,
                xticklabels=presence_matrix.columns,
                yticklabels=display_names,
                cmap='YlOrRd',
                cbar_kws={'label': 'Rank (0 = not in top 10)'},
                ax=ax,
                linewidths=0.5,
                vmin=0, vmax=top_n)

    ax.set_title(f'Top {top_n} Taxa Ranking Across Parameter Sets\n(Lower number = higher rank)',
                 fontsize=14, fontweight='bold', pad=20)
    ax.set_xlabel('Parameter Set', fontsize=12)
    ax.set_ylabel('Taxon', fontsize=12)

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_top{top_n}_heatmap.png", dpi=300, bbox_inches='tight')
    print(f"\nSaved heatmap to: {output_prefix}_top{top_n}_heatmap.png")
    plt.close()

    # Figure 2: Line plot of abundances for most consistent taxa
    # Get taxa that appear in at least 50% of parameter sets
    consistent_taxa = [taxa for taxa, count in taxa_counter.items()
                       if count >= len(data_dict) * 0.5]

    if len(consistent_taxa) > 0:
        # Limit to top 10 most frequent
        consistent_taxa = consistent_taxa[:min(10, len(consistent_taxa))]

        fig, ax = plt.subplots(figsize=(14, 8))

        # Prepare data for line plot
        x_positions = range(len(data_dict))
        colors = plt.cm.tab20(np.linspace(0, 1, len(consistent_taxa)))

        for i, taxa in enumerate(consistent_taxa):
            abundances = []
            for filename in data_dict.keys():
                df = data_dict[filename]
                if taxa in df.index:
                    abundances.append(df.loc[taxa].mean())
                else:
                    abundances.append(0)

            ax.plot(x_positions, abundances,
                    marker='o', linewidth=2, markersize=8,
                    label=shorten_taxa_name(taxa, max_length=50),
                    color=colors[i], alpha=0.7)

        ax.set_xticks(x_positions)
        ax.set_xticklabels(list(data_dict.keys()), rotation=45, ha='right')
        ax.set_xlabel('Parameter Set', fontsize=12, fontweight='bold')
        ax.set_ylabel('Mean Relative Abundance', fontsize=12, fontweight='bold')
        ax.set_title(
            f'Abundance of Most Consistent Taxa Across Parameter Sets\n(Taxa appearing in ≥50% of parameter sets)',
            fontsize=14, fontweight='bold')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(f"{output_prefix}_abundance_lines.png", dpi=300, bbox_inches='tight')
        print(f"Saved line plot to: {output_prefix}_abundance_lines.png")
        plt.close()
    else:
        print("\nNo consistent taxa found across parameter sets for line plot.")

    # Figure 3: Bar chart showing taxa frequency
    fig, ax = plt.subplots(figsize=(12, 6))

    taxa_freq_sorted = taxa_counter.most_common()
    taxa_names = [shorten_taxa_name(t[0]) for t in taxa_freq_sorted]
    frequencies = [t[1] for t in taxa_freq_sorted]

    colors_bar = ['darkgreen' if f == len(data_dict) else 'steelblue' if f >= len(data_dict) * 0.5 else 'lightcoral'
                  for f in frequencies]

    bars = ax.barh(range(len(taxa_names)), frequencies, color=colors_bar, alpha=0.7)
    ax.set_yticks(range(len(taxa_names)))
    ax.set_yticklabels(taxa_names, fontsize=9)
    ax.set_xlabel('Number of Parameter Sets Where Taxon is in Top 10', fontsize=12)
    ax.set_title(f'Consistency of Top {top_n} Taxa Across Parameter Sets',
                 fontsize=14, fontweight='bold')
    ax.axvline(x=len(data_dict), color='green', linestyle='--', alpha=0.5,
               label=f'Present in all {len(data_dict)} sets')
    ax.axvline(x=len(data_dict) * 0.5, color='orange', linestyle='--', alpha=0.5,
               label='Present in ≥50% of sets')
    ax.legend()
    ax.invert_yaxis()

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_frequency_bars.png", dpi=300, bbox_inches='tight')
    print(f"Saved frequency bar chart to: {output_prefix}_frequency_bars.png")
    plt.close()

    # Save summary table
    summary_df = pd.DataFrame({
        'Taxon': [t for t, _ in taxa_freq_sorted],
        'Appearances': [c for _, c in taxa_freq_sorted],
        'Percentage': [c / len(data_dict) * 100 for _, c in taxa_freq_sorted]
    })

    summary_file = f"{output_prefix}_top{top_n}_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"Saved summary table to: {summary_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Compare top 10 most abundant taxa across QIIME2 parameter sets',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '-d', '--directory',
        required=True,
        help='Directory containing level-6.csv files'
    )

    parser.add_argument(
        '-o', '--output',
        default='taxa_comparison',
        help='Output prefix for plots and tables (default: taxa_comparison)'
    )

    parser.add_argument(
        '-n', '--top-n',
        type=int,
        default=10,
        help='Number of top taxa to analyze (default: 10)'
    )

    parser.add_argument(
        '-p', '--pattern',
        default='*.csv',
        help='File pattern to match (default: *.csv)'
    )

    parser.add_argument(
        '-a', '--aggregate',
        choices=['mean', 'median', 'max'],
        default='mean',
        help='How to aggregate abundances across samples (default: mean)'
    )

    args = parser.parse_args()

    # Find all level-6 CSV files
    data_dir = Path(args.directory)
    if not data_dir.exists():
        print(f"Error: Directory {data_dir} does not exist")
        return 1

    csv_files = list(data_dir.glob(args.pattern))
    if not csv_files:
        print(f"Error: No files matching '{args.pattern}' found in {data_dir}")
        return 1

    print(f"Found {len(csv_files)} level-6.csv files to analyze")

    # Load and process all files
    data_dict = {}
    for csv_file in csv_files:
        print(f"\nLoading: {csv_file.name}")
        try:
            df_norm = load_and_normalize(csv_file)
            data_dict[csv_file.stem] = df_norm  # Use stem to remove .csv
            print(f"  Loaded {len(df_norm)} taxa across {len(df_norm.columns)} samples")
        except Exception as e:
            print(f"  Error loading {csv_file.name}: {e}")
            continue

    if not data_dict:
        print("\nNo data loaded successfully")
        return 1

    # Create comparison plots
    print("\n" + "=" * 80)
    print("GENERATING COMPARISON PLOTS")
    print("=" * 80)

    create_comparison_plots(data_dict, args.output, args.top_n, args.aggregate)

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)

    return 0


if __name__ == '__main__':
    exit(main())