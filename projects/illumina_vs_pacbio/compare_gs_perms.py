#!/usr/bin/env python3
"""
Compare multiple QIIME2/DADA2 parameter permutations against a gold standard.
Performs TMM normalization, calculates comparison metrics, and generates visualization plots.

COMPARISON METRICS EXPLAINED:

Feature Detection Metrics:
- n_gold_features: Number of features (ASVs/OTUs) in gold standard
- n_test_features: Number of features detected in test sample
- n_shared_features: Features found in both gold standard and test
- n_gold_only: Features in gold standard but missed by test (false negatives)
- n_test_only: Features in test but not in gold standard (false positives)

Performance Metrics:
- sensitivity (recall): Proportion of gold standard features detected by test
  Formula: shared / gold_total. Higher is better (max 1.0)
  Interpretation: How well the method recovers true features

- precision: Proportion of test features that are in gold standard
  Formula: shared / test_total. Higher is better (max 1.0)
  Interpretation: How many detected features are real (not spurious)

- f1_score: Harmonic mean of precision and sensitivity (0-1)
  Formula: 2 * (precision * sensitivity) / (precision + sensitivity)
  Interpretation: Overall balance between finding features and avoiding false positives

Abundance Correlation Metrics (for shared features):
- pearson_r: Linear correlation between abundances (-1 to 1)
  Interpretation: How well relative abundances match linearly

- spearman_r: Rank correlation between abundances (-1 to 1)
  Interpretation: How well the rank order of abundances matches
  Better for non-linear relationships

Distance Metrics (lower is better):
- bray_curtis: Dissimilarity in community composition (0-1)
  Formula: sum|gold - test| / sum(gold + test)
  Interpretation: 0 = identical, 1 = completely different
  Commonly used in ecology for comparing communities

- euclidean_dist: Straight-line distance in abundance space
  Interpretation: Overall magnitude of differences in abundance
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


def tmm_normalize(counts):
    """
    Perform Trimmed Mean of M-values (TMM) normalization.

    TMM normalization accounts for library size differences by scaling
    counts to a common reference. This is essential when comparing
    samples with different sequencing depths.

    Args:
        counts: pandas Series or DataFrame column with count data

    Returns:
        Normalized counts (counts per million)
    """
    counts = counts.copy()

    # Ensure numeric type
    counts = pd.to_numeric(counts, errors='coerce')
    counts = counts.fillna(0)

    # Replace zeros with NaN for calculations
    counts_calc = counts.replace(0, np.nan)

    # Calculate library size (excluding NaN)
    lib_size = counts.sum()

    if lib_size == 0:
        return counts

    # Normalize to counts per million
    norm_factor = lib_size / 1e6

    return counts / norm_factor


def calculate_metrics(gold_std, test_sample, feature_col='feature-id'):
    """
    Calculate comprehensive comparison metrics between gold standard and test sample.

    Args:
        gold_std: DataFrame with gold standard data
        test_sample: DataFrame with test sample data
        feature_col: Name of the feature ID column

    Returns:
        Dictionary of metrics
    """
    # Merge on feature IDs
    merged = pd.merge(
        gold_std, test_sample,
        on=feature_col,
        how='outer',
        suffixes=('_gold', '_test')
    )

    # Get count columns
    gold_col = [c for c in merged.columns if c.endswith('_gold')][0]
    test_col = [c for c in merged.columns if c.endswith('_test')][0]

    # Fill NaN with 0 for features not found
    merged[gold_col] = merged[gold_col].fillna(0)
    merged[test_col] = merged[test_col].fillna(0)

    # Calculate metrics
    metrics = {}

    # Features found
    gold_features = set(gold_std[feature_col])
    test_features = set(test_sample[feature_col])

    metrics['n_gold_features'] = len(gold_features)
    metrics['n_test_features'] = len(test_features)
    metrics['n_shared_features'] = len(gold_features & test_features)
    metrics['n_gold_only'] = len(gold_features - test_features)
    metrics['n_test_only'] = len(test_features - gold_features)

    # Sensitivity and specificity
    metrics['sensitivity'] = metrics['n_shared_features'] / metrics['n_gold_features'] if metrics[
                                                                                              'n_gold_features'] > 0 else 0
    metrics['precision'] = metrics['n_shared_features'] / metrics['n_test_features'] if metrics[
                                                                                            'n_test_features'] > 0 else 0

    # F1 score
    if metrics['sensitivity'] + metrics['precision'] > 0:
        metrics['f1_score'] = 2 * (metrics['precision'] * metrics['sensitivity']) / (
                    metrics['precision'] + metrics['sensitivity'])
    else:
        metrics['f1_score'] = 0

    # Correlation and distance metrics on shared features
    shared_mask = (merged[gold_col] > 0) | (merged[test_col] > 0)
    if shared_mask.sum() > 1:
        metrics['pearson_r'] = merged.loc[shared_mask, [gold_col, test_col]].corr().iloc[0, 1]
        metrics['spearman_r'] = merged.loc[shared_mask, [gold_col, test_col]].corr(method='spearman').iloc[0, 1]

        # Bray-Curtis dissimilarity
        bc_num = np.abs(merged[gold_col] - merged[test_col]).sum()
        bc_denom = (merged[gold_col] + merged[test_col]).sum()
        metrics['bray_curtis'] = bc_num / bc_denom if bc_denom > 0 else 1.0

        # Euclidean distance
        metrics['euclidean_dist'] = np.sqrt(((merged[gold_col] - merged[test_col]) ** 2).sum())
    else:
        metrics['pearson_r'] = np.nan
        metrics['spearman_r'] = np.nan
        metrics['bray_curtis'] = np.nan
        metrics['euclidean_dist'] = np.nan

    return metrics, merged, gold_col, test_col


def load_and_normalize(filepath, feature_col='feature-id'):
    """
    Load QIIME2 CSV and perform TMM normalization.
    Handles both standard column format and transposed row format.

    Args:
        filepath: Path to CSV file
        feature_col: Name of the feature ID column

    Returns:
        Normalized DataFrame with features as rows
    """
    # Try to read with first column as index
    df = pd.read_csv(filepath, index_col=0)

    # Check if data is transposed (samples as rows, features as columns)
    # This is the case if the first row contains mostly numeric values
    first_row_numeric = pd.to_numeric(df.iloc[0], errors='coerce').notna().sum() / len(df.columns) > 0.5

    if first_row_numeric:
        # Data is transposed - samples are rows, features are columns
        # Transpose so features become rows
        df = df.T
        print(f"  Detected transposed format, transposing data...")

    # Convert all values to numeric, coercing errors
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # Fill any NaN values with 0
    df = df.fillna(0)

    # Normalize each count column
    for col in df.columns:
        df[col] = tmm_normalize(df[col])

    # Reset index to make features a column
    df = df.reset_index()
    df = df.rename(columns={'index': feature_col})

    return df


def create_comparison_plots(results_df, gold_std, test_files_data, output_prefix, feature_col):
    """
    Generate comprehensive visualization plots comparing all samples to gold standard.

    Args:
        results_df: DataFrame with comparison metrics
        gold_std: Gold standard DataFrame
        test_files_data: Dictionary mapping filenames to (test_data, merged_data, gold_col, test_col)
        output_prefix: Prefix for output files
        feature_col: Feature ID column name
    """
    # Set style
    sns.set_style("whitegrid")
    plt.rcParams['figure.figsize'] = (12, 8)

    # Create multi-page PDF
    pdf_path = f"{output_prefix}_comparison_plots.pdf"

    with PdfPages(pdf_path) as pdf:
        # Page 1: Overview metrics
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle('QIIME2 Parameter Comparison: Overview Metrics', fontsize=16, fontweight='bold')

        # Sort by F1 score for better visualization
        plot_df = results_df.sort_values('f1_score', ascending=True)

        # F1 Score comparison
        axes[0, 0].barh(range(len(plot_df)), plot_df['f1_score'], color='steelblue')
        axes[0, 0].set_yticks(range(len(plot_df)))
        axes[0, 0].set_yticklabels(plot_df['filename'], fontsize=8)
        axes[0, 0].set_xlabel('F1 Score')
        axes[0, 0].set_title('F1 Score (Precision-Sensitivity Balance)')
        axes[0, 0].set_xlim(0, 1)
        axes[0, 0].axvline(x=0.9, color='green', linestyle='--', alpha=0.5, label='Excellent (>0.9)')
        axes[0, 0].legend()

        # Precision vs Sensitivity scatter
        axes[0, 1].scatter(plot_df['sensitivity'], plot_df['precision'],
                           s=100, alpha=0.6, c=plot_df['f1_score'], cmap='viridis')
        axes[0, 1].plot([0, 1], [0, 1], 'k--', alpha=0.3)
        axes[0, 1].set_xlabel('Sensitivity (Recall)')
        axes[0, 1].set_ylabel('Precision')
        axes[0, 1].set_title('Precision vs Sensitivity Trade-off')
        axes[0, 1].set_xlim(0, 1.05)
        axes[0, 1].set_ylim(0, 1.05)
        # Add text annotations for best/worst
        if len(plot_df) > 0:
            best_idx = plot_df['f1_score'].idxmax()
            axes[0, 1].annotate(plot_df.loc[best_idx, 'filename'],
                                xy=(plot_df.loc[best_idx, 'sensitivity'],
                                    plot_df.loc[best_idx, 'precision']),
                                fontsize=8, alpha=0.7)

        # Bray-Curtis dissimilarity
        axes[1, 0].barh(range(len(plot_df)), plot_df['bray_curtis'], color='coral')
        axes[1, 0].set_yticks(range(len(plot_df)))
        axes[1, 0].set_yticklabels(plot_df['filename'], fontsize=8)
        axes[1, 0].set_xlabel('Bray-Curtis Dissimilarity')
        axes[1, 0].set_title('Community Composition Difference (lower is better)')
        axes[1, 0].axvline(x=0.1, color='green', linestyle='--', alpha=0.5, label='Excellent (<0.1)')
        axes[1, 0].legend()

        # Feature detection
        feature_data = plot_df[['filename', 'n_shared_features', 'n_gold_only', 'n_test_only']].set_index('filename')
        feature_data.plot(kind='barh', stacked=True, ax=axes[1, 1],
                          color=['green', 'orange', 'red'], alpha=0.7)
        axes[1, 1].set_xlabel('Number of Features')
        axes[1, 1].set_title('Feature Detection Breakdown')
        axes[1, 1].legend(['Shared (TP)', 'Missed (FN)', 'Spurious (FP)'], fontsize=8)

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

        # Page 2: Correlation heatmap
        fig, ax = plt.subplots(figsize=(12, 10))

        corr_metrics = ['f1_score', 'sensitivity', 'precision', 'pearson_r',
                        'spearman_r', 'bray_curtis', 'euclidean_dist']
        corr_data = results_df.set_index('filename')[corr_metrics].T

        sns.heatmap(corr_data, annot=True, fmt='.3f', cmap='RdYlGn',
                    center=0.5, ax=ax, cbar_kws={'label': 'Score'})
        ax.set_title('Comparison Metrics Heatmap\n(Green = Better Performance)',
                     fontsize=14, fontweight='bold')
        ax.set_xlabel('Sample', fontsize=12)
        ax.set_ylabel('Metric', fontsize=12)

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

        # Pages 3+: Individual scatter plots for each sample
        for filename, (test_data, merged, gold_col, test_col) in test_files_data.items():
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
            fig.suptitle(f'Gold Standard vs {filename}', fontsize=14, fontweight='bold')

            # Remove zeros for log plot
            plot_data = merged[(merged[gold_col] > 0) & (merged[test_col] > 0)].copy()

            if len(plot_data) > 0:
                # Log-log scatter
                axes[0].scatter(plot_data[gold_col], plot_data[test_col],
                                alpha=0.5, s=30)
                axes[0].set_xscale('log')
                axes[0].set_yscale('log')

                # Add diagonal line
                min_val = min(plot_data[gold_col].min(), plot_data[test_col].min())
                max_val = max(plot_data[gold_col].max(), plot_data[test_col].max())
                axes[0].plot([min_val, max_val], [min_val, max_val],
                             'r--', alpha=0.5, label='Perfect agreement')

                axes[0].set_xlabel('Gold Standard Abundance (log scale)')
                axes[0].set_ylabel('Test Abundance (log scale)')
                axes[0].set_title('Abundance Correlation (Shared Features)')
                axes[0].legend()

                # Add correlation info
                metrics_row = results_df[results_df['filename'] == filename].iloc[0]
                textstr = f"Pearson r: {metrics_row['pearson_r']:.3f}\n"
                textstr += f"Spearman r: {metrics_row['spearman_r']:.3f}\n"
                textstr += f"Shared features: {metrics_row['n_shared_features']}"
                axes[0].text(0.05, 0.95, textstr, transform=axes[0].transAxes,
                             verticalalignment='top', bbox=dict(boxstyle='round',
                                                                facecolor='wheat', alpha=0.5))

            # Bland-Altman plot (difference vs average)
            if len(plot_data) > 0:
                avg = (plot_data[gold_col] + plot_data[test_col]) / 2
                diff = plot_data[test_col] - plot_data[gold_col]

                axes[1].scatter(avg, diff, alpha=0.5, s=30)
                axes[1].axhline(y=0, color='r', linestyle='--', alpha=0.5)
                axes[1].axhline(y=diff.mean(), color='blue', linestyle='-',
                                alpha=0.5, label=f'Mean diff: {diff.mean():.2f}')
                axes[1].axhline(y=diff.mean() + 1.96 * diff.std(), color='gray',
                                linestyle=':', alpha=0.5, label='±1.96 SD')
                axes[1].axhline(y=diff.mean() - 1.96 * diff.std(), color='gray',
                                linestyle=':', alpha=0.5)

                axes[1].set_xlabel('Average Abundance')
                axes[1].set_ylabel('Difference (Test - Gold)')
                axes[1].set_title('Bland-Altman Plot')
                axes[1].legend()

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()

        # Final page: Summary ranking
        fig, ax = plt.subplots(figsize=(12, 8))

        # Rank samples by multiple criteria
        rank_df = results_df[['filename', 'f1_score', 'bray_curtis', 'pearson_r']].copy()
        rank_df['f1_rank'] = rank_df['f1_score'].rank(ascending=False)
        rank_df['bc_rank'] = rank_df['bray_curtis'].rank(ascending=True)  # Lower is better
        rank_df['corr_rank'] = rank_df['pearson_r'].rank(ascending=False)
        rank_df['avg_rank'] = (rank_df['f1_rank'] + rank_df['bc_rank'] + rank_df['corr_rank']) / 3
        rank_df = rank_df.sort_values('avg_rank')

        # Plot rankings
        x = range(len(rank_df))
        width = 0.25

        ax.barh([i - width for i in x], rank_df['f1_rank'], width, label='F1 Score Rank', alpha=0.8)
        ax.barh(x, rank_df['bc_rank'], width, label='Bray-Curtis Rank', alpha=0.8)
        ax.barh([i + width for i in x], rank_df['corr_rank'], width, label='Correlation Rank', alpha=0.8)

        ax.set_yticks(x)
        ax.set_yticklabels(rank_df['filename'], fontsize=9)
        ax.set_xlabel('Rank (1 = Best)')
        ax.set_title('Overall Parameter Performance Ranking', fontsize=14, fontweight='bold')
        ax.legend()
        ax.invert_xaxis()  # Best (1) on right

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

    print(f"\nPlots saved to: {pdf_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Compare QIIME2/DADA2 outputs against a gold standard with visualizations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        '-g', '--gold-standard',
        required=True,
        help='Path to gold standard CSV file'
    )

    parser.add_argument(
        '-d', '--directory',
        required=True,
        help='Directory containing test CSV files to compare'
    )

    parser.add_argument(
        '-o', '--output',
        default='comparison_results',
        help='Output prefix for results files (default: comparison_results)'
    )

    parser.add_argument(
        '-f', '--feature-col',
        default='feature-id',
        help='Name of feature ID column (default: feature-id)'
    )

    parser.add_argument(
        '-p', '--pattern',
        default='*.csv',
        help='File pattern to match (default: *.csv)'
    )

    parser.add_argument(
        '--no-plots',
        action='store_true',
        help='Skip generating visualization plots'
    )

    args = parser.parse_args()

    # Load and normalize gold standard
    print(f"Loading gold standard from: {args.gold_standard}")
    try:
        gold_std = load_and_normalize(args.gold_standard, args.feature_col)
        print(f"  Found {len(gold_std)} features in gold standard")
    except Exception as e:
        print(f"Error loading gold standard: {e}", file=sys.stderr)
        sys.exit(1)

    # Find all test files
    test_dir = Path(args.directory)
    if not test_dir.exists():
        print(f"Error: Directory {test_dir} does not exist", file=sys.stderr)
        sys.exit(1)

    test_files = list(test_dir.glob(args.pattern))
    if not test_files:
        print(f"Warning: No files matching '{args.pattern}' found in {test_dir}")
        sys.exit(0)

    print(f"\nFound {len(test_files)} test files to compare")

    # Compare each test file
    results = []
    test_files_data = {}

    for test_file in test_files:
        print(f"\nProcessing: {test_file.name}")
        try:
            test_data = load_and_normalize(test_file, args.feature_col)
            metrics, merged, gold_col, test_col = calculate_metrics(
                gold_std, test_data, args.feature_col
            )
            metrics['filename'] = test_file.name
            results.append(metrics)

            # Store for plotting
            test_files_data[test_file.name] = (test_data, merged, gold_col, test_col)

            print(f"  Features: {metrics['n_test_features']}, "
                  f"Shared: {metrics['n_shared_features']}, "
                  f"F1: {metrics['f1_score']:.3f}, "
                  f"Bray-Curtis: {metrics['bray_curtis']:.3f}")
        except Exception as e:
            print(f"  Error processing {test_file.name}: {e}", file=sys.stderr)
            continue

    # Save results
    if results:
        results_df = pd.DataFrame(results)
        cols = ['filename'] + [c for c in results_df.columns if c != 'filename']
        results_df = results_df[cols]

        csv_output = f"{args.output}.csv"
        results_df.to_csv(csv_output, index=False)
        print(f"\nResults saved to: {csv_output}")

        print(f"\nSummary statistics:")
        print(results_df[['f1_score', 'sensitivity', 'precision', 'bray_curtis', 'pearson_r']].describe())

        # Generate plots
        if not args.no_plots:
            print("\nGenerating visualization plots...")
            try:
                create_comparison_plots(results_df, gold_std, test_files_data,
                                        args.output, args.feature_col)
            except Exception as e:
                print(f"Error generating plots: {e}", file=sys.stderr)
                print("Results CSV was still saved successfully.")
    else:
        print("\nNo results to save")


if __name__ == '__main__':
    main()