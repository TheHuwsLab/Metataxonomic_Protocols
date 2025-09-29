#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_directory> -o <output_directory>"
    echo "  -i, --input     Input manifest file path where seqs are located"
    echo "  -t, --threads   Number of threads (default: 10)"
    echo "  -c, --classifier Classifier location"
    echo "  -o, --output    Output directory path"
    echo "  -h, --help      Display this help message"
    exit 1
}

# Initialise variables
MANIFEST=""
BASE_OUTPUT_DIR=""
THREADS=10
CLASSIFIER=""


# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            MANIFEST="$2"
            shift 2
            ;;
        -o|--output)
            BASE_OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -c|--classifier)
            CLASSIFIER="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Error: Unknown option $1"
            usage
            ;;
    esac
done

# Check if required arguments are provided
if [[ -z "$MANIFEST" || -z "$BASE_OUTPUT_DIR" ]]; then
    echo "Error: Both input and output paths must be specified."
    usage
fi

# Validate input directory exists
if [[ ! -f "$MANIFEST" ]]; then
    echo "Error: Manifest path '$MANIFEST' does not exist."
    exit 1
fi

# Create output directory if it doesn't exist
if [[ ! -d "$BASE_OUTPUT_DIR" ]]; then
    echo "Creating output directory: $BASE_OUTPUT_DIR"
    mkdir -p "$BASE_OUTPUT_DIR"
    if [[ $? -ne 0 ]]; then
        echo "Error: Failed to create output directory '$BASE_OUTPUT_DIR'"
        exit 1
    fi
fi

# Check if classifier argument is provided and exists
if [[ -z "$CLASSIFIER" ]]; then
    echo "Error: Classifier path must be specified."
    usage
fi
if [[ ! -f "$CLASSIFIER" ]]; then
    echo "Error: Classifier file '$CLASSIFIER' does not exist."
    exit 1
fi


# Convert to absolute paths
MANIFEST=$(realpath "$MANIFEST")
BASE_OUTPUT_DIR=$(realpath "$BASE_OUTPUT_DIR")


## Activate env
source activate qiime2-amplicon-2024.10

# Define parameter combinations
# Format: "trim_left_f,trim_left_r,trunc_len_f,trunc_len_r"
PARAMETER_COMBINATIONS=(
    "10,10,240,240"
    "10,10,250,250"
    "10,10,230,230"
    "15,15,240,240"
    "15,15,250,250"
    "20,20,240,240"
    "5,5,240,240"
    "10,15,240,230"
    "15,10,250,240"
)

# make output dir
mkdir -p $BASE_OUTPUT_DIR/qiime

# import data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $MANIFEST \
  --output-path $BASE_OUTPUT_DIR/qiime/demux.qza \
  --input-format PairedEndFastqManifestPhred33V2



echo "Starting QIIME DADA2 parameter testing..."
echo "Input file: $MANIFEST"
echo "Base output directory: $BASE_OUTPUT_DIR"
echo "Number of threads: $THREADS"
echo "Testing ${#PARAMETER_COMBINATIONS[@]} parameter combinations"
echo "----------------------------------------"

# Counter for tracking progress
counter=1
total=${#PARAMETER_COMBINATIONS[@]}

# Loop through parameter combinations
for params in "${PARAMETER_COMBINATIONS[@]}"; do
    # Parse parameters
    IFS=',' read -r trim_left_f trim_left_r trunc_len_f trunc_len_r <<< "$params"

    # Create output directory name based on parameters
    output_subdir="trim_${trim_left_f}_${trim_left_r}_trunc_${trunc_len_f}_${trunc_len_r}"
    full_output_dir="${BASE_OUTPUT_DIR}/${output_subdir}"

    echo "[$counter/$total] Processing parameters: trim-left-f=$trim_left_f, trim-left-r=$trim_left_r, trunc-len-f=$trunc_len_f, trunc-len-r=$trunc_len_r"
    echo "Output directory: $full_output_dir"

    # Create output directory
    mkdir -p "$full_output_dir"

    # Run QIIME DADA2 denoise-paired
    echo "Running QIIME DADA2..."
    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs "$BASE_OUTPUT_DIR/qiime/demux.qza" \
        --p-trim-left-f "$trim_left_f" \
        --p-trim-left-r "$trim_left_r" \
        --p-trunc-len-f "$trunc_len_f" \
        --p-trunc-len-r "$trunc_len_r" \
        --p-n-threads "$THREADS" \
        --o-representative-sequences "${full_output_dir}/rep-seqs.qza" \
        --o-denoising-stats "${full_output_dir}/denoising-stats.qza" \
        --o-table "${full_output_dir}/table.qza" \
        --verbose > "${full_output_dir}/denoising_error.txt" 2>&1

    # Check if command was successful
    if [[ $? -eq 0 ]]; then
        echo "✓ Successfully completed parameter set $counter"


        if command -v qiime &> /dev/null; then
            echo "Generating quick summary..."
            qiime feature-table summarize \
                --i-table "${full_output_dir}/table.qza" \
                --o-visualization "${full_output_dir}/table-summary.qzv" 2>/dev/null

            qiime feature-table tabulate-seqs \
                --i-data "${full_output_dir}/rep-seqs.qza" \
                --o-visualization "${full_output_dir}/rep-seqs-summary.qzv" 2>/dev/null

            qiime tools export \
                --input-path "${full_output_dir}/denoising-stats.qza" \
                --output-path "${full_output_dir}/denoising-stats_exported" 2>/dev/null

            qiime tools export \
                --input-path "${full_output_dir}/rep-seqs.qza"  \
                --output-path "${full_output_dir}/rep-seqs_export" 2>/dev/null

            qiime tools export \
                --input-path "${full_output_dir}/table.qza" \
                 --output-path "${full_output_dir}/table_export" 2>/dev/null
            biom convert \
                -i "${full_output_dir}/table_export/feature-table.biom" \
                -o "${full_output_dir}/table_export/feature-table.tsv" \
                --to-tsv --header-key taxonomy 2>/dev/null

            ### Classify
            qiime feature-classifier classify-sklearn \
            --i-classifier $CLASSIFIER  \
            --i-reads "${full_output_dir}/rep-seqs.qza" \
            --o-classification "${full_output_dir}/sklearn-taxonomy_GG2_V4_ASV.qza" \
            --p-n-jobs $THREADS \
            --verbose 2> "${full_output_dir}/classification_error.txt"

            qiime metadata tabulate \
                --m-input-file "${full_output_dir}/sklearn-taxonomy_GG2_V4_ASV.qza" \
                --o-visualization "${full_output_dir}/sklearn-taxonomy_GG2_V4_ASV.qzv" 2>/dev/null

            qiime taxa barplot \
                --i-table "${full_output_dir}/table.qza" \
                --i-taxonomy "${full_output_dir}/sklearn-taxonomy_GG2_V4_ASV.qza" \
                --o-visualization "${full_output_dir}/sklearn-taxonomy_GG2_V4_ASV_Bar" 2>/dev/null

        fi

    else
        echo "✗ Failed parameter set $counter - check ${full_output_dir}/denoising_error.txt"
    fi

    echo "----------------------------------------"
    ((counter++))
done

echo "Parameter testing completed!"
echo "Results saved in: $BASE_OUTPUT_DIR"
echo ""
echo "Summary of runs:"
ls -la "$BASE_OUTPUT_DIR"

# Optional: Create a summary report
summary_file="${BASE_OUTPUT_DIR}/parameter_summary.txt"
echo "QIIME DADA2 Parameter Testing Summary" > "$summary_file"
echo "Generated on: $(date)" >> "$summary_file"
echo "Input file: $BASE_OUTPUT_DIR/qiime/demux.qza" >> "$summary_file"
echo "Total combinations tested: $total" >> "$summary_file"
echo "" >> "$summary_file"
echo "Parameter combinations:" >> "$summary_file"

for params in "${PARAMETER_COMBINATIONS[@]}"; do
    IFS=',' read -r trim_left_f trim_left_r trunc_len_f trunc_len_r <<< "$params"
    output_subdir="trim_${trim_left_f}_${trim_left_r}_trunc_${trunc_len_f}_${trunc_len_r}"

    if [[ -f "${BASE_OUTPUT_DIR}/${output_subdir}/rep-seqs.qza" ]]; then
        status="SUCCESS"
    else
        status="FAILED"
    fi

    echo "$output_subdir: $status" >> "$summary_file"
done

echo "Summary report saved to: $summary_file"





























