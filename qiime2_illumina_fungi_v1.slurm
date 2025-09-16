#!/bin/sh
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --job-name=qiime_&&
#SBATCH --error=/users/&&/jobs/Qiime2_&&.err
#SBATCH --output=/users/&&/jobs/Qiime2_&&.txt
#SBATCH --partition=&& # Replace with system-specific partition.
#SBATCH --time=60:00:00
#SBATCH --nodes=1
#SBATCH --mail-user=$$
#SBATCH --mail-type=BEGIN,END,FAIL

# If you have a high number of samples, consider increasing --mem and --time.


num_threads=10 # This is the number of CPU threads to use.

#Directory containing project data.
cd .../

## Make sure fastq files are in a directory called 'data'

# Make a directory 'qiime' for all output files
mkdir -p qiime

# Initialise the output file
echo -e "sampleID\tforward-absolute-filepath\treverse-absolute-filepath" > ./qiime/illumina-manifest.tsv

# Loop through each forward (R1) and reverse (R2) file
for forward in data/*_R1.fastq.gz; do
    # Extract the sampleID from the filename (assuming the sample ID is the part before _R1)
    sampleID=$(basename $forward _R1.fastq.gz)
    # Find the corresponding reverse (R2) file
    reverse="data/${sampleID}_R2.fastq.gz"
    # Check if the reverse file exists before adding to the manifest
    if [[ -f "$reverse" ]]; then
        echo -e "$sampleID\t$(realpath $forward)\t$(realpath $reverse)" >> ./qiime/illumina-manifest.tsv
    else
        echo "Warning: Reverse file for $sampleID not found. Skipping..." >&2
    fi
done


for forward in /data/*/*_R1_*.fastq.gz; do
    filename=$(basename "$forward")
    sampleID=$(echo "$filename" | sed 's/_R1_001\.fastq\.gz$//')
    forward_dir=$(dirname "$forward")
    reverse="${forward_dir}/${sampleID}_R2_001.fastq.gz"
    if [[ -f "$reverse" ]]; then
        echo -e "$sampleID\t$(realpath "$forward")\t$(realpath "$reverse")" >> ./qiime/illumina-manifest.tsv
    else
        echo "Warning: Reverse file for $sampleID not found. Skipping..." >&2
    fi
done

# Activate the version of QIIME2 we want to use - I suggest 'qiime2-amplicon-2024.10'
source activate .../conda-envs/qiime2-amplicon-2024.10

# import data

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./qiime/illumina-manifest.tsv \
  --output-path ./qiime/demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

# Validate demux file
qiime tools validate ./qiime/demux.qza

# export demux summary
qiime demux summarize --i-data ./qiime/demux.qza --o-visualization qiime/demux_summary

# denoise - set the trim and truncation lengths according to your data quality.
#removed primer sequences not by cut adapt but by trimming front

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./qiime/demux.qza \
  --p-trim-left-f 21 \
  --p-trim-left-r 22 \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 240 \
  --p-n-threads 6 \
  --o-representative-sequences ./qiime/rep-seqs.qza \
  --o-denoising-stats ./qiime/denoising-stats.qza \
  --o-table ./qiime/table.qza \
  --verbose > ./qiime/denoising_error.txt 2>&1


# export  table

qiime tools export --input-path ./qiime/table.qza --output-path ./qiime/table_export

# export denoising

qiime tools export --input-path ./qiime/denoising-stats.qza --output-path ./qiime/denoising-stats_export



# summary of table and rep-seqs
qiime feature-table summarize --i-table ./qiime/table.qza --o-visualization ./qiime/table.qzv

qiime feature-table tabulate-seqs --i-data ./qiime/rep-seqs.qza --o-visualization ./qiime/rep-seqs.qzv

# Optional filtering of the feature table to remove low-abundance features
##########################################################################
# Filter features that only appear in a single sample
qiime feature-table filter-features \
  --i-table ./qiime/table.qza \
  --p-min-samples 2 \
  --o-filtered-table ./qiime/table_min_2samples.qza

#Can remove all features with a total abundance (summed across all samples) of less than 10.
qiime feature-table filter-features \
  --i-table ./qiime/table_min_2samples.qza \
  --p-min-frequency 10 \
  --o-filtered-table ./qiime/table_min_2samples_mintotal10.qza

# Because we filtered some ASVs out, the rep-seq file has more sequences than dada2 table.
# We will now filter rep-seq file so that only IDs present in table will remain.
qiime feature-table filter-seqs \
 --i-data ./qiime/rep-seqs.qza \
 --i-table ./qiime/table_min_2samples_mintotal10.qza \
 --o-filtered-data ./qiime/rep-seqs_filtered.qza


#let's have a look what we retain
qiime feature-table tabulate-seqs \
      --i-data ./qiime/rep-seqs_filtered.qza \
      --o-visualization ./qiime/rep-seqs_filtered.qzv
##########################################################################


# classify - here we are using the 'Anaerobic Fungi Network' (AFN) classifier provided as part of this repo.
# Change 'Filtered' naming below if you did not filter your table/rep-seq files.

qiime feature-classifier classify-sklearn \
  --i-classifier ./datasets/sklearn_models/AGF_Fungi_Network/AGF_Full_region_v2_classifier.qza   \
  --i-reads ./qiime/rep-seqs_filtered.qza \
  --o-classification ./qiime/sklearn-taxonomy_AGF_Full_v2_Filtered_ASV.qza \
  --p-n-jobs $num_threads \
  --verbose

# metadata
qiime metadata tabulate \
  --m-input-file ./qiime/sklearn-taxonomy_Fungi_AFN_ASV.qza \
  --o-visualization ./qiime/sklearn-taxonomy_Fungi_AFN_ASV.qzv

# visualisation
qiime taxa barplot \
 --i-table ./qiime/table_min_2samples_mintotal10.qza \
 --i-taxonomy ./qiime/sklearn-taxonomy_AGF_Full_v2_Filtered_ASV.qza \
 --o-visualization ./qiime/sklearn-taxonomy_AGF_Full_v2_Filtered_ASV_Bar


# export representative sequences

qiime tools export \
  --input-path ./qiime/rep-seqs.qza \
    --output-path ./qiime/rep-seqs_export
