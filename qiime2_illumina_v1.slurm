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

# Activate the version of QIIME2 we want to use - I suggest 'qiime2-amplicon-2024.10'
source activate .../conda-envs/qiime2-amplicon-2024.10  

# import data

$COMMAND qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./qiime/illumina-manifest.tsv \
  --output-path ./qiime/demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

# denoise

$COMMAND qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ./qiime/demux.qza \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --p-trunc-len-f 245 \
  --p-trunc-len-r 245 \
  --p-n-threads $num_threads \
  --o-representative-sequences ./qiime/rep-seqs.qza \
  --o-denoising-stats ./qiime/denoising-stats.qza \
  --o-table ./qiime/table.qza \
  --verbose > ./qiime/denoising_error.txt 2>&1


# export  table

$COMMAND qiime tools export --input-path ./qiime/table.qza --output-path ./qiime/table_export

# export denoising

$COMMAND qiime tools export --input-path ./qiime/denoising-stats.qza --output-path ./qiime/denoising-stats_export

# export summary

$COMMAND qiime demux summarize --i-data qiime/demux.qza --o-visualization qiime/demux_summary



# classify - Download sklearn model from https://library.qiime2.org/data-resources#naive-bayes-classifiers

qiime feature-classifier classify-sklearn \
  --i-classifier .../2024.09.backbone.v4.nb.sklearn-1.4.2.qza  \
  --i-reads ./qiime/rep-seqs.qza \
  --o-classification ./qiime/sklearn-taxonomy_GG2_ASV.qza \
  --p-n-jobs $num_threads \
  --verbose

# metadata
$COMMMAND qiime metadata tabulate \
  --m-input-file ./qiime/sklearn-taxonomy_GG2_ASV.qza \
  --o-visualization ./qiime/sklearn-taxonomy_GG2_ASV.qzv

# visualisation
$COMMAND qiime taxa barplot \
 --i-table ./qiime/table.qza \
 --i-taxonomy ./qiime/sklearn-taxonomy_GG2_ASV.qza \
 --o-visualization ./qiime/sklearn-taxonomy_GG2_ASV_Bar


# export representative sequences

$COMMAND qiime tools export \
  --input-path ./qiime/rep-seqs.qza \
    --output-path ./qiime/rep-seqs_export
