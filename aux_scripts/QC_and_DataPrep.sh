#!/bin/sh
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --job-name=qiime_&&
#SBATCH --error=/users/&&/jobs/Qiime2_&&.err
#SBATCH --output=/users/&&/jobs/Qiime2_&&.txt
#SBATCH --partition=&& # Replace with system-specific partition.
#SBATCH --time=60:00:00 
#SBATCH --nodes=1
#SBATCH --mail-user=$$
#SBATCH --mail-type=BEGIN,END,FAIL

# If you have a high number of samples, consider increasing --mem and --time.


num_threads=20 # This is the number of CPU threads to use. 

#Directory containing project data.
cd .../

## Make sure fastq files are in a directory called 'data'

# Make a directory 'qiime' for all output files
mkdir -p qiime

# Initialise the output file
echo -e "sampleID\tforward-absolute-filepath\treverse-absolute-filepath" > ./qiime/illumina-manifest.tsv

############ Pick one of the following two options ############
############ depending on how your data is organised ##########
#Option 1:
# Use if samples are all in a single directory: Loop through each forward (R1) and reverse (R2) file

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

#Option 2:

# Use if samples are all in individual directories: Loop through each forward (R1) and reverse (R2) file
# for forward in /data/*/*_R1*.fastq.gz; do
#     filename=$(basename "$forward")
#     sampleID=$(echo "$filename" | sed 's/_R1\.fastq\.gz$//')
#     forward_dir=$(dirname "$forward")
#     reverse="${forward_dir}/${sampleID}_R2.fastq.gz"
#     if [[ -f "$reverse" ]]; then
#         echo -e "$sampleID\t$(realpath "$forward")\t$(realpath "$reverse")" >> ./qiime/illumina-manifest.tsv
#     else
#         echo "Warning: Reverse file for $sampleID not found. Skipping..." >&2
#     fi
# done
##############################################################

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



# denoise - set the trim and truncation lengths according to your data quality.

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

# export denoising

qiime tools export --input-path ./qiime/denoising-stats.qza --output-path ./qiime/denoising-stats_export

# export summary

qiime demux summarize --i-data qiime/demux.qza --o-visualization qiime/demux_summary


# export  table

qiime tools export --input-path ./qiime/table.qza --output-path ./qiime/table_export

biom convert -i qiime/table_export/feature-table.biom  -o qiime/table_export/feature-table.tsv --to-tsv --header-key taxonomy


# export representative sequences

qiime tools export \
  --input-path ./qiime/rep-seqs.qza \
    --output-path ./qiime/rep-seqs_export


# If you want to skip dada2 denoising 
#############
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ./qiime/illumina-manifest.tsv \
  --output-path ./qiime/demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime quality-filter q-score \
  --i-demux demux-seqs.qza \
  --o-filtered-sequences filtered-seqs.qza \
  --o-filter-stats filtered-stats.qza
#  --p-min-quality 20 \
#  --p-quality-window 3 \
#  --p-min-length-fraction 0.75


# Convert filtered sequences to the format needed for classification
qiime vsearch dereplicate-sequences \
  --i-sequences filtered-seqs.qza \
  --o-dereplicated-table derep-table.qza \
  --o-dereplicated-sequences derep-seqs.qza