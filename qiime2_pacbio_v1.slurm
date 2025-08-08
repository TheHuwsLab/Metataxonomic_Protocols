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


num_threads=24 # This is the number of CPU threads to use. 

#Directory containing project data.
cd .../


## Make sure fastq files are in a directory called 'data'

# Make a directory 'qiime' for all output files
mkdir -p qiime


# Create the manifest file - the below code expects that all reads will be in a dir called 'data'
echo -e "sample-id\tabsolute-filepath" > ./qiime/pacbio-manifest.tsv
for x in  data/*.gz; do
    echo -e `basename $x .fastq.gz`"\t"`realpath $x` >> pacbio-manifest.tsv
done


# Activate the version of QIIME2 we want to use.
source activate .../conda-envs/qiime2-amplicon-2024.10



# Import

$COMMAND qiime tools import \
  --type SampleData[SequencesWithQuality] \
  --input-path ./qiime/pacbio-manifest.tsv \
  --input-format SingleEndFastqManifestPhred33V2 \
  --output-path ./qiime/demux-pacbio.qza

# Summary

$COMMAND qiime demux summarize \
    --i-data ./qiime/demux-pacbio.qza \
    --o-visualization ./qiime/pacbio_reads_summary.qzv

# Denoise/Trimming - adjust accordingly

$COMMAND qiime dada2 denoise-ccs \
   --i-demultiplexed-seqs ./qiime/demux-pacbio.qza \
   --p-front AGRGTTYGATYMTGGCTCAG \
   --p-adapter  RGYTACCTTGTTACGACTT \
   --p-n-threads $num_threads \
   --o-representative-sequences ./qiime/seqs_pacbio.qza \
   --o-denoising-stats ./qiime/denoising-stats_pacbio.qza \
   --o-table ./qiime/table_pacbio.qza \
   --verbose > ./qiime/denoising_error.txt 2>&1

# # OTU Clustering - Change 0.97/0.98/0.99 etc

# $COMMAND qiime vsearch cluster-features-de-novo \
#  --i-table ./table_pacbio.qza \
#  --i-sequences ./seqs_pacbio.qza \
#  --p-perc-identity 0.99 \
#  --o-clustered-table ./table-dn-99_pacbio.qza \
#  --o-clustered-sequences ./rep-seqs-dn-99_pacbio.qza

# export  table

$COMMAND qiime tools export --input-path ./qiime/table_pacbio.qza --output-path ./qiime/table_pacbio_export

# export denoising

$COMMAND qiime tools export --input-path ./qiime/denoising-stats_pacbio.qza --output-path ./qiime/denoising-stats_pacbio_export



# Classification

## ASV

$COMMAND qiime feature-classifier classify-sklearn  \
 --i-classifier .../16S_datasets/2024.09.backbone.full-length.nb.sklearn-1.4.2.qza  \
 --i-reads  ./qiime/seqs_pacbio.qza  \
 --p-n-jobs -1 \
 --o-classification ./qiime/seqs-pacbio-class_ASV_GG_sklearn_NB.qza  \
 --verbose
  ## metadata
$COMMMAND qiime metadata tabulate \
  --m-input-file ./qiime/seqs-pacbio-class_ASV_GG_sklearn_NB.qza \
  --o-visualization ./qiime/seqs-pacbio-class_ASV_GG_sklearn_NB.qzv
### visualisation
$COMMAND qiime taxa barplot \
 --i-table ./qiime/table_pacbio.qza \
 --i-taxonomy ./qiime/seqs-pacbio-class_ASV_GG_sklearn_NB.qza \
 --o-visualization ./qiime/seqs-pacbio-class_ASV_GG_sklearn_NB_bar # This is the main output file.
####

# Export representative sequences

$COMMAND qiime tools export \
  --input-path ./qiime/seqs_pacbio.qza \
    --output-path ./qiime/seqs_pacbio_export

#$COMMAND qiime tools export \
#  --input-path ./qiime/seqs-dn-99_pacbio.qza \
#    --output-path ./qiime/seqs-dn-99_pacbio_export


echo "done"


# #OTU
# $COMMAND qiime feature-classifier classify-consensus-blast \
#  --i-query ./qiime/rep-seqs-dn-99_pacbio.qza \
#  --i-reference-reads .../16S_datasets/silva-138-99-seqs.qza \
#  --i-reference-taxonomy .../16S_datasets/silva-138-99-tax.qza \
#  --o-classification ./rep-seqs-pacbio-class_OTU_99_SILVA_BLAST.qza \
#  --verbose
# ## metadata
# $COMMMAND qiime metadata tabulate \
#   --m-input-file ./rep-seqs-pacbio-class_OTU_99_SILVA_BLAST.qza \
#   --o-visualization ./rep-seqs-pacbio-class_OTU_99_SILVA_BLAST.qzv
# ### visualisation
# $COMMAND qiime taxa barplot \
#  --i-table ./table_pacbio.qza \
#  --i-taxonomy ./rep-seqs-pacbio-class_OTU_99_SILVA_BLAST.qza \
#  --o-visualization ./rep-seqs-pacbio-class_OTU_99_SILVA_BLAST_bar
# ####

#  $COMMAND qiime feature-classifier classify-consensus-blast \
#  --i-query ./rep-seqs-dn-99_pacbio.qza \
#  --i-reference-reads .../16S_datasets/2024.09.backbone.full-length.nb.sklearn-1.4.2.qza \
#  --i-reference-taxonomy .../16S_datasets/2024.09.backbone.tax.qza \
#  --o-classification ./rep-seqs-pacbio-class_OTU_99_GG_BLAST.qza \
#  --verbose
# ## metadata
# $COMMMAND qiime metadata tabulate \
#   --m-input-file ./rep-seqs-pacbio-class_OTU_99_GG_BLAST.qza \
#   --o-visualization ./rep-seqs-pacbio-class_OTU_99_GG_BLAST.qzv
# ### visualisation
# $COMMAND qiime taxa barplot \
#  --i-table ./table_pacbio.qza \
#  --i-taxonomy ./rep-seqs-pacbio-class_OTU_99_GG_BLAST.qza \
#  --o-visualization ./rep-seqs-pacbio-class_OTU_99_GG_BLAST_bar
# ####


# $COMMAND qiime feature-classifier classify-sklearn  \
#  --i-classifier .../16S_datasets/silva-138-99-nb-classifier.qza  \
#  --i-reads  ./rep-seqs-dn-99_pacbio.qza  \
#  --p-n-jobs -1 \
#  --o-classification ./rep-seqs-pacbio-class_OTU_99_Silva_sklearn_NB.qza  \
#  --verbose
#  ## metadata
# $COMMMAND qiime metadata tabulate \
#   --m-input-file ./rep-seqs-pacbio-class_OTU_99_Silva_sklearn_NB.qza \
#   --o-visualization ./rep-seqs-pacbio-class_OTU_99_Silva_sklearn_NB.qzv
# ### visualisation
# $COMMAND qiime taxa barplot \
#  --i-table ./table_pacbio.qza \
#  --i-taxonomy ./rep-seqs-pacbio-class_OTU_99_Silva_sklearn_NB.qza \
#  --o-visualization ./rep-seqs-pacbio-class_OTU_99_Silva_sklearn_NB_bar
# ####

# $COMMAND qiime feature-classifier classify-sklearn  \
#  --i-classifier .../16S_datasets/2024.09.backbone.full-length.nb.sklearn-1.4.2.qza  \
#  --i-reads  ./rep-seqs-dn-99_pacbio.qza  \
#  --p-n-jobs -1 \
#  --o-classification ./rep-seqs-pacbio-class_OTU_99_GG_sklearn_NB.qza  \
#  --verbose
#  ## metadata
# $COMMMAND qiime metadata tabulate \
#   --m-input-file ./rep-seqs-pacbio-class_OTU_99_GG_sklearn_NB.qza \
#   --o-visualization ./rep-seqs-pacbio-class_OTU_99_GG_sklearn_NB.qzv
# ### visualisation
# $COMMAND qiime taxa barplot \
#  --i-table ./table_pacbio.qza \
#  --i-taxonomy ./rep-seqs-pacbio-class_OTU_99_GG_sklearn_NB.qza \
#  --o-visualization ./rep-seqs-pacbio-class_OTU_99_GG_sklearn_NB_bar
####

# #ASV

# $COMMAND qiime feature-classifier classify-consensus-blast \
#  --i-query ./rep-seqs_pacbio.qza \
#  --i-reference-reads .../16S_datasets/silva-138-99-seqs.qza \
#  --i-reference-taxonomy .../16S_datasets/silva-138-99-tax.qza \
#  --o-classification ./rep-seqs-pacbio-class_ASV_SILVA_BLAST.qza \
#  --verbose
#  ## metadata
# $COMMMAND qiime metadata tabulate \
#   --m-input-file ./rep-seqs-pacbio-class_ASV_SILVA_BLAST.qza \
#   --o-visualization ./rep-seqs-pacbio-class_ASV_SILVA_BLAST.qzv
# ### visualisation
# $COMMAND qiime taxa barplot \
#  --i-table ./table_pacbio.qza \
#  --i-taxonomy ./rep-seqs-pacbio-class_ASV_SILVA_BLAST.qza \
#  --o-visualization ./rep-seqs-pacbio-class_ASV_SILVA_BLAST_bar
# ####

#  $COMMAND qiime feature-classifier classify-consensus-blast \
#  --i-query ./rep-seqs_pacbio.qza \
#  --i-reference-reads .../16S_datasets/2024.09.backbone.full-length.nb.sklearn-1.4.2.qza \
#  --i-reference-taxonomy .../16S_datasets/2024.09.backbone.tax.qza \
#  --o-classification ./rep-seqs-pacbio-class_ASV_GG_BLAST.qza \
#  --verbose
#  ## metadata
# $COMMMAND qiime metadata tabulate \
#   --m-input-file ./rep-seqs-pacbio-class_ASV_GG_BLAST.qza \
#   --o-visualization ./rep-seqs-pacbio-class_ASV_GG_BLAST.qzv
# ### visualisation
# $COMMAND qiime taxa barplot \
#  --i-table ./table_pacbio.qza \
#  --i-taxonomy ./rep-seqs-pacbio-class_ASV_GG_BLAST.qza \
#  --o-visualization ./rep-seqs-pacbio-class_ASV_GG_BLAST_bar
# ####

# $COMMAND qiime feature-classifier classify-sklearn  \
#  --i-classifier .../16S_datasets/silva-138-99-nb-classifier.qza  \
#  --i-reads  ./rep-seqs_pacbio.qza  \
#  --p-n-jobs -1 \
#  --o-classification ./rep-seqs-pacbio-class_ASV_Silva_sklearn_NB.qza  \
#  --verbose
# ## metadata
# $COMMMAND qiime metadata tabulate \
#   --m-input-file ./rep-seqs-pacbio-class_ASV_Silva_sklearn_NB.qza \
#   --o-visualization ./rep-seqs-pacbio-class_ASV_Silva_sklearn_NB.qzv
# ### visualisation
# $COMMAND qiime taxa barplot \
#  --i-table ./table_pacbio.qza \
#  --i-taxonomy ./rep-seqs-pacbio-class_ASV_Silva_sklearn_NB.qza \
#  --o-visualization ./rep-seqs-pacbio-class_ASV_Silva_sklearn_NB_bar
####



