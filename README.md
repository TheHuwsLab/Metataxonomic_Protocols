# This repo contains the basic amplicon-based metataxonomic protocols used by the Huws Lab at Queen's University Belfast.

## ERRORS ARE EXPECTED!  - This is a work in progress and is not guaranteed to work on all systems or with all data types. Please report any issues you encounter - https://github.com/TheHuwsLab/Metataxonomic_Protocols/issues.

### Where possible, we use the QIIME2 framework and the provided sklearn models for GreenGenes2 for our analyses.
* Sklearn prebuilt models are downloaded from : https://library.qiime2.org/data-resources#naive-bayes-classifiers

### The protocols are designed to be run on a SLURM based HPC Linux system, but the commands are provided in a format suitable for a standard bash shell.

### The current version of the Illumina (short paired-end) protocol and PacBio (long single-end) protocol are setup to use version `qiime2-amplicon-2024.10'` of the QIIME2 pipeline (https://docs.qiime2.org/2024.10/install/index.html#).

# QIIME2 Illumina/PacBio Metataxonomic Protocol

These protocols automates the processing of Illumina paired-end and PacBio long-read sequencing data for metataxonomic analysis using QIIME2. It is designed to run on an HPC cluster with SLURM job scheduling.

## Features

- SLURM resource allocation 
- Automated manifest file creation for paired-end FASTQ files
- QIIME2 environment activation and data import
- Quality control, denoising, and feature table generation
- Taxonomic classification using a pre-trained sklearn model
- Export of key results and visualisations

## Usage

1. **Prepare Data**  
   Place your paired-end FASTQ files in a directory named `data/`. Files should be named as `sampleID_R1.fastq.gz` and `sampleID_R2.fastq.gz` for Illumina and `sampleID.fastq.gz` for PacBio. 

   
2. **Edit SLURM Directives**  
   Update SLURM parameters (partition, paths, email) to match your system.
   - Can also be run in a standard bash shell by removing the SLURM directives at the top of the script (and $COMMAND) .

3. **Activate QIIME2 Environment**  
   Ensure the specified QIIME2 conda environment exists.

4. **Run the Script**  
   Submit the script to your SLURM scheduler:

## Output

- `qiime/illumina-manifest.tsv`: Manifest file for QIIME2 import.
- `qiime/demux.qza`, `qiime/demux.qzv`: Imported and summarised data.
- `qiime/rep-seqs.qza`, `qiime/table.qza`: Denoised sequences and feature table.
- `qiime/sklearn-taxonomy_GG2_V4_ASV.qza`, `.qzv`: Taxonomic assignments and visualisation.
- Output includes feature tables, representative sequences, taxonomic assignments, and barplots.

## Notes

- Adjust memory, CPU, and time parameters for large datasets.
- Download the appropriate sklearn classifier from [QIIME2 Library](https://library.qiime2.org/data-resources#naive-bayes-classifiers).
- Review and modify trimming/truncation parameters based on your data quality.

## References

- [QIIME2 Documentation](https://docs.qiime2.org/)
- [QIIME2 Library](https://library.qiime2.org/)

---
This protocol is maintained by The Huws Lab.