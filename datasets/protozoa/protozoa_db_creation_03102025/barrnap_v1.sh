#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --job-name=barrnap_ciliates_v1
#SBATCH --error=/users/3057233/ciliates/barrnap_ciliates_v1.error
#SBATCH --output=/users/3057233/ciliates/barrnap_ciliates_v1.txt
#SBATCH --partition=k2-bioinf,k2-hipri
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mail-user=k.lawther@qub.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL

# Activate the conda environment with barrnap installed
module load apps/anaconda3/2024.06/bin
source activate /mnt/scratch2/igfs-anaconda/conda-envs/barrnap

# Define directory containing genome files
INPUT_DIR="/users/3057233/sharedscratch/ciliates/assemblies"

# Loop through all genome files
for input in "$INPUT_DIR"/*_genomic.fna; do
    # Skip if no files are found
    [[ -e "$input" ]] || { echo "No assembly files found in $INPUT_DIR"; break; }

    # Define output filenames in same directory
    base=$(basename "$input" _genomic.fna)
    output_fasta="${INPUT_DIR}/${base}_rrna.fasta"
    output_gff="${INPUT_DIR}/${base}_rrna.gff"

    echo "Processing: $(basename "$input")"
    
    barrnap --kingdom euk --threads 20 --outseq "$output_fasta" "$input" > "$output_gff"

    # Extract just 18S sequences
    grep -A 1 "18S_rRNA" "$output_fasta" > "${INPUT_DIR}/${base}_18S.fasta"

    echo "Finished: $(basename "$input")"

done

echo "✅ All genome files in $INPUT_DIR processed!"