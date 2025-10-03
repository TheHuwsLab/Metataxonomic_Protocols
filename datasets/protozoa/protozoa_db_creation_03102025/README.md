# Removing redundancy in ciliate 18S dataset

This is the README for the directory that contains the code and intermediary files used to get v1.0.0 of the ciliate 18S dataset.

## Data acquisition

Data was acquired XX. 


## Step 1: Adding unique IDs

Unique IDs were added as lineage-based names were repeated in the dataset.
- Input: `ciliate_tax_ref_2025-04-01.fasta`

## Step 2: Clustering with CD-HIT

CD-HIT was used to cluster sequences at **100% identity and 100%
length** to identify duplicates.
- Output: `ciliate_tax_ref_ID100.clstr`

## Step 3: Manual inspection of clusters

The `.clstr` file was examined and for clusters containing multiple
sequences only the first sequence in the cluster was kept.
- The remaining identical sequences were flagged for removal (31 sequences in total).
- The sequences that were kept then underwent 1 of 3 options regarding their lineages

### Case 1: Where the full lineages were identical the name remains unchanged

    >Cluster 18
    Kept    0   1294nt, >090_k__Eukaryota;...;s__
    Removed 1   1294nt, >091_k__Eukaryota;...;s__

### Case 2: If only the SAG numbers differed (the rest of the lineage is identical) then the SAG number was removed

    >Cluster 2
    Renamed 0   1640nt, >186_k__Eukaryota;...Isotricha_intestinalis_SAG3 
    Removed 1   1640nt, >187_k__Eukaryota;...Isotricha_intestinalis_SAGT1
    Removed 2   1640nt, >188_k__Eukaryota;...Isotricha_intestinalis_SAGT2

Kept sequence renamed from:

    >186_k__...;s__Isotricha_intestinalis_SAG3

to:

    >186_k__...;s__Isotricha_intestinalis_SAG

### Case 3: Where the species name differed it was removed and generalised to `s__`

    >Cluster 9
    Renamed 0   1488nt, >123_k__Eukaryota;...Ostracodinium_gracile *
    Removed 1   1488nt, >128_k__Eukaryota;...Ostracodinium_trivesiculatum

Kept sequence renamed from:

    >123_k__...;s__Ostracodinium_gracile

to:

    >123_k__...;s__

## Step 4: Remove the redundant sequences

-   Mapping file: `remove.txt`
-   Output: `ciliate_tax_ref_ID_filtered.fasta`

## Step 5: Apply name changes

-   Mapping file: `rename.txt`
-   Output: `ciliate_tax_ref_ID_filtered_renamed.fasta`

## Step 6: Remove the unique IDs added in Step 1

-   Final output: `ciliate_tax_ref_ID_filtered_renamed_noID.fasta`

