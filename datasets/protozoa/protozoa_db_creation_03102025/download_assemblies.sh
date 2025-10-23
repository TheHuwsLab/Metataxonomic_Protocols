#Download ciliate genome assemblies

while read -r accession; do
    echo "Downloading $accession..."
    
    # Extract the numeric part: GCA_023806825.1 -> 023806825
    num=$(echo $accession | sed 's/GCA_//;s/\..*//')
    
    # Split into directory structure: 023/806/825
    dir1=$(echo $num | cut -c1-3)
    dir2=$(echo $num | cut -c4-6)
    dir3=$(echo $num | cut -c7-9)
    
    # Build the correct URL
    base_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/${dir1}/${dir2}/${dir3}"
    
    # Find the exact directory name (includes assembly name after accession)
    full_path=$(curl -s "${base_url}/" | grep -o "${accession}_[^/\"]*" | head -1)
    
    if [ -n "$full_path" ]; then
        wget "${base_url}/${full_path}/${full_path}_genomic.fna.gz"
    else
        echo "Could not find assembly for $accession"
    fi
    
done < assemblies.txt

#extract 
gunzip *.fna.gz