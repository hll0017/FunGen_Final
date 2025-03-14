# !/bin/bash

# Goal: Make list of SRR IDs for the bioinformatic pipeline
# Goal: Make a table of SRR IDs, Breed Names, and Breed Size

# Extract Big Breeds' SRR IDs  and add to SRR_IDs.txt (exclude individuals marked with *)
awk -F'\t' ' NR > 1 && $0 !~ /\*/ { print $1 }' Big_Breeds_List.txt > SRR_IDs.txt

# Extract Small Breeds' IDs and append to SRR_IDs.txt (exclude individuals marked with *)
awk -F'\t' ' NR > 1 && $0 !~ /\*/ { print $1 }' Small_Breeds_List.txt >> SRR_IDs.txt

# Make table with SRR IDs, Breed Size, Breed Name, and Tissue Type
# Print header
printf "SRR\tName\tSize\tTissue\n" > Breed_Table.tsv

# Append Big Dog Breeds to Breed_Table.tsv
awk -F'\t' ' NR > 1 && $0 !~ /\*/ { print $1 FS $2 FS "Big" FS "Testes" }' Big_Breeds_List.txt >> Breed_Table.tsv

# Append Small Dog Breeds to Breed_Table.tsv
awk -F'\t' ' NR > 1 && $0 !~ /\*/ { print $1 FS $2 FS "Small" FS "Testes" }' Small_Breeds_List.txt >> Breed_Table.tsv
