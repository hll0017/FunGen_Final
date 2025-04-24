# FunGen_Final
Final Project for Functional Genomics (BIOL 5850/6850) - RNAseq differential expression experiment in testes across dog breeds (large vs. small breeds)

## Data
Raw single-end RNAseq data from testicular tissue were retrieved from the NCBI Sequence Read Archive (SRA) under BioProject PRJNA69086. A total of 12 datasets from 12 different dog breeds were selected. These datasets were divided into two size groups based on breed size according to the [American Kennel Club](https://www.akc.org/dog-breeds/)'s definition of Large/XLarge and Small/XSmall dog breeds. Samples used are listed in the table below with associated SRR number, breed name and size group.

| SRR | Breed | Size Group |
| --- | --- | --- |
| SRR13389818 | Havanese | Small |
| SRR13389805 | Papillon | Small |
| SRR13389807 | Bichon Frise | Small |
| SRR13389809 | Maltese | Small |
| SRR13389841 | Border Terrier | Small |
| SRR13389834 | Yorkshire Terrier | Small |
| SRR13389820 | Golden Retriever | Large |
| SRR13389817 | Irish Wolfhound | Large |
| SRR13389815 | Labrador Retriever | Large |
| SRR13389849 | Rottweiler | Large |
| SRR13389840 | Weimaraner | Large |
| SRR13389824 | German Shepherd | Large |

## Download, Trimming, and Mapping Data
### Download
Reads were downloaded from NCBI Sequence Read Archive and raw read quality was assessed using the RNAseq_scripts/1_Download_QC.sh script on the Alabama Supercomputer (ASC). See script for recommended job parameters.

### Trimming
Raw read were trimmed and cleaned using the RNAseq_scripts/2_Trim_QC.sh script on the ASC. See script for recommended job parameters.

### Mapping
Trimmed reads were then mapped to the annotated dog reference genome Dog10K_Boxer_Tasha_1.0 GCF_000002285.5 and reads were counted using the RNAseq_scripts/3_Map_Count.sh on the ASC. See script for recommended job parameters. **This is an extremely memory intensive script, so parameters may need to be adjusted for script to run properly**

## Differential Gene Expression Analysis
Read counts (gene_count_matrix.csv) were imported into RStudio. The script DEG_Script.R was used to calculate differential expression and create a ranked gene list (DGErankName.rnk) which was used for subsequent GSEA and Cytoscape analyses. 

## Analysis of Insulin and Insulin-like Signaling (IIS) Pathway
### GSEA


### Cytoscape

## Analysis of Cancer-Associated Gene Sets
### GSEA

### Cytoscape

