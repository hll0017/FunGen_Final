# FunGen_Final
Final Project for Functional Genomics (BIOL 5850/6850) - RNAseq differential expression experiment in testes across dog breeds (large vs. small breeds)

## Data
Raw single-end RNAseq data from testicular tissue were retrieved from the NCBI Sequence Read Archive (SRA) under BioProject PRJNA69086. A total of 12 datasets from 12 different dog breeds were selected. These datasets were divided into two size groups based on breed size according to the [American Kennel Club](https://www.akc.org/dog-breeds/)'s definition of Large/XLarge and Small/XSmall dog breeds. Samples used are listed in the table below with associated SRR number, breed name and size group.

| SRR | Breed | Size Group |
| --- | --- | --- |


## Download, Trimming, and Mapping Data
### Download
Datasets were downloaded from NCBI Sequence Read Archive and raw read quality was assessed using the RNAseq_scripts/1_Download_QC.sh script on the Alabama Supercomputer (ASC). See script for recommended job parameters.

### Trimming
Raw datasets were trimmed and cleaned using the RNAseq_scripts/2_Trim_QC.sh script on the ASC. See script for recommended job parameters.

###

## Differential Gene Expression Analysis

## GSEA
### Insulin and Insulin-like Signaling (IIS)

### Cancer
