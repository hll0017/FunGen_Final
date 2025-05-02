# FunGen_Final
Final Project for Functional Genomics (BIOL 5850/6850) - RNAseq differential expression and sequence variation in testes across dog breeds (large vs. small breeds)

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
Read counts (gene_count_matrix.csv) were imported into RStudio (v4.4.1). The script DEG_Script.R was used to calculate differential expression and create a ranked gene list (DGErankName.rnk) which was used for subsequent GSEA and Cytoscape analyses. 

## Analysis Insulin and Insulin-like Signaling (IIS)
### Gene Set Enrichment Analysis (GSEA)  
#### KEGG Pathways
GSEA was performed using ranked differentially expressed preranked gene list (DEGrankName.rnk) and the CP:KEGG_LEGACY Gene Set Collection (c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt) with the No_Collapse option selected and 1000 permutations. The rest of the default options were left unchanged. The resulting folder is GSEA_Resluts/my_GSEA_KEGG_whole_analysis.GseaPreranked.1744228114343.

### Insulin and Insulin-like Signaling (IIS) Pathway
GSEA was specifically performed using using ranked differentially expressed preranked gene list (DEGrankName.rnk) and a manually curated IIS pathway gene sets, with the No_Collapse option selected and 1000 permutations. The rest of the default options were left unchanged. The resulting folder is GSEA_Resluts/my_GSEA_Insul_GeneEntry_analysis.GseaPreranked.1744336406675

### Gene Ontology Enrichment 
GO enrichment analysis was performed using the ranked gene list (DEGrankName.rnk) with the gseaGO function from the clusterProfiler R package on R Studio (v4.4.1). The analysis was conducted with a p-value cutoff of < 0.25 and default settings, using the Canis familiaris annotation database (org.Cf.eg.db). Enrichment analysis was performed separately for the Biological Process (BP) and Molecular Function (MF) GO categories.

### Cytoscape
The results from GSEA of KEGG Pathways were imported into Cytoscape. Gene sets that met the threshold FDR < 0.1 were used to create a network representing gene set enrichment and similarity between gene sets.

### Sequence Variation 
Trimmed reads were mapped to a subset of IIS genes (IIS_CDS.fasta) using the script Map_VariantCall_DogIIS.sh on the ASC. See script for recommended job parameters

### Geneious
The resulting consensus sequences obtained from sequence variation were imported into the program Geneious. Consensus sequences were grouped based on gene and consensus sequence alignments were built for each gene. Each alignment was analyzed manually to identify shared sequence variations. For an SNP or amino acid change to be considered it must be shared by at least 3 samples and not contain ambiguous or missing nucleotides or amino acids.

## Analysis of Cancer-Associated Gene Sets
### GSEA
GSEA was performed using ranked differentially expressed gene list (DEGrankName.rnk), a treatment class file (treatment_final.cls), and the C2: Curated Gene Set Collection (c2.all.v2024.1.Hs.symbols.gmt) with the No_Collapse option selected and 1000 permutations. Gene sets smaller than 5 genes were excluded from GSEA. The resulting folder was C2.GseaPreranked.1744769312897.

### Cytoscape
The results from GSEA were imported into Cytoscape. A gene set network of C2 GSEA gene sets at an FDR < 0.1 was first created (C2.cys). From this file two additional networks were created. In the first, all nodes that had three or less neighbors within three connections were removed to identify connected gene sets (C2_Connections.cys, C2_UP.svg, and C2_DOWN.svg). In the second network, nodes were filtered out leaving only gene sets whose name contained the regular expression "CANCER" (C2_CANCER.cys and C2_CANCER.svg). For each network, all nodes were selected and all differentially expressed genes with a rank score > 12 were marked for further analysis.
