####  DESeq2 Script for Differential Gene Expression Analysis 


#### Install the DESeq2 package 
install.packages("BiocManager")

#load required packages
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(ggplot2)
library(ggrepel)
library(DESeq2). ## Load the DESeq2 library 


## Confirm working directory - Source File Directory
getwd()

########## Input data   ##############

### Input the count data, the gene(/transcript) count matrix and labels
countdata <- as.matrix(read.csv("RNAseq_Counts/RNAseq/Counts_H_S_2025/gene_count_matrix.csv", row.names="gene_id"))
dim(countdata)
head(countdata)


### Input the meta data or phenotype data
##  Make sure the individual names match between the count data and the metadata
coldata <-(read.table("DEG_files/PHENODATA.txt", header=TRUE, row.names=1))
dim(coldata)
head(coldata)


#Check all sample IDs in colData are also in CountData and match their orders
all(rownames(coldata) %in% colnames(countdata))
countdata <- countdata[, rownames(coldata)]
all(rownames(coldata) == colnames(countdata))


## Create the DESEQ dataset and define the statistical model (page 6 of the manual)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData=coldata,  design = ~size)
#look at it
dds

#Check Library Size (Total Reads per Sample)
colSums(counts(dds))


#Visualize Library Size Distribution - optional
jpeg("Total_reads.jpg", width = 6, height = 6, units = "in", res = 300)
barplot(colSums(counts(dds)), las=2, main="Total Read Counts per Sample", col="steelblue")
dev.off()

#####   Prefiltering    
# Here we perform a minimal pre-filtering to remove rows that have less than 20 reads mapped.
dds <- dds[ rowSums(counts(dds)) > 20, ]
dds


## set factors for statistical analyses
###### The levels to our treatment names in the PHENO_DATA: small is the control, large is the comparison group
dds$condition <- factor(dds$size, levels=c("Small_breed","Large_breed"))



###### Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res


# We can order our results table by the smallest adjusted p value:
resOrdered <- res[order(res$padj),]
resOrdered

# We can summarize some basic tallies using the summary function 
summary(res)
#How many adjusted p-values were less than 0.1
sum(res$padj < 0.1, na.rm=TRUE)

#How many adjusted p-values were less than 0.05
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)



###    MA-plot
##plotMA shows the log2 fold changes attributable to a given variable over the meanof normalized counts. 

jpeg("MA_plot01.jpg", width = 6, height = 6, units = "in", res = 300)  ## using p < 0.01
plotMA(res, main="DESeq2", ylim=c(-8,8))
dev.off()

jpeg("MA_plot05.jpg", width = 6, height = 6, units = "in", res = 300)  ## using p < 0.05
plotMA(res05, main="DESeq2", ylim=c(-8,8))
dev.off()




#After calling plotMA, one can use the function identify to interactively detect the row number of individual genes by clicking on the plot. 
# One can then recover the gene identifiers by saving the resulting indices:
plotMA(res, main="DESeq2", ylim=c(-8,8))

idx <- identify(res$baseMean, res$log2FoldChange)
# after selecting a gene. You need to press escape to move on
rownames(res)[idx]


## Plot counts - sanity check!

# You can select the gene to plot by rowname or by numeric index.
plotCounts(dds, gene="gene-PTEN|PTEN", intgroup="size")
# You can plot the gene with th lowest adjusted P-value
plotCounts(dds, gene=which.min(res$padj), intgroup="size")
dds

##  Write your results to a file 
write.csv(as.data.frame(resOrdered), file="DGESeq_results.csv")  

## Extracting transformed values
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
head(assay(rld), 3)

### Heatmap of the count matrix
#library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)

library("pheatmap")
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("size", "tissue")])
df <- as.data.frame(colData(dds)[,c("size", "tissue")])
pheatmap(mat, annotation_col = anno)

#To save heatmap
jpeg("heatmapCM_Main.jpg", width = 10, height = 10, units = "in", res = 300)
pheatmap(mat, annotation_col = anno)
dev.off()


#check if factors follow desired order
#Reorder the factor levels so Small_breed comes before Big_breed if needed
rld$size <- factor(rld$size, levels = c("Small_breed", "Large_breed"))
print(rld$size)
str(rld$size)

# Check distant to distance bewteen samples
sampleDists <- dist(t(assay(rld)))

# Heatmap of the sample-to-sample distances
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$size)
colnames(sampleDistMatrix) <- paste(rld$size)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#To save heatmap
jpeg("heatmapDist.jpg", width = 10, height = 10, units = "in", res = 300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()




#  Principal component plot of the samples
plotPCA(rld, intgroup=c("size"))

#To save Principal component plot 
jpeg("PCAPlot.jpg", width = 10, height = 10, units = "in", res = 300)
plotPCA(rld, intgroup=c("size"))
dev.off()


# Create PCA plot object
p <- plotPCA(rld, intgroup = c("size"), returnData = TRUE)
percentVar <- round(100 * attr(p, "percentVar"))

# Create custom ggplot version of PCA
pca_plot <- ggplot(p, aes(PC1, PC2, color = size)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal(base_size = 16) +  # Adjust font size here
  theme(
    axis.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  ggtitle("Principal Component Analysis")

# Save to file
jpeg("PCAPlot.jpg", width = 10, height = 10, units = "in", res = 300)
print(pca_plot)
dev.off()


############ Preparing Data for GSEA and Cytoscape.  #############

### Merge 'gene names' with DGE results by Gene Model

## Import the DGE results file 
DGEresults <- read.csv("DEG_files/DGESeq_results.csv", stringsAsFactors = FALSE)
summary(DGEresults)
dim(DGEresults)

## Rename first column so it matches "gene_id" in annotation file
names(DGEresults)[1]<- "gene_id" 



############################# Make ranked list for GSEA ####################

DGE_Anno_Rank <-  within(DGEresults, rank <- sign(log2FoldChange) * -log10(pvalue))
DGE_Anno_Rank 

#subset the results so only Gene Name and rank
DGErank = subset(DGE_Anno_Rank, select = c(gene_id,rank) )
DGErank



### ...remove the "gene-" from row names
## https://stackoverflow.com/questions/39897155/replace-and-remove-part-of-string-in-rownames/39897315  "URS000075AF9C-snoRNA_GTATGTGTGGACAGCACTGAGACTGAGTCT"    to   "snoRNA"
## We can use gsub to match one of more characters that are not a - ([^-]+) from the start (^) of the string followed by 

# Extract gene names using gsub 
DGErank$gene_name <- gsub("gene-(.*)\\|.*", "\\1", DGErank$gene_id)
# Check the result
head(DGErank)


# Replace gene_name with NA where gene_id starts with "gene-LOC" or "rna-NC" or "rna-NR" to remove unnamed, non-coding or microRNAs
DGErank$gene_name[grepl("^gene-LOC|rna-NC|rna-NR", DGErank$gene_id)] <- NA
# View the result
head(DGErank)
sum(is.na(DGErank$gene_name))


#omit all NAs
DGErank_withName <- na.omit(DGErank)
DGErank_withName
dim(DGErank_withName)
head(DGErank_withName)


# Select only the gene_name and rank columns
DGErank_withName <- DGErank_withName[, c("gene_name", "rank")]

# View the result
head(DGErank_withName)
print(DGErank_withName$gene_name)


write.table(as.data.frame(DGErank_withName), file="DGErankName.rnk", quote=FALSE, row.names=FALSE, sep = "\t")  

##############  We also need the normalized expression DATA
## Obtain the transformed normalized count matrix
nt <- normTransform(dds) # defaults to log2(x+1)
head(assay(nt))
# compare to original count data
head(countdata)
# make the transformed normalized count matrix a new dataframe and make column 1 gene name
NormTransExp<-assay(nt)
summary(NormTransExp)
head(NormTransExp)
NormTransExp <- data.frame(gene_name = rownames(NormTransExp), NormTransExp, row.names = NULL)


# # Replace gene_name with NA where gene_id starts with "gene-LOC" or "rna-NC" or "rna-NR" to remove unnamed, non-coding or microRNAs
NormTransExp$gene_name[grepl("^gene-LOC|rna-NC|rna-NR", NormTransExp$gene_name)] <- NA
sum(is.na(NormTransExp$gene_name))

NormTransExp_Anno_withName <- na.omit(NormTransExp)
dim(NormTransExp_Anno_withName)

# Extract gene names using gsub 
NormTransExp_Anno_withName$gene_name <- gsub("gene-(.*)\\|.*", "\\1", NormTransExp_Anno_withName$gene_name)


## Write the transformed normalized count matrix with Gene Names to a tab delimited text file that can be imported into Cytoscape
write.table(as.data.frame(NormTransExp_Anno_withName), file="NormTransExp_Anno_Names.txt", quote=FALSE, row.names=FALSE, sep = "\t")  




# List of known top regulators IIS pathway related to growth
IIS_genes <- c("IGF1", "IGF2", "INS", "INSR", "IGF1R", "IGF2R",
               "IGFBP1", "IGFBP2", "IGFBP3", "IGFBP4", "IGFBP5", "IGFBP6",
               "IRS1", "IRS2", "IRS4")

# plotting list of top regulator of IIS pathway related to growth found in gene rank list
subset_IIS <- sorted_gene_data[sorted_gene_data$gene %in% IIS_genes, ]
ggplot(subset_IIS, aes(x = reorder(gene, -rank), y = rank)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Ranks of TOP Regulators of IIS Pathway", y = "Rank Score", x = "Gene")

## IIS genes enriched in GSEA
Insul_KEGG <- read.delim("GSEA_Results/my_GSEA_Insul_GeneEntry_analysis.GseaPreranked.1744336406675/from_text_entry_.tsv", stringsAsFactors = FALSE)
head(Insul_KEGG)
as.data.frame(Insul_KEGG)

#plot IIS genes enriched in GSEA
library(ggrepel)

ggplot(Insul_KEGG, aes(x = RANK.IN.GENE.LIST, y = RUNNING.ES)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(data = subset(Insul_KEGG, CORE.ENRICHMENT == "Yes"),
             aes(x = RANK.IN.GENE.LIST, y = RUNNING.ES),
             color = "red", size = 2) +
  geom_text_repel(data = subset(Insul_KEGG, CORE.ENRICHMENT == "Yes"),
                  aes(label = SYMBOL),
                  size = 3, max.overlaps = 20) +
  labs(title = "KEGG INSULIN SIGNALING PATHWAY",
       x = "Rank in Gene List",
       y = "Running ES") +
  theme_minimal()


####combination of insulin top regulators in ranked list and Insulin pathway################

# First plot (IIS Genes Rank)
plot1 <- ggplot(subset_IIS, aes(x = reorder(gene, -rank), y = rank)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = " ", y = "Rank Score", x = "Gene") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Second plot (KEGG Insulin Pathway)
plot2 <- ggplot(Insul_KEGG, aes(x = RANK.IN.GENE.LIST, y = RUNNING.ES)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(data = subset(Insul_KEGG, CORE.ENRICHMENT == "Yes"),
             aes(x = RANK.IN.GENE.LIST, y = RUNNING.ES),
             color = "red", size = 2) +
  geom_text_repel(data = subset(Insul_KEGG, CORE.ENRICHMENT == "Yes"),
                  aes(label = SYMBOL),
                  size = 3, max.overlaps = 20) +
  labs(title = " ",
       x = "Rank in Gene List",
       y = "Running ES") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Add A and B tags
plot1_labeled <- arrangeGrob(plot1, top = textGrob("A", x = unit(0, "npc"),
                                                   y = unit(1, "npc"), just = c("left", "top"),
                                                   gp = gpar(fontsize = 16, fontface = "bold")))
plot2_labeled <- arrangeGrob(plot2, top = textGrob("B", x = unit(0, "npc"),
                                                   y = unit(1, "npc"), just = c("left", "top"),
                                                   gp = gpar(fontsize = 16, fontface = "bold")))

# Combine plots
combined_IIS <- grid.arrange(plot1_labeled, plot2_labeled, ncol = 2)

# Create directory first
dir.create("GSEA_Results", showWarnings = FALSE)

# Combine plots
combined_IIS <- grid.arrange(plot1_labeled, plot2_labeled, ncol = 2)

# Save with ggsave instead
ggsave("GSEA_Results/enriched_IIS_genes.jpg", combined_IIS, width = 8, height = 6, dpi = 300)




########## gsea KEGG WHOLE- Result visualization  ########################

KEGG_Whole_up <- read.delim("GSEA_Results/my_GSEA_KEGG_whole_analysis.GseaPreranked.1744228114343/gsea_report_for_na_pos_1744228114343.tsv", stringsAsFactors = FALSE)
KEGG_Whole_dn <- read.delim("GSEA_Results/my_GSEA_KEGG_whole_analysis.GseaPreranked.1744228114343/gsea_report_for_na_neg_1744228114343.tsv", stringsAsFactors = FALSE)

# Select top 15 upregulated and top 30 downregulated with FDR.q.val < 0.05
top_up <- KEGG_Whole_up %>%
  filter(FDR.q.val < 0.25) %>%
  slice_max(NES, n = 20)

top_dn <- KEGG_Whole_dn %>%
  filter(FDR.q.val < 0.25) %>%
  slice_min(NES, n = 20)

# Add a direction column
top_up$Direction <- "Upregulated"
top_dn$Direction <- "Downregulated"

# Combine datasets
top_combined <- bind_rows(top_up, top_dn)

# Sort by NES within direction and update factor levels
top_combined <- top_combined %>%
  arrange(Direction, NES) %>%
  mutate(NAME = factor(NAME, levels = NAME))


# Plot
KEGG_up_dn <- ggplot(top_combined, aes(x = NAME, y = NES, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  #geom_text(aes(label = label),
  #position = position_dodge(width = 0.8),
  #hjust = ifelse(top_combined$Direction == "Up", -0.2, 1.2),
  #size = 4) +
  coord_flip() +
  labs(title = " ",
       x = NULL,
       y = "Normalized Enrichment Score (NES)") +
  scale_fill_manual(values = c("Upregulated" = "steelblue", "Downregulated" = "red")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5))


jpeg("GSEA_Results/KEGG_up_dn.jpg", width = 8, height = 6, units = "in", res = 300)
print(KEGG_up_dn)
dev.off()









######Part 2: Individual Component #################
######################## Gene Ontology Analysis 1 #######################################################

####################################################################################################



###################Gene Ontology using gseGO ###############
# Install required packages (only run once)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Cf.eg.db")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Cf.eg.db", "DOSE"))

# Load libraries
library(AnnotationDbi)
library(GO.db)
library(clusterProfiler)
library(org.Cf.eg.db)  # Dog gene annotation
library(DOSE)
library(dplyr)
library(tibble)

# Convert the gene list into a named vector
gene_list <- sorted_gene_data$rank
names(gene_list) <- sorted_gene_data$gene  # Assign gene names as the vector names

# Ensure the gene list is sorted in decreasing order
gene_list_sorted <- sort(gene_list, decreasing = TRUE)

# Run GSEA -BP
gsea_BP <- gseGO(geneList = gene_list_sorted, 
                      OrgDb = org.Cf.eg.db, 
                      keyType = "SYMBOL", 
                      ont = "BP", 
                      pvalueCutoff = 0.25)

# View results
head(gsea_BP)
as.data.frame(gsea_BP)
str(gsea_BP)

# to know number of positively and negatively enriched molecular functions

str(gsea_BP@result[gsea_BP@result$NES < 0, ])
sum(gsea_BP@result$NES < 0)


str(gsea_BP@result[gsea_BP@result$NES > 0, ])
sum(gsea_BP@result$NES > 0)

# Now plot
dotplot(gsea_BP, showCategory = 15) + ggtitle("GO Enrichment - Biological Process")



#save dotplots to files
jpeg("Gene_Ontology/GO_BP_dotplot.jpg", width = 8, height = 6.5, units = "in", res = 300)
dotplot(gsea_BP, showCategory = 15, title = "GO Enrichment - Biological Process")
dev.off()


# Check the actual term names first
head(gsea_BP@result$Description, 10)  # Or View(gsea_BP@result$Description)
length(gsea_BP@result$Description)

# Extract all genes enriched in GSEA insulin pathway analysis
insulin_genes <- Insul_KEGG$SYMBOL
core_insulin_genes <- Insul_KEGG[Insul_KEGG$CORE.ENRICHMENT == "Yes", "SYMBOL"]

# Extract core genes enriched in gseGO from gsea_BP
all_core_genes <- unique(unlist(strsplit(gsea_BP@result$core_enrichment, "/")))
print(all_core_genes)
length(all_core_genes)

IIS_gsea_BP <- gsea_BP@result[grep("insulin", gsea_BP@result$Description, ignore.case = TRUE), ]
IIS_gsea_BP
# You can add this information in a supplementary table

# Known IIS top regulators
IIS_genes <- c(
  "IGF1", "IGF2", "INS", "INSR", "IGF1R", "IGF2R",
  "IGFBP1", "IGFBP2", "IGFBP3", "IGFBP4", "IGFBP5", "IGFBP6",
  "IRS1", "IRS2", "IRS4"
)

# 1. Intersect BP core genes with GSEA-enriched IIS genes
intersect_gsea_iis <- intersect(all_core_genes, core_insulin_genes)

# 2. Intersect BP core genes with known IIS top regulators
intersect_known_iis <- intersect(all_core_genes, IIS_genes)

# Then for both intersections, extract Biological Process descriptions
map_genes_to_BP <- function(genes, gsea_result) {
  gene_to_bp <- list()
  
  for (i in seq_len(nrow(gsea_result))) {
    core_genes <- unlist(strsplit(gsea_result$core_enrichment[i], "/"))
    matched_genes <- intersect(genes, core_genes)
    
    if (length(matched_genes) > 0) {
      for (gene in matched_genes) {
        gene_to_bp[[gene]] <- c(gene_to_bp[[gene]], gsea_result$Description[i])
      }
    }
  }
  
  return(gene_to_bp)
}

# Apply function
BP_for_gsea_iis <- map_genes_to_BP(intersect_gsea_iis, gsea_BP@result)
BP_for_known_iis <- map_genes_to_BP(intersect_known_iis, gsea_BP@result)


# Make the tidy BP table first
tidy_BP_for_gsea_iis <- tibble(
  Gene = names(BP_for_gsea_iis),
  Biological_Process = vapply(BP_for_gsea_iis, function(x) paste(unique(x), collapse = "; "), character(1))
)

tidy_BP_for_gsea_iis <- tidy_BP_for_gsea_iis %>%
  mutate(Group = case_when(
    Gene %in% all_core_genes & Gene %in% core_insulin_genes & Gene %in% IIS_genes ~ "Top IIS regulators, GSEA-enriched biological process and IIS pathway",
    Gene %in% all_core_genes & Gene %in% core_insulin_genes ~ "GSEA-enriched biological process and IIS pathway",
    Gene %in% all_core_genes & Gene %in% IIS_genes ~ "GSEA-enriched biological process and top IIS regulators",
    Gene %in% core_insulin_genes & Gene %in% IIS_genes ~ "GSEA-enriched insulin pathway and top IIS regulators",
    Gene %in% all_core_genes ~ "GSEA-enriched biological process only",
    Gene %in% core_insulin_genes ~ "GSEA-enriched IIS pathway only",
    Gene %in% IIS_genes ~ "Top IIS regulators only",
    TRUE ~ "Other"
  ))

write.csv(tidy_BP_for_gsea_iis, "Gene_Ontology/IIS_comparison_description.csv", row.names = FALSE)



# Updated group list of Enriched BP genes, Enriched IIS pathway genes and Top IIS regulators related to growth
venn_list <- list(
  "Enriched BP genes" = all_core_genes,
  "Enriched IIS pathway genes" = core_insulin_genes,
  "Top IIS regulators" = IIS_genes
)


# Plot intersections to see overlap

library(ggvenn)
IIS_comparison <- ggvenn(
  venn_list, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
IIS_comparison

#save 
jpeg("Gene_Ontology/IIS_comparison.jpg", width = 8, height = 6.5, units = "in", res = 300)
print(IIS_comparison)
dev.off()



# Run GSEA -MF
gsea_MF <- gseGO(geneList = gene_list_sorted, 
                 OrgDb = org.Cf.eg.db, 
                 keyType = "SYMBOL", 
                 ont = "MF", 
                 pvalueCutoff = 0.25)

# View results
head(gsea_MF)
as.data.frame(gsea_MF)
print(gsea_MF)
# Now plot
dotplot(gsea_MF, showCategory = 15) + ggtitle("GO Enrichment - Molecular Function")
print(gsea_MF)

# to know number of positively and negatively enriched molecular functions
str(gsea_MF@result[gsea_MF@result$NES < 0, ])
sum(gsea_MF@result$NES < 0)


str(gsea_MF@result[gsea_MF@result$NES > 0, ])
sum(gsea_MF@result$NES > 0)

#save dotplots to files
jpeg("Gene_Ontology/GO_MF_dotplot.jpg", width = 8, height = 6.5, units = "in", res = 300)
dotplot(gsea_MF, showCategory = 15, title = "GO Enrichment - Molecular Function")
dev.off()



# Run GSEA -CC
gsea_CC <- gseGO(geneList = gene_list_sorted, 
                 OrgDb = org.Cf.eg.db, 
                 keyType = "SYMBOL", 
                 ont = "CC", 
                 pvalueCutoff = 0.25)

# View results
head(gsea_CC)
as.data.frame(gsea_CC)
print(gsea_CC)

# to know number of positively and negatively enriched molecular functions
str(gsea_CC@result[gsea_CC@result$NES < 0, ])
sum(gsea_CC@result$NES < 0)


str(gsea_CC@result[gsea_CC@result$NES > 0, ])
sum(gsea_CC@result$NES > 0)


# Now plot
dotplot(gsea_CC, showCategory = 15) + ggtitle("GO Enrichment - Cellular Component")
print(gsea_CC)

#save dotplots to files
jpeg("Gene_Ontology/GO_CC_dotplot.jpg", width = 8, height = 6.5, units = "in", res = 300)
dotplot(gsea_CC, showCategory = 15, title = "GO Enrichment - Cellular Component")
dev.off()



# Load required packages
library(gridExtra)

# Create the plots of combined  results
plot_BP <- dotplot(gsea_BP, showCategory = 15) + ggtitle("Biological Process")
plot_MF <- dotplot(gsea_MF, showCategory = 15) + ggtitle("Molecular Function")
plot_CC <- dotplot(gsea_CC, showCategory = 15) + ggtitle("Cellular Component")


# Add labels to the plots
plot_BP <- plot_BP + annotation_custom(grid::textGrob("A", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold")))
plot_MF <- plot_MF + annotation_custom(grid::textGrob("B", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold")))
plot_CC <- plot_CC + annotation_custom(grid::textGrob("C", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold")))

# Combine the plots into one image (2x2 grid, you can adjust as needed)
combined_plot <- grid.arrange(plot_BP, plot_MF, plot_CC, nrow = 1)  # Adjust nrow/ncol as needed

# Save the combined plot to a file
jpeg("Gene_Ontology/Combined_GO_Enrichment.jpg", width = 20, height = 10, units = "in", res = 300)
grid.arrange(plot_BP, plot_MF, plot_CC, nrow = 1)
dev.off()




# Get top 15 results as tables
top15_BP <- gsea_BP@result %>% 
  arrange(p.adjust) %>% 
  head(15)

top15_MF <- gsea_MF@result %>% 
  arrange(p.adjust) %>% 
  head(15)

top15_CC <- gsea_CC@result %>% 
  arrange(p.adjust) %>% 
  head(15)
write.csv(top15_BP, "Gene_Ontology/Top15_BP.csv", row.names = FALSE)
write.csv(top15_MF, "Gene_Ontology/Top15_MF.csv", row.names = FALSE)
write.csv(top15_CC, "Gene_Ontology/Top15_CC.csv", row.names = FALSE)









############# Optional- USing enrichGO for Overrepresented genes (to compare but not reported) ########

library(ggplot2)

# Make a copy of your results
volcano_data <- DGEresults

# Remove NA p-values
volcano_data <- volcano_data[!is.na(volcano_data$pvalue), ]

# Classify genes
volcano_data$Direction <- "Not Significant"
volcano_data$Direction[volcano_data$log2FoldChange > 1 & volcano_data$pvalue < 0.05] <- "Upregulated"
volcano_data$Direction[volcano_data$log2FoldChange < -1 & volcano_data$pvalue < 0.05] <- "Downregulated"

# Volcano plot
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue), color = Direction)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black")) +
  theme_minimal() +
  labs(title = " ",
       x = "Log2 Fold Change",
       y = "-log10(p-value)") +
  theme(plot.title = element_text(hjust = 0.5))



# 1. Install required packages (only run once)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Cf.eg.db", "DOSE"))

# 2. Load libraries
library(clusterProfiler)
library(org.Cf.eg.db)  # Dog gene annotation
library(DOSE)

#3. Load DESeq2 results
DEGRESULTS <- read.csv("DEG_files/DGESeq_results.csv", stringsAsFactors = FALSE)
head(DEGRESULTS)
DEGRESULTS$gene_name <- gsub(".*-(.*)\\|.*", "\\1", DEGRESULTS$X)
DEGRESULTS$gene_name[grepl("^LOC|rna-NC|rna-NR", DEGRESULTS$gene_name)] <- NA
sum(is.na(DEGRESULTS$gene_name))


# Remove rows with NA gene names
DEGRESULTS_clean <- na.omit(DEGRESULTS)


sig_genes <- DEGRESULTS_clean$gene_name[
  DEGRESULTS_clean$pvalue < 0.05 & abs(DEGRESULTS_clean$log2FoldChange) > 1
]
length(unique(sig_genes))
print(sig_genes)

up_genes <- DEGRESULTS_clean$gene_name[
  DEGRESULTS_clean$pvalue < 0.05 & DEGRESULTS_clean$log2FoldChange > 1
]
length(unique(up_genes))

down_genes <- DEGRESULTS_clean$gene_name[
  DEGRESULTS_clean$pvalue < 0.05 & DEGRESULTS_clean$log2FoldChange < -1
]
length(unique(down_genes))


# 6. Map SYMBOL to ENTREZ ID
up_gene_entrez <- bitr(up_genes, fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Cf.eg.db)

# 7. Use unique Entrez IDs
up_gene_ids <- unique(up_gene_entrez$ENTREZID)



# 8. Run GO enrichment for BP, MF, CC
up_go_bp <- enrichGO(gene         = up_gene_ids,
                     OrgDb        = org.Cf.eg.db,
                     keyType      = "ENTREZID",
                     ont          = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.25,
                     readable      = TRUE)
as.data.frame(up_go_bp)

# 9. Map SYMBOL to ENTREZ ID
down_gene_entrez <- bitr(down_genes, fromType = "SYMBOL",
                         toType = "ENTREZID",
                         OrgDb = org.Cf.eg.db)

# 10. Use unique Entrez IDs
down_gene_ids <- unique(down_gene_entrez$ENTREZID)



# 11. Run GO enrichment for BP, MF, CC
down_go_bp <- enrichGO(gene         = down_gene_ids,
                       OrgDb        = org.Cf.eg.db,
                       keyType      = "ENTREZID",
                       ont          = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.25,
                       readable      = TRUE)
as.data.frame(down_go_bp)

# Filter for IIS-related GO terms
IIS_up_go <- up_go_bp@result[grep("insulin", up_go_bp@result$Description, ignore.case = TRUE), ]
# View IIS-related terms
print(IIS_up_go)
as.data.frame(IIS_up_go)


up_go_mf <- enrichGO(gene         = up_gene_ids,
                     OrgDb        = org.Cf.eg.db,
                     keyType      = "ENTREZID",
                     ont          = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

down_go_mf <- enrichGO(gene         = down_gene_ids,
                       OrgDb        = org.Cf.eg.db,
                       keyType      = "ENTREZID",
                       ont          = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)

up_go_cc <- enrichGO(gene         = up_gene_ids,
                     OrgDb        = org.Cf.eg.db,
                     keyType      = "ENTREZID",
                     ont          = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

down_go_cc <- enrichGO(gene         = down_gene_ids,
                       OrgDb        = org.Cf.eg.db,
                       keyType      = "ENTREZID",
                       ont          = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)

# 9. Visualize top 15 GO terms with dotplots
barplot(up_go_bp, showCategory = 15, title = "GO Enrichment - Biological Process")
barplot(up_go_mf, showCategory = 15, title = "GO Enrichment - Molecular Function")
barplot(up_go_cc, showCategory = 15, title = "GO Enrichment - Cellular Component")

# 10. Save enrichment results as CSV
write.csv(as.data.frame(go_bp), "GO_BP_results.csv", row.names = FALSE)
write.csv(as.data.frame(go_mf), "GO_MF_results.csv", row.names = FALSE)
write.csv(as.data.frame(go_cc), "GO_CC_results.csv", row.names = FALSE)



# 11. Optionally save dotplots to files
jpeg("GO_BP_dotplot.jpg", width = 8, height = 6, units = "in", res = 300)
dotplot(go_bp, showCategory = 15, title = "GO Enrichment - Biological Process")
dev.off()

jpeg("GO_MF_dotplot.jpg", width = 8, height = 6, units = "in", res = 300)
dotplot(go_mf, showCategory = 15, title = "GO Enrichment - Molecular Function")
dev.off()

jpeg("GO_CC_dotplot.jpg", width = 8, height = 6, units = "in", res = 300)
dotplot(go_cc, showCategory = 15, title = "GO Enrichment - Cellular Component")
dev.off()













