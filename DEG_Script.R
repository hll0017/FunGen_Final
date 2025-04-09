####  DESeq2 Script for Differential Gene Expression Analysis 


#### Install the DESeq2 package 
install.packages("BiocManager")

## Load the DESeq2 library 
library(DESeq2)


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
###### Note you need to change condition to treatment (to match our design above)
#  and levels to our treatment names in the PHENO_DATA: small is the null, big is the alternative group
dds$condition <- factor(dds$size, levels=c("Small_breed","Big_breed"))



###### Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res


# We can order our results table by the smallest adjusted p value:
resOrdered <- res[order(res$padj),]
resOrdered

# We can summarize some basic tallies using the summary function the default is p<0.1.
summary(res)
#How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)
#If the adjusted p value will be a value other than 0.1, alpha should be set to that value:
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)



###    MA-plot
##plotMA shows the log2 fold changes attributable to a given variable over the meanof normalized counts. 
## Points will be colored red if the adjusted p value is less than 0.1. 
## Points which fall out of the window are plotted as open triangles pointing either up or down

jpeg("MA_plot01.jpg", width = 6, height = 6, units = "in", res = 300)
plotMA(res, main="DESeq2", ylim=c(-8,8))
dev.off()

jpeg("MA_plot05.jpg", width = 6, height = 6, units = "in", res = 300)
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
jpeg("heatmapCM.jpg", width = 10, height = 10, units = "in", res = 300)
pheatmap(mat, annotation_col = anno)
dev.off()


#Reorder the factor levels so Small_breed comes before Big_breed
rld$size <- factor(rld$size, levels = c("Small_breed", "Big_breed"))
print(rld$size)
str(rld$size)
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


############## reorder for sample to sample distance Check the breed sizes and their structure
print(rld$size)
str(rld$size)

# Step 1: Reorder the factor levels so Small_breed comes before Big_breed
rld$size <- factor(rld$size, levels = c("Small_breed", "Big_breed"))

# Step 2: Get the correct order based on the newly ordered levels
breed_order <- order(rld$size)  # Now this should give the correct order: Small_breed first, then Big_breed
print(breed_order)

# Step 3: Reorder the sample distance matrix rows and columns according to breed_order
sampleDists <- dist(t(assay(rld)))  # Calculate the distances again if needed
sampleDistMatrix <- as.matrix(sampleDists)


# Reorder rows and columns based on breed_order
sampleDistMatrix <- sampleDistMatrix[breed_order, breed_order]

# Step 4: Set the row and column names to match the breed size order
rownames(sampleDistMatrix) <- rld$size[breed_order]
colnames(sampleDistMatrix) <- rld$size[breed_order]

# Step 5: Create the heatmap with the reordered distance matrix
library("RColorBrewer")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,  # Distance for rows
         clustering_distance_cols = sampleDists,  # Distance for columns
         col = colors)



# Reorder factor levels (if needed for other analyses)
rld$size <- factor(rld$size, levels = c("Small_breed", "Big_breed"))
str(rld$size)

# Calculate sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

# Create annotation data frame for the breeds
sample_names <- rownames(sampleDistMatrix)
sample_annotation <- data.frame(
  Breed = rld$size,
  row.names = sample_names
)

# Create color palette for breeds
ann_colors <- list(
  Breed = c(Small_breed = "#E69F00", Big_breed = "#56B4E9")
)

# Create the heatmap with clustering enabled and annotations
library("RColorBrewer")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,  # Enable clustering
         clustering_distance_cols = sampleDists,  # Enable clustering
         col = colors,
         annotation_row = sample_annotation,
         annotation_col = sample_annotation,
         annotation_colors = ann_colors)




# Create a data frame for column annotations
anno <- data.frame(Breed = rld$size)
rownames(anno) <- colnames(rld)

# Create the sample-to-sample distance matrix
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)

# Add proper names
rownames(sampleDistMatrix) <- colnames(rld)
colnames(sampleDistMatrix) <- colnames(rld)

# Plot with clustering but color bar to show breed
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         annotation_col = anno,
         annotation_row = anno,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))




########################### Reorder the annotation data based on breed size- CM
anno <- anno[order(anno$size), ]  # Reorder annotations so Big_breed is on one side and Small_breed on the other

# Reorder the matrix columns to match the annotation
mat <- mat[, rownames(anno)]  # Reorder mat columns to match the order of samples in anno

# Create the heatmap without specifying custom colors (it will auto-assign colors)
jpeg("heatmapCM_Modified.jpg", width = 10, height = 10, units = "in", res = 300)
pheatmap(mat, 
         annotation_col = anno,    # Add breed size annotation for columns
         cluster_rows = TRUE,      # Cluster rows based on gene expression
         cluster_cols = FALSE,     # Do not cluster columns (samples), since we already reordered them
         show_rownames = TRUE,     # Show row names (genes)
         show_colnames = TRUE)     # Show column names (samples)

dev.off()



#  Principal component plot of the samples
plotPCA(rld, intgroup=c("size"))

#To save Principal component plot 
jpeg("PCAPlot.jpg", width = 10, height = 10, units = "in", res = 300)
plotPCA(rld, intgroup=c("size"))
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
DGErank$gene_name[grepl("^gene-LOC", DGErank$gene_id)] <- NA
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
NormTransExp$gene_name[grepl("^gene-LOC", NormTransExp$gene_name)] <- NA
sum(is.na(NormTransExp$gene_name))

NormTransExp_Anno_withName <- na.omit(NormTransExp)
dim(NormTransExp_Anno_withName)

# Extract gene names using gsub 
NormTransExp_Anno_withName$gene_name <- gsub("gene-(.*)\\|.*", "\\1", NormTransExp_Anno_withName$gene_name)


## Write the transformed normalized count matrix with Gene Names to a tab delimited text file that can be imported into Cytoscape
write.table(as.data.frame(NormTransExp_Anno_withName), file="NormTransExp_Anno_Names.txt", quote=FALSE, row.names=FALSE, sep = "\t")  



######################## Gene Ontology Analysis #######################################################
library(AnnotationDbi)
library(GO.db)
library(clusterProfiler)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Cf.eg.db")

genes_to_test <- rownames(resOrdered[resOrdered$log2FoldChange > 0.5,])
head(genes_to_test)
# Clean the names to get only gene symbols
clean_gene_names <- gsub(".*-(.*)\\|.*", "\\1", genes_to_test)
head(clean_gene_names)

#1. Biological process
GO_BP_results <- enrichGO(gene = clean_gene_names, OrgDb = "org.Cf.eg.db", keyType = "SYMBOL",
                       ont = "BP",  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
as.data.frame(GO_BP_results)

library(enrichplot)
library(ggplot2)

dotplot(GO_BP_results , showCategory = 20) + ggtitle("Top 20 GO Biological Processes")

#2. molecular process
GO_MF_results <- enrichGO(gene = clean_gene_names, OrgDb = "org.Cf.eg.db", keyType = "SYMBOL",
                          ont = "MF",  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
as.data.frame(GO_MF_results)

dotplot(GO_MF_results , showCategory = 20) + ggtitle("Top 20 GO Molecular function")

#3.Cellular 
GO_CC_results <- enrichGO(gene = clean_gene_names, OrgDb = "org.Cf.eg.db", keyType = "SYMBOL",
                          ont = "CC",  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
as.data.frame(GO_CC_results)

dotplot(GO_CC_results , showCategory = 20) + ggtitle("Top 20 GO Cellular function")


####### Filter GO terms containing "insulin" in the description
insulin_related_go <- GO_results@result[grep("insulin", GO_results@result$Description, ignore.case = TRUE), ]

# View the filtered insulin-related GO terms
print(insulin_related_go)


# Extract insulin-related term IDs
insulin_ids <- grep("insulin", GO_results@result$Description, ignore.case = TRUE)
insulin_terms <- GO_results@result$ID[insulin_ids]

# Subset the enrichResult object while preserving its class
insulin_enrich <- GO_results[GO_results@result$ID %in% insulin_terms]

# Now plot
dotplot(insulin_enrich, showCategory = 10) + ggtitle("Top Insulin-Related GO Terms")



