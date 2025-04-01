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


# Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$size)
colnames(sampleDistMatrix) <- NULL
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
