# script to perform differential gene expression analysis using DESeq2 package
# setwd("~/Desktop/demo/DESeq2_tutorial/data")

# load libraries
library(DESeq2)
library(tidyverse)
library(airway)

# Step 1: preparing count data ----------------

# read in counts data
counts_data <- read.csv('counts_data.csv')
head(counts_data)


# read in sample info
colData <- read.csv('sample_info.csv')


# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# are they in the same order?
all(colnames(counts_data) == rownames(colData))


# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ dexamethasone)

dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds

# set the factor level
dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated")

# NOTE: collapse technical replicates

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)

res



# Explore Results ----------------

summary(res)

res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# contrasts
resultsNames(dds)


# Check available comparisons
resultsNames(dds)

# Get results for treated vs untreated (default)
res <- results(dds)

# OR explicitly specify the contrast
res <- results(dds, contrast = c("dexamethasone", "treated", "untreated"))



# MA plot
plotMA(res)


# Filter for significant genes (padj < 0.05 and |log2FC| > 1)
res_sig <- subset(res0.01, padj < 0.01 & abs(log2FoldChange) > 1)
res_sig

# Order by adjusted p-value
res_ordered <- res[order(res$padj),]
head(res_ordered, 20)

# Export results
write.csv(as.data.frame(res_ordered), file="DESeq2_results_all.csv")
write.csv(as.data.frame(res_sig), file="DESeq2_results_significant.csv")



library(EnhancedVolcano)

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Treated vs Untreated',
                pCutoff = 0.01,
                FCcutoff = 1)


library(pheatmap)

# Get normalized counts
vsd <- vst(dds, blind=FALSE)

# Select top 50 genes by adjusted p-value
top_genes <- head(order(res$padj), 50)
mat <- assay(vsd)[top_genes,]

# Create heatmap
pheatmap(mat, 
         cluster_rows=TRUE, 
         cluster_cols=TRUE,
         show_rownames=FALSE,
         scale="row")


plotPCA(vsd, intgroup="dexamethasone")




# remove this line (it will always error for this dataset)
# results(dds, contrast = c("dexamethasone", "treated_4hrs", "untreated"))

# correct explicit contrast (optional, same as your current res)
res <- results(dds, contrast = c("dexamethasone", "treated", "untreated"))

# filter and export
res_sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(res_sig), "DESeq2_significant_genes.csv")


# install and load org.Hs.eg.db in one go

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)


library(clusterProfiler)
library(org.Hs.eg.db)

# 1) Prepare a significant gene list (adjust cutoffs as you like)
res_sig <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)

# 2) Map Ensembl IDs to Entrez IDs and symbols
anno <- bitr(rownames(res_sig),
             fromType = "ENSEMBL",
             toType   = c("ENTREZID","SYMBOL"),
             OrgDb    = org.Hs.eg.db)

res_sig_anno <- merge(as.data.frame(res_sig),
                      anno,
                      by.x = "row.names",
                      by.y = "ENSEMBL")

# 3) GO enrichment (Biological Process)
ego <- enrichGO(gene          = res_sig_anno$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05)

# 4) View top enriched terms
head(as.data.frame(ego))


