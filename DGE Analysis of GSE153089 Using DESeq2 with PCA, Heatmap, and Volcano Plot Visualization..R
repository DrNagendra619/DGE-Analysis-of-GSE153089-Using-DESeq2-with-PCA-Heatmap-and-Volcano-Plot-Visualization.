# -------------------------------------------------------------------
# Project: Differential Gene Expression Analysis
# Dataset: GSE153089
# Title:   A microRNA cluster in the DLK1-DIO3 imprinted region on 
#          chromosome 14q32.2 in hepatoblastoma.
# -------------------------------------------------------------------

# Define the GEO Accession ID
geo_id <- "GSE153089"


#  Step 1: Package Installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "DESeq2", "pheatmap", "EnhancedVolcano", "ggplot2"), update = FALSE)


#  Step 2: Load Libraries
library(GEOquery)
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(ggplot2)


# Step 3: Fetch GEO Dataset
gse <- getGEO(geo_id, GSEMatrix = TRUE)
gse <- gse[[1]]

# Extract expression matrix & metadata
exprSet <- exprs(gse)
pdata <- pData(gse)

cat("Expression matrix dimensions:\n")
print(dim(exprSet))


# Step 4: Define Sample Condition (Corrected Logic)
#
# *** THIS IS THE FINAL FIX ***
# We use the 'characteristics_ch1' column and the exact text 
# "tissue: normal liver" to identify the normal samples.
# All other types ("fetal subtype", etc.) will be grouped as "Cancer".
#
pdata$condition <- ifelse(pdata$characteristics_ch1 == "tissue: normal liver", "Normal", "Cancer")

# Convert the condition to a factor
pdata$condition <- as.factor(pdata$condition)

# CRITICAL CHECK: This table MUST show two groups.
cat("\nSample Group Distribution (Must show 2 groups):\n")
print(table(pdata$condition))


# Step 5: Prepare DESeq2 Dataset
# (This step will now work because the table above shows 2 groups)

# Ensure matching order between expression and metadata
if (!all(colnames(exprSet) == rownames(pdata))) {
  stop("Sample names in exprSet and pdata do not match!")
}
exprSet <- round(exprSet) # Convert expression values to integer counts
# This line will now work
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = pdata,
                              design = ~ condition)
# Filter out low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]
cat("\nDESeqDataSet object created successfully.\n")


# Step 6: Run DESeq2 Pipeline
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj), ]
summary(res)


# Step 7: Variance Stabilizing Transformation (VST)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
cat("\nVST data dimensions:\n")
print(dim(assay(vsd)))


# Step 8: PCA Plot (Dimensional Reduction)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal(base_size = 14) +
  ggtitle(paste("PCA:", geo_id, "- Cancer vs Normal"))


# Step 9: Heatmap of Top 50 Variable Genes
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)

pheatmap(assay(vsd)[topVarGenes, ],
         cluster_rows = TRUE,
         show_rownames = FALSE,
         cluster_cols = TRUE,
         annotation_col = as.data.frame(colData(vsd)[, "condition", drop = FALSE]),
         main = paste(geo_id, ": Top 50 Variable Genes"))


# Step 10: Volcano Plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                title = paste(geo_id, ': Volcano Plot - Cancer vs Normal'),
                subtitle = 'Differential Gene Expression',
                legendLabels = c('NS', 'Log2FC', 'p-value', 'p-value & Log2FC'))
}