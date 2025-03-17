# scp -i ~/.ssh/March13.pem *R ubuntu@54.163.197.71:/home/ubuntu
# scp -i ~/.ssh/March13.pem ubuntu@54.163.197.71:/home/ubuntu/*png ./

library(Seurat)
library(dplyr)
library(ggplot2)
library(edgeR)
library(limma)
library(clusterProfiler)
library(ggrepel)
# library(AnnotationDbi)
# library(org.Hs.eg.db)
# library(parallel)
library(EnsDb.Hsapiens.v86)
library(foreach)
library(doParallel)

###################################################################################
###################################################################################

packageVersion("Seurat")

###################################################################################
###################################################################################

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

result <- mclapply(1:100, function(x) {
  Sys.sleep(1)  
  x^2
}, mc.cores = num_cores)

# print(result)
# stopCluster(cl)  # Stop cluster after completion
# print(result)

###################################################################################
###################################################################################

# Load the RDS file
test <- readRDS("dataset-test.RDS")

# Check metadata and object class
print(test)
class(test)

# Set the prefix for filenames
NAME <- "gs"

###################################################################################
###################################################################################

# str(test)
colnames(test@meta.data)

# Count unique sample IDs
unique(test@meta.data$Sample_ID)
# [1] Sample_1  Sample_15 Sample_2  Sample_3  Sample_4  Sample_5  Sample_6 
# [8] Sample_16 Sample_7  Sample_8  Sample_9  Sample_10 Sample_11 Sample_12
# [15] Sample_13 Sample_14
# 16 Levels: Sample_1 Sample_2 Sample_3 Sample_4 Sample_5 Sample_6 ... Sample_16

# Number of cells in each sample:
cell_counts <- as.data.frame(table(test@meta.data$Sample_ID))
colnames(cell_counts) <- c("Sample_ID", "Cell_Count")
print(cell_counts)

#   Sample_ID Cell_Count
# 1   Sample_1       3073
# 2   Sample_2       2995
# 3   Sample_3       2462
# 4   Sample_4       1453
# 5   Sample_5       2858
# 6   Sample_6       3878
# 7   Sample_7       2186
# 8   Sample_8        521
# 9   Sample_9        630
# 10 Sample_10       1496
# 11 Sample_11        945
# 12 Sample_12       2628
# 13 Sample_13       2514
# 14 Sample_14       3053
# 15 Sample_15       2835
# 16 Sample_16       1782

###################################################################################
###################################################################################

# GetAssayData(test, assay = "RNA", slot = "counts")  # Raw counts
# GetAssayData(test, assay = "RNA", slot = "data")    # Normalized data
# GetAssayData(test, assay = "RNA", slot = "scale.data")  # Scaled data 

# We detected only the raw counts, not the normalized or scaled data.

###################################################################################
###################################################################################

# There are several Ensembl ID Conversion Methods
# org.Hs.eg.db (Method 1)	
# AnnotationDbi (Method 2)
# biomaRt (Method 3)
# EnsDb.Hsapiens.v86 (Method 4)

###################################################################################
################################################################################### convert Ensembl ID into Gene Symbol

# Convert Ensembl IDs to Gene Symbols
gene_conversion <- ensembldb::select(EnsDb.Hsapiens.v86,
                                     keys = rownames(test),
                                     keytype = "GENEID",
                                     columns = c("SYMBOL"))

gene_map <- setNames(gene_conversion$SYMBOL, gene_conversion$GENEID)
updated_gene_names <- gene_map[rownames(test)]

# Keep original Ensembl ID if no mapping exists
rownames(test) <- ifelse(!is.na(updated_gene_names), updated_gene_names, rownames(test))

# Identify duplicate gene symbols
duplicate_genes <- duplicated(rownames(test))

# Print warning if duplicates exist
if (sum(duplicate_genes) > 0) {
  cat("Warning:", sum(duplicate_genes), "duplicated gene symbols found after conversion!\n")
  duplicate_rows <- duplicated(rownames(test)) | duplicated(rownames(test), fromLast = TRUE)
  test <- test[!duplicate_rows, ]   # Remove the duplicated rows
  cat("Duplicates removed. Total genes remaining:", nrow(test), "\n")
} else {
  cat("No duplicated gene symbols found after conversion.\n")
}

# Ensure remaining Ensembl IDs are converted to a a name : Unmapped___
# ensembl_ids_remaining <- grepl("^ENSG", rownames(test))
# rownames(test)[ensembl_ids_remaining] <- paste0("Unmapped_", rownames(test)[ensembl_ids_remaining])

# Final Check : Check for remaining duplicate gene names
final_duplicate_genes <- sum(duplicated(rownames(test)))
if (final_duplicate_genes > 0) {
  cat("Warning: There are still", final_duplicate_genes, "duplicated gene symbols remaining!\n")
} else {
  cat("Final check passed: No duplicated gene symbols remain.\n")
}

# Check if any genes remain in Ensembl format
num_ensembl_ids <- sum(grepl("^ENSG", rownames(test)))
total_genes <- nrow(test)
cat("Number of genes still using Ensembl IDs:", num_ensembl_ids, "\n")
cat("Total number of genes in dataset after cleanup:", total_genes, "\n")
head(rownames(test))

# Number of genes still using Ensembl IDs: 2598 
# Total number of genes in dataset after cleanup: 36465 

###################################################################################
################################################################################### check if normalized and scaled data available 

# Check if normalized data is available
if ("data" %in% names(test@assays$RNA)) {
  cat("Normalized data is available.\n")
  print(dim(GetAssayData(test, assay = "RNA", slot = "data")))
} else {
  cat("Normalized data is NOT available.\n")
}

# Check if scale.data is available
if ("scale.data" %in% names(test@assays$RNA)) {
  cat("Scaled data is available.\n")
  print(dim(GetAssayData(test, assay = "RNA", slot = "scale.data")))
} else {
  cat("Scaled data is NOT available.\n")
}

# Normalized data is NOT available.
# Scaled data is NOT available.

###################################################################################
###################################################################################

# Search for mitochondrial genes with different capitalization patterns
mt_genes <- grep("^MT-|^mt-|^Mt-", rownames(test), value = TRUE)

# Check if mitochondrial genes were found and print appropriate message
if (length(mt_genes) > 0) {
  cat("Yes, I have found these mitochondrial genes:\n")
  print(mt_genes)
} else {
  cat("I did not find any mitochondrial genes.\n")
}
###################################################################################
################################################################################### remove the genes with an expression that is 0

# Get the counts matrix from the RNA assay using the new `layer` argument
counts <- GetAssayData(test, assay = "RNA", layer = "counts")

# Identify genes that have non-zero expression in at least one cell
expressed_genes <- rownames(counts)[Matrix::rowSums(counts) > 0]

# Subset the Seurat object to keep only these expressed genes
test <- subset(test, features = expressed_genes)

cat("Number of genes remaining (with non-zero expression):", length(expressed_genes), "\n")
# Number of genes remaining (with non-zero expression): 32507

###################################################################################
###################################################################################

print("As the normalized data and the scaled data is not available, we perform the analysis from the raw counts data")

# Compute quantiles for nFeature_RNA
quantile_nFeature <- quantile(test@meta.data$nFeature_RNA, probs = c(0.25, 0.5, 0.75))
# Compute quantiles for nCount_RNA
quantile_nCount <- quantile(test@meta.data$nCount_RNA, probs = c(0.25, 0.5, 0.75))
print(paste("nFeature_RNA - 25th percentile:", quantile_nFeature[1], 
            "| Median (50th percentile):", quantile_nFeature[2], 
            "| 75th percentile:", quantile_nFeature[3]))
print(paste("nCount_RNA - 25th percentile:", quantile_nCount[1], 
            "| Median (50th percentile):", quantile_nCount[2], 
            "| 75th percentile:", quantile_nCount[3]))

# "nFeature_RNA - 25th percentile: 2259 | Median (50th percentile): 3414 | 75th percentile: 5701"
# "nCount_RNA - 25th percentile: 4680 | Median (50th percentile): 8639 | 75th percentile: 19212"

###################################################################################
################################################################################### percent MT

# MT genes 
test[["percent.mt"]] <- PercentageFeatureSet(test, pattern = "^MT-")

# QC metrics
png(paste0(NAME, "_VlnPlot_QC.png"), width = 1800, height = 600)
print(VlnPlot(test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
dev.off()

###################################################################################
################################################################################### QC

# Scatter plot for outliers (UMIs vs. Genes per Cell)
png(paste0(NAME, "_ScatterPlot_UMIs_vs_Genes.png"), width = 2000, height = 600)
plot1 <- FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

# Check the new dimensions before filtering
print(dim(test))

###################################################################################
################################################################################### filtering 

# Define filtering thresholds : I have selected the filtering thresholds based on 
# 25th percentile and 75th percentile.

test = subset(test, 
                        subset = nFeature_RNA > 1000 & nFeature_RNA < 9000 & 
                                 nCount_RNA > 1000 & nCount_RNA < 40000 )

# Check the new dimensions after filtering
print(dim(test))
# 32507 18745

###################################################################################
################################################################################### normalization

test <- NormalizeData(test, normalization.method = "LogNormalize", scale.factor = 10000)
test <- FindVariableFeatures(test, selection.method = "vst", nfeatures = 1000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(test), 10)
write.table(VariableFeatures(test), file = paste0(NAME, "_variable_features.txt"), 
            quote = FALSE, row.names = TRUE, col.names = TRUE)

plot1 <- VariableFeaturePlot(test)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
combined_plot <- plot1 + plot2

ggsave(filename = paste0(NAME, "_VariableFeatures_top10.png"), 
                         plot = combined_plot, 
                         width = 20, height = 12, dpi = 300)

###################################################################################
################################################################################### standardization

all.genes <- rownames(test)
test <- ScaleData(test, features = all.genes)
print(length(rownames(test)))

###################################################################################
################################################################################### PCA

test <- RunPCA(test, features = VariableFeatures(object = test))
print(test[["pca"]], dims = 1:5, nfeatures = 5)

# Save the Elbow Plot 
png(paste0(NAME, "_ElbowPlot.png"), width = 600, height = 600)
print(ElbowPlot(test, ndims = 50))
dev.off()

# Save PCA loadings
png(paste0(NAME, "_PCA_Loadings.png"), width = 600, height = 600)
print(VizDimLoadings(test, dims = 1:2, reduction = "pca"))
dev.off()

# Save PCA DimPlot
png(paste0(NAME, "_PCA_DimPlot.png"), width = 600, height = 600)
print(DimPlot(test, reduction = "pca") + NoLegend())
dev.off()

################################################################################### Neighbours
################################################################################### Clusters

# Find neighbors and clusters
test <- FindNeighbors(test, dims = 1:10)
test <- FindClusters(test, resolution = 0.5)
print(table(Idents(test)))

################################################################################### UMAP
################################################################################### TSNE

test <- RunUMAP(test, dims = 1:10)
test <- RunTSNE(test, dims = 1:10)

# Save UMAP clustering plot
png(paste0(NAME, "_UMAP_Clustering.png"), width = 800, height = 600)
print(DimPlot(test, reduction = "umap", label = TRUE) + ggtitle("UMAP Clustering"))
dev.off()

# Save t-SNE clustering plot
png(paste0(NAME, "_TSNE_Clustering.png"), width = 800, height = 600)
print(DimPlot(test, reduction = "tsne", label = TRUE) + ggtitle("t-SNE Clustering"))
dev.off()

################################################################################### 
################################################################################### FindAllMarkers
 
test_markers <- FindAllMarkers(test, only.pos = TRUE, 
                                     min.pct = 0.25, 
                                     logfc.threshold = 0.25, 
                                     # test.use = "DESeq2",
                                     test.use = "wilcox")

head(test_markers)
dim(test_markers)
test_markers$name <- rownames(test_markers)

# Top markers per cluster
top_10markers <- test_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
print(top_10markers)

write.table(top_10markers, 
            file = paste0(NAME, "_test_cluster_markers_top10.txt"), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

write.table(test_markers, 
            file = paste0(NAME, "_test_cluster_markers.txt"), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Summarize the number of significant genes per cluster (p_val_adj < 0.01)
sig_counts <- test_markers %>%
  group_by(cluster) %>%
  summarize(NumSigGenes = sum(p_val_adj < 0.01, na.rm = TRUE))

# Write the summary to a tab-delimited text file, prefixed by the NAME variable
write.table(sig_counts, file = paste0(NAME, "_significant_genes_per_cluster.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Summary of significant genes per cluster saved to file:", 
    paste0(NAME, "_significant_genes_per_cluster.txt"), "\n")

###################################################################################
################################################################################### Markers per cluster

# Saving the files containing the markers per cluster
uniquec <- unique(test_markers$cluster)

# Loop through each cluster and save markers separately
for (acluster in uniquec) {
  cluster_markers <- subset(test_markers, cluster == acluster)
  
  write.table(cluster_markers, 
              file = paste0(NAME, "_Cluster_", acluster, "_Markers.txt"), 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE)
}

################################################################################### Genes of interest
################################################################################### Violin, Feature Plots 

png(paste0(NAME, "_ViolinPlot_genes_interest.png"), width = 1200, height = 1200)
print(VlnPlot(test, features = c("GABRQ", "ADRA1A", "FEZF2", "VAT1L"), ncol = 2))
dev.off()

png(paste0(NAME, "_FeaturePlot_genes_interest.png"), width = 1200, height = 1200)
print(FeaturePlot(test, features = c("GABRQ", "ADRA1A", "FEZF2", "VAT1L")))
dev.off()

png(paste0(NAME, "_Heatmap_genes_interest.png"), width = 1200, height = 400)
print(DoHeatmap(test, features = c("GABRQ", "ADRA1A", "FEZF2", "VAT1L")) + NoLegend())
dev.off()

###################################################################################
###################################################################################

# Optionally, save your entire workspace to a file
# save.image(file = paste0(NAME, "_test_gdrive.RData"))

###################################################################################
###################################################################################

# Create a table of the number of cells per cluster
cluster_counts <- table(Idents(test))

write.table(cluster_counts, file = paste0(NAME, "_cluster_cell_counts.txt"), 
            quote = FALSE, col.names = NA)

#   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
# 7044 3700 3483 3100 2043 1815 1610 1506 1304 1184 1142  923  885  697  664  589 
#  16   17   18   19   20   21 
# 388  362  210  174  158   28 

################################################################################### to select the clusters for downstream analysis
################################################################################### function of "GABRQ", "ADRA1A", "FEZF2", "VAT1L"

# features : c("GABRQ", "ADRA1A", "FEZF2", "VAT1L")

# to count, for each cluster, the number of cells that express each of these genes 
# (an expression value greater than zero) ; print and save the results to a file.

features <- c("GABRQ", "ADRA1A", "FEZF2", "VAT1L")

# Fetch expression data for these genes and add the cluster identities
test_df <- FetchData(test, vars = features)
test_df$Cluster <- Idents(test)

# For each cluster, count the number of cells with expression > 0 for each gene
cluster_4genes_summary <- test_df %>%
  group_by(Cluster) %>%
  summarize(
    n_GABRQ = sum(GABRQ > 0),
    n_ADRA1A = sum(ADRA1A > 0),
    n_FEZF2 = sum(FEZF2 > 0),
    n_VAT1L = sum(VAT1L > 0),
    total_cells = sum(GABRQ > 0 | ADRA1A > 0 | FEZF2 > 0 | VAT1L > 0)
  ) %>%
  arrange(desc(total_cells))

# Print the result to the console
print(cluster_4genes_summary)

#    Cluster n_GABRQ n_ADRA1A n_FEZF2 n_VAT1L total_cells
#        1       1     2074     196     124        2205
#        3      26     1464       9    1080        1980
#        4       8      392     571     895        1316
#        8       1      816       2     356         958
#        6       3      610       4     304         757
#       11       0      273     387     522         739
#        2       5      197      14     547         730
#       14       0      509       1     393         589
#       12      11      140       4     383         460
#       7       7       31      39     391         439
#      15       4       25     215     322         426
#      10       3       17      42     353         380
#       9       2       72       2     307         361
#      13       0      342      19       5         355
#      16      31      123     204     163         308
#      18       5      152      94     165         173
#       0      12       39      13      95         158
#      17       0        6       3      48          56
#       5       0        7       6      14          27
#      19       0        1       0       1           2
#      20       0        0       1       0           1
#      21       0        1       0       0           1

# Save the results to a file prefixed with NAME (e.g., "gs_gene_cluster_counts.csv")
write.table(cluster_4genes_summary, 
            file = paste0(NAME, "_cluster_counts_for_4markers.txt"),
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

################################################################################### Clusters of interest 
################################################################################### 1, 3, 4

# There are few cells that express all 4 markers markers in the same cell.
# Because the clusters 1, 3 and 4 have more than 1000 cells, 
# we decide to progress the analysis with clusters 1, 3, 4. 

# Define the clusters of interest as numeric
selected_clusters <- c(1, 3, 4)

###################################################################################
###################################################################################
###################################################################################
################################################################################### aggregate expression

pseudo <- AggregateExpression(
  object = test,
  assays = "RNA",
  return.seurat = TRUE,
  group.by = "seurat_clusters"
)

pseudo_df <- AggregateExpression(
  object = test,
  assays = "RNA",
  return.seurat = FALSE,
  group.by = "seurat_clusters"
)

aggregated_matrix <- pseudo_df[["RNA"]]
head(aggregated_matrix)
dim(aggregated_matrix)

write.table(aggregated_matrix, file = paste0(NAME, "_test_aggregated_expression.txt", sep = "\t", row.names = TRUE))

# 6 x 22 sparse Matrix of class "dgCMatrix"
#  [[ suppressing 22 column names 'g0', 'g1', 'g2' ... ]]
#                                                                              
# MIR1302-2       .   .   .  .   .   .   1   .  1   1   .  .  .  .  .   1  .  .
# RP11-34P13.7  136 106 391 57 278  88 118 221 42 153 139 87 32 14 41 117 39 46
# RP11-34P13.8    1   2   1  .   .   .   1   .  .   1   3  .  1  .  1   .  2  .
# RP11-34P13.14   .   1   .  .   .   .   .   .  .   1   .  .  .  .  .   .  .  .
# RP11-34P13.13  78  60  69 20  78 280  31  40 11  24  43 24 10 12 10  12 13 11
# RP11-34P13.9    .   .   1  .   .   1   .   .  .   .   .  .  .  .  .   .  .  .
                      
# MIR1302-2      . . . .
# RP11-34P13.7  35 3 2 1
# RP11-34P13.8   . . . .
# RP11-34P13.14  . . . .
# RP11-34P13.13  9 . . .
# RP11-34P13.9   . . . .

################################################################################### differential expression based on
################################################################################### pseudo-bulk (cell aggregates)

# There are multiple strategies to compute the differential expression based on aggregated values.

# A common approach to simulate pseudo-replicates when we only have one sample is to randomly split the cells in each cluster 
# into several groups and then aggregate the counts for each group. 
# We may generate multiple “pseudo‐bulk” profiles per cluster, which can be used as replicates in differential expression analysis.

# In our approach, we preferred to take a simpler approach by using statistical tests such as wilcox_limma.

# Set cluster identities
Idents(pseudo) <- "seurat_clusters"  

table(Idents(pseudo)) 
# g0  g1  g2  g3  g4  g5  g6  g7  g8  g9 g10 g11 g12 g13 g14 g15 g16 g17 g18

# Find markers for Cluster 1 vs. all other clusters (pseudo-bulk expression)
cluster1_markers_pseudobulk <- FindMarkers(pseudo, 
                                ident.1 = "g1", 
                                ident.2 = NULL, 
                                only.pos = TRUE,  
                                # min.pct = 0.25, 
                                # logfc.threshold = 0.25, 
                                # test.use = "wilcox",
                                min.pct = 0.01, 
                                logfc.threshold = 0.1, 
                                # test.use = "t", 
                                # test.use = "negbinom", 
                                # test.use = "poisson", 
                                # test.use = "MAST", 
                                # test.use = "wilcox_limma", 
                                # test.use = "DESeq2", 
                                test.use = "wilcox_limma",
                                min.cells.group = 1) 


num1_significant_genes <- sum(cluster1_markers_pseudobulk$p_val < 0.1, na.rm = TRUE)
num1_significant_genes
write.table(cluster1_markers_pseudobulk, 
            file = paste0(NAME, "_Markers_Pseudobulk_pval01.Cluster1.txt"),
            sep = "\t", row.names = TRUE, quote = FALSE)
# 878 DEG

# Find markers for Cluster 3 vs. all other clusters (pseudo-bulk expression)
cluster3_markers_pseudobulk <- FindMarkers(pseudo, 
                                ident.1 = "g3", 
                                ident.2 = NULL, 
                                only.pos = TRUE,  
                                # min.pct = 0.25, 
                                # logfc.threshold = 0.25, 
                                # test.use = "wilcox",
                                min.pct = 0.01, 
                                logfc.threshold = 0.1, 
                                # test.use = "t", 
                                # test.use = "negbinom", 
                                # test.use = "poisson", 
                                # test.use = "MAST", 
                                # test.use = "wilcox_limma", 
                                # test.use = "DESeq2", 
                                test.use = "wilcox_limma",
                                min.cells.group = 1) 

num3_significant_genes <- sum(cluster3_markers_pseudobulk$p_val < 0.1, na.rm = TRUE)
num3_significant_genes
write.table(cluster3_markers_pseudobulk, 
            file = paste0(NAME, "_Markers_Pseudobulk_pval01.Cluster3.txt"),
            sep = "\t", row.names = TRUE, quote = FALSE)
# 629 DEG 

# Find markers for Cluster 4 vs. all other clusters (pseudo-bulk expression)
cluster4_markers_pseudobulk <- FindMarkers(pseudo, 
                                ident.1 = "g4", 
                                ident.2 = NULL, 
                                only.pos = TRUE,  
                                # min.pct = 0.25, 
                                # logfc.threshold = 0.25, 
                                # test.use = "wilcox",
                                min.pct = 0.01, 
                                logfc.threshold = 0.1, 
                                # test.use = "t", 
                                # test.use = "negbinom", 
                                # test.use = "poisson", 
                                # test.use = "MAST", 
                                # test.use = "wilcox_limma", 
                                # test.use = "DESeq2", 
                                test.use = "wilcox_limma",
                                min.cells.group = 1) 


num4_significant_genes <- sum(cluster4_markers_pseudobulk$p_val < 0.1, na.rm = TRUE)
num4_significant_genes
write.table(cluster4_markers_pseudobulk, 
            file = paste0(NAME, "_Markers_Pseudobulk_pval01.Cluster4.txt"),
            sep = "\t", row.names = TRUE, quote = FALSE)

# 604 DEG

###################################################################################
###################################################################################

stopCluster(cl)  # Stop cluster after completion

###################################################################################
###################################################################################