# Article : https://pubmed.ncbi.nlm.nih.gov/38211712/
# Data were downloaded from :
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE209832

library(Seurat)
library(Matrix)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(data.table)
library(vioplot)
library(harmony)
library(cowplot)
library(patchwork)
library(future)
library(foreach)
library(doParallel)
options(future.globals.maxSize = 20 * 1024^3) 
packageVersion("Seurat")



# CD1 : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6401749
#       GSM6401749_GEX_RST10031plus2.filtered_feature_bc_matrix.h5

cd1 = Read10X_h5("GSM6401749_GEX_RST10031plus2.filtered_feature_bc_matrix.h5")
cd1 = CreateSeuratObject(counts = cd1, project = "cd1")

# CD2 : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6401751
#       GSM6401751_GEX_RST10147.filtered_feature_bc_matrix.h5

cd2 = Read10X_h5("GSM6401751_GEX_RST10147.filtered_feature_bc_matrix.h5")
cd2 = CreateSeuratObject(counts = cd2, project = "cd2")

# CD3 : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6401752
#       GSM6401752_GEX_RST10148.filtered_feature_bc_matrix.h5

cd3 = Read10X_h5("GSM6401752_GEX_RST10148.filtered_feature_bc_matrix.h5")
cd3 = CreateSeuratObject(counts = cd3, project = "cd3")

# CD4 : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6401753
#       GSM6401753_GEX_RST10150.filtered_feature_bc_matrix.h5

cd4 = Read10X_h5("GSM6401753_GEX_RST10150.filtered_feature_bc_matrix.h5")
cd4 = CreateSeuratObject(counts = cd4, project = "cd4")

# CD5  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6401756
#      GSM6401756_GEX_RST10175.filtered_feature_bc_matrix.h5

cd5 = Read10X_h5("GSM6401756_GEX_RST10175.filtered_feature_bc_matrix.h5")
cd5 = CreateSeuratObject(counts = cd5, project = "cd5")

dim(cd1@meta.data)[1]
dim(cd2@meta.data)[1]
dim(cd3@meta.data)[1]
dim(cd4@meta.data)[1]
dim(cd5@meta.data)[1]

# Healthy 1 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6401750
#            GSM6401750_GEX_RST10142.filtered_feature_bc_matrix.h5

hd1 = Read10X_h5("GSM6401750_GEX_RST10142.filtered_feature_bc_matrix.h5")
hd1 = CreateSeuratObject(counts = hd1, project = "hd1")

# Healthy 2 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6401754
#           GSM6401754_GEX_RST10153.filtered_feature_bc_matrix.h5

hd2 = Read10X_h5("GSM6401754_GEX_RST10153.filtered_feature_bc_matrix.h5")
hd2 = CreateSeuratObject(counts = hd2, project = "hd2")

# Healthy 3 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6401755
#            GSM6401755_GEX_RST10168.filtered_feature_bc_matrix.h5

hd3 = Read10X_h5("GSM6401755_GEX_RST10168.filtered_feature_bc_matrix.h5")
hd3 = CreateSeuratObject(counts = hd3, project = "hd3")

# Healthy 4 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6401757 
#            GSM6401757_GEX_RST12171.filtered_feature_bc_matrix.h5

hd4 = Read10X_h5("GSM6401757_GEX_RST12171.filtered_feature_bc_matrix.h5")
hd4 = CreateSeuratObject(counts = hd4, project = "hd4")

# Healthy 5 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6401758
#             GSM6401758_GEX_RST12173.filtered_feature_bc_matrix.h5

hd5 = Read10X_h5("GSM6401758_GEX_RST12173.filtered_feature_bc_matrix.h5")
hd5 = CreateSeuratObject(counts = hd5, project = "hd5")

# Healthy 6 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6401759
#            GSM6401759_GEX_RST12179.filtered_feature_bc_matrix.h5

hd6 = Read10X_h5("GSM6401759_GEX_RST12179.filtered_feature_bc_matrix.h5")
hd6 = CreateSeuratObject(counts = hd6, project = "hd6")

dim(hd1@meta.data)[1]
dim(hd2@meta.data)[1]
dim(hd3@meta.data)[1]
dim(hd4@meta.data)[1]
dim(hd5@meta.data)[1]

# healthy donors :  Merge all Seurat objects
healthy <- merge(hd1, y = list(hd2, hd3, hd4, hd5, hd6), 
                          add.cell.ids = c("hd1", "hd2", "hd3", "hd4", "hd5", "hd6"))

healthy@meta.data$status <- "hd"
# unique(healthy@meta.data$orig.ident)
# head(healthy@meta.data, 1)
print("number of cells : hd")
dim(healthy@meta.data)[1]

# CD donors :  Merge all Seurat objects
disease <- merge(cd1, y = list(cd2, cd3, cd4, cd5), 
                          add.cell.ids = c("cd1", "cd2", "cd3", "cd4", "cd5"))

disease@meta.data$status <- "cd"
# unique(disease@meta.data$orig.ident)
# head(disease@meta.data, 1)
print("number of cells : cd")
dim(disease@meta.data)[1]

# to free the memory 

rm(hd1)
rm(hd2)
rm(hd3)
rm(hd4)
rm(hd5)
rm(hd6)

rm(cd1)
rm(cd2)
rm(cd3)
rm(cd4)
rm(cd5)

combined <- merge(healthy, y = disease, add.cell.ids = c("healthy", "disease"))

print("the number of cells in each sample:")
table(combined@meta.data$orig.ident)
# head(combined@meta.data, 1)
# tail(combined@meta.data, 1)
print("the number of cells in each condition:")
table(combined@meta.data$status)

# Search for mitochondrial genes with different capitalization patterns
print("mitochondrial genes")
mt_genes <- grep("^MT-|^mt-|^Mt-", rownames(combined), value = TRUE)
print(mt_genes)

# to free the memory 

rm(healthy)
rm(disease)

print("the number of cells in each condition:")
table(combined@meta.data$status)

# compute the percent of mithochondrial genes
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-|^MT-|^Mt") 

options(repr.plot.width = 20, repr.plot.height = 8)
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

options(repr.plot.width = 20, repr.plot.height = 8)
CombinePlots(plots = list(plot1, plot2))

# Compute summary statistics
summary(combined@meta.data)

RESOLUTION = 0.3

combined <- NormalizeData(combined) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)

combined <- RunHarmony(combined, group.by.vars = "orig.ident")

combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30)
combined <- RunTSNE(combined, reduction = "harmony", dims = 1:30)

Reductions(combined)

RESOLUTION = 0.3

combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = RESOLUTION)

# save.image("cd.v1.RData")

head(combined@meta.data,3)
print("number of clusters")
table(combined@meta.data$seurat_clusters)

print("Visualization using PCA :")
options(repr.plot.width = 16, repr.plot.height = 8)
p1 <- DimPlot(combined, reduction = "pca", group.by = "orig.ident")
p2 <- DimPlot(combined, reduction = "pca", group.by = c("seurat_clusters"))
plot_grid(p1, p2)

# print("Visualization using HARMONY :")
# options(repr.plot.width = 22, repr.plot.height = 8)
# p1 <- DimPlot(combined, reduction = "harmony", group.by = "orig.ident")
# p2 <- DimPlot(combined, reduction = "harmony", group.by = "seurat_clusters")
# plot_grid(p1, p2)

# Visualization using tSNE :

options(repr.plot.width = 20, repr.plot.height = 10)
p1 <- DimPlot(combined, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(combined, reduction = "tsne", group.by = "seurat_clusters")
plot_grid(p1, p2)

# Visualization using UMAP :

options(repr.plot.width = 20, repr.plot.height = 10)
p1 <- DimPlot(combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label=TRUE)
plot_grid(p1, p2)

print("T-reg specific markers")

VlnPlot(combined, features = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT", "LRRC32"), pt.size = 1)

fp <- FeaturePlot(
  combined,
  features = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT", "LRRC32"),
  reduction = "umap",
  pt.size = 1,
  order = TRUE,
  min.cutoff = "q10",
  max.cutoff = "q90"
)

fp[[1]] + fp[[2]] + fp[[3]] +
fp[[4]] + fp[[5]] + fp[[6]] +
          plot_layout(ncol = 3)

options(repr.plot.width = 20, repr.plot.height = 10)
p1 <- DimPlot(combined, reduction = "tsne", group.by = "seurat_clusters", label = TRUE)
p2 <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 6)
plot_grid(p1, p2)

head(combined@meta.data, 2)



print("Based on the expression of markers such as FOXP3, IL2RA, CTLA4, IKZF2, TIGIT, and LRRC32, we infer that cluster 4 likely 
      represents regulatory T cells (Tregs). However, a more comprehensive analysis should include automated cell type annotation 
      methods such as SingleR, CellTypist, Azimuth, or scmap for greater confidence.")

print("Computing the differential expression bewteen HD and CD only for cluster 4")

# selecting Cluster 4
Idents(combined) <- "seurat_clusters"

cluster4 <- WhichCells(combined, idents = "4")
cluster4o <- subset(combined, cells = cluster4)
Idents(cluster4o) <- "status"

cluster4o <- JoinLayers(cluster4o)
# ident.1 = "cd" is the test group
# ident.2 = "hd" is the reference group
cluster4_markers <- FindMarkers(
  cluster4o,
  ident.1 = "cd",
  ident.2 = "hd",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

head(cluster4_markers, 4)

# Define significance threshold and subset the significant genes

sig <- 0.05
logfc <- 0.25

sig_genes <- subset(cluster4_markers, p_val_adj < sig & abs(avg_log2FC) > logfc)
up <- sum(sig_genes$avg_log2FC > logfc)
down <- sum(sig_genes$avg_log2FC < -logfc)

cat("Significantly upregulated genes in cd vs hd:", up, "\n")
cat("Significantly downregulated genes in cd vs hd:", down, "\n")

# cat("top 10 up-regulated genes")
head(sig_genes[order(-sig_genes$avg_log2FC), ], 5)

# cat("top 10 down-regulated genes")
head(sig_genes[order(sig_genes$avg_log2FC), ], 5)

# re-order the table print the results

sig_genes$gene <- rownames(sig_genes)
head(sig_genes)

sig_genes$direction <- ifelse(sig_genes$avg_log2FC > 0, "up sj", "down sj")
sig_table <- sig_genes[, c("gene", "avg_log2FC", "p_val_adj", "pct.1", "pct.2", "direction")]
sig_table <- sig_table[order(sig_table$avg_log2FC, decreasing = TRUE), ]

write.csv(sig_table, "cluster4_DE_genes_cd_vs_hd.csv", row.names = FALSE)

print("The number of differentially expressed genes in cluster 4 bewteen Healthy Donors and Crohn Patients is :")
dim(sig_genes)

print(paste("Over-representation and Gene Set Enrichment analysis for", nrow(sig_genes), "genes"))

print("The analysis with the package clusterprofiler")



library(clusterProfiler)
library(org.Hs.eg.db) 
library(GO.db)         
library(DO.db)         
library(KEGGREST)      
library(ReactomePA)    
library(enrichplot)    
library(dplyr)
library(msigdbr)
library(msigdb)
library(msigdf)
library(msigdbdf)

gene_list = sig_genes$gene
# head(gene_list)
filtered = sig_genes
# head(filtered)

# Convert gene names to Entrez IDs
gene_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# head(gene_ids, 2)
gene_list_merged <- merge(gene_ids, filtered, by.x = "SYMBOL", by.y = "gene")
# head(gene_list_merged, 2)
fin_name = "cluster4_DE_genes_cd_vs_hd.clusterprofiler"



gene_list2 <- setNames(gene_list_merged$avg_log2FC, gene_list_merged$ENTREZID)
gene_list2 <- sort(gene_list2, decreasing = TRUE)
head(gene_list2, 2)

# GO Over-Representation Analysis

result <- tryCatch({

  ego <- enrichGO(gene          = gene_ids$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  keyType       = "ENTREZID",  
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 1,
                  readable      = TRUE)
  
  # Save results (using the result slot for consistency)
  write.table(ego@result, file = paste0(fin_name, "_GO_OverRepresentation_Results.txt"), 
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  png(paste0(fin_name, "_GO_OverRepresentation.png"), width = 1000, height = 800)
  print(dotplot(ego, showCategory = 20))
  dev.off()
  
  # Display plot interactively
  # print(dotplot(ego, showCategory = 20))
  
}, error = function(e) {
  cat("An error occurred in GO over-representation analysis:", conditionMessage(e), "\n")
})

options(repr.plot.width = 8, repr.plot.height = 8)
print(dotplot(ego, showCategory = 20))

# GO Enrichment Analysis (GSEA)

result2 <- tryCatch({

  ego2 <- gseGO(gene = gene_list2,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",  
                ont = "BP",      # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)
  
  write.table(ego2@result, file = paste0(fin_name, "_GO_Enrichment_Results.txt"), 
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  png(paste0(fin_name, "_GO_Enrichment_Plot.png"), width = 1000, height = 800)
  print(dotplot(ego2, showCategory = 20))
  dev.off()
  
  # print(dotplot(ego2, showCategory = 20))
  
}, error = function(e) {
  cat("An error occurred in GO enrichment analysis:", conditionMessage(e), "\n")
})

options(repr.plot.width = 8, repr.plot.height = 8)
print(dotplot(ego2, showCategory = 20))

# KEGG Over-Representation Analysis

result <- tryCatch({

  kegg_enrich <- enrichKEGG(gene = gene_ids$ENTREZID,
                            organism = "hsa",  
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05)
    
  if (nrow(kegg_enrich@result) > 0) {
    # print(dotplot(kegg_enrich, showCategory = 20))
    write.table(kegg_enrich@result, 
                file = paste0(fin_name, "_KEGG_OverRepresentation_Results.txt"), 
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    png(paste0(fin_name, "_KEGG_OverRepresentation_Plot.png"), width = 600, height = 600)
    print(dotplot(kegg_enrich, showCategory = 20))
    dev.off()
    
  } else {
    cat("No enriched KEGG terms found. Please adjust your parameters.\n")
  }
  
}, error = function(e) {
  cat("An error occurred during KEGG over-representation analysis:", conditionMessage(e), "\n")
})

options(repr.plot.width = 8, repr.plot.height = 8)
print(dotplot(kegg_enrich, showCategory = 20))

# KEGG Gene Set Enrichment Analysis

result2 <- tryCatch({
    
  kegg_gse <- gseKEGG(geneList   = gene_list2,
                      organism     = 'hsa',
                      minGSSize    = 10,
                      pvalueCutoff = 0.05,
                      verbose      = FALSE)
  
  if (nrow(kegg_gse@result) > 0) {
    # print(dotplot(kegg_gse, showCategory = 20))
    write.table(kegg_gse@result, "KEGG_Enrichment_Results.txt", sep = "\t", row.names = FALSE, quote = FALSE)
    
    png(paste0(fin_name, "_KEGG_Enrichment_Plot.png"), width = 600, height = 600)
    print(dotplot(kegg_gse, showCategory = 20))
    dev.off()
    
  } else {
    cat("No enriched KEGG terms found. Please adjust your parameters.\n")
  }
  
}, error = function(e) {
  cat("An error occurred during KEGG GSEA analysis:", conditionMessage(e), "\n")
})

options(repr.plot.width = 8, repr.plot.height = 8)
print(dotplot(kegg_gse, showCategory = 20))

################################################################################ ORA
# WikiPathways Over-Representation Analysis

result <- tryCatch({

  wikipathways_enrich <- enrichWP(
    gene = gene_ids$ENTREZID,
    organism = "Homo sapiens",
    pvalueCutoff = 0.05
  )
  
  if (nrow(wikipathways_enrich@result) > 0) {
    write.table(wikipathways_enrich@result,
                file = paste0(fin_name, "_WikiPathways_ORA_Results.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    png(paste0(fin_name, "_WikiPathways_ORA_Plot.png"), width = 1000, height = 800)
    print(dotplot(wikipathways_enrich, showCategory = 20))
    dev.off()
  } else {
    cat("No enriched WikiPathways terms found. Please adjust your parameters.\n")
  }
  
}, error = function(e) {
  cat("An error occurred during WikiPathways over-representation analysis:", conditionMessage(e), "\n")
})

options(repr.plot.width = 8, repr.plot.height = 8)
print(dotplot(wikipathways_enrich, showCategory = 20))

# WikiPathways Gene Set Enrichment Analysis

result2 <- tryCatch({

  wikipathways_gse <- gseWP(
    gene = gene_list2,
    organism = "Homo sapiens",
    pvalueCutoff = 0.05
  )
  
  if (nrow(wikipathways_gse@result) > 0) {
    write.table(wikipathways_gse@result,
                file = paste0(fin_name, "_WikiPathways_Enrichment_Results.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    png(paste0(fin_name, "_WikiPathways_Enrichment_Plot.png"), width = 1000, height = 800)
    print(dotplot(wikipathways_gse, showCategory = 20))
    dev.off()
  } else {
    cat("No enriched WikiPathways terms found. Please adjust your parameters.\n")
  }
  
}, error = function(e) {
  cat("An error occurred during WikiPathways GSEA analysis:", conditionMessage(e), "\n")
})

options(repr.plot.width = 8, repr.plot.height = 8)
print(dotplot(wikipathways_gse, showCategory = 20))

# WikiPathways Gene Set Enrichment Analysis

result2 <- tryCatch({

  wikipathways_gse <- gseWP(
    gene = gene_list2,
    organism = "Homo sapiens",
    pvalueCutoff = 0.05
  )
  
  if (nrow(wikipathways_gse@result) > 0) {
    write.table(wikipathways_gse@result,
                file = paste0(fin_name, "_WikiPathways_Enrichment_Results.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    png(paste0(fin_name, "_WikiPathways_Enrichment_Plot.png"), width = 1000, height = 800)
    print(dotplot(wikipathways_gse, showCategory = 20))
    dev.off()
  } else {
    cat("No enriched WikiPathways terms found. Please adjust your parameters.\n")
  }
  
}, error = function(e) {
  cat("An error occurred during WikiPathways GSEA analysis:", conditionMessage(e), "\n")
})

options(repr.plot.width = 8, repr.plot.height = 8)
print(dotplot(wikipathways_gse, showCategory = 20))

################################################################################ GSEA
# Reactome Gene Set Enrichment Analysis using fgsea

result2 <- tryCatch({

  reactome_gse <- gsePathway(
    gene = gene_list2,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  if (nrow(reactome_gse@result) > 0) {
    write.table(reactome_gse@result,
                file = paste0(fin_name, "_Reactome_Enrichment_Results.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    png(paste0(fin_name, "_Reactome_Enrichment_Plot.png"), width = 1000, height = 800)
    print(dotplot(reactome_gse, showCategory = 20))
    dev.off()
  } else {
    cat("No enriched Reactome terms found. Please adjust your parameters.\n")
  }
  
}, error = function(e) {
  cat("An error occurred during Reactome GSEA analysis:", conditionMessage(e), "\n")
})

options(repr.plot.width = 6, repr.plot.height = 6)
print(dotplot(reactome_gse, showCategory = 20))

# save.image("cd.v2.RData")


