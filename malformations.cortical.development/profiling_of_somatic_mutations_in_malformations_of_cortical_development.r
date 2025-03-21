library(Seurat)
library(Matrix)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(data.table)
library(vioplot)

# Fastq files from single-nucleus libraries were processed through Cell Ranger (v6.0.2) analysis pipeline with –include-introns 
# option and hg19 reference genome. Pooled library was demultiplexed and singlets were taken by demuxlet (v1.0). 
# Seurat (v4.0.5) package was used to handle single nuclei data objects. 
# Nuclei passed a control filter (number of genes > 500, number of reads >1000, percentage of mitochondrial gene < 10%) 
# was used for downstream analysis. Protein coding genes were used for further downstream analysis. 
# Data were normalized and scaled with the most variable 3000 features using the ‘SCTransform’ functions. 
# Dimensionality reduction by PCA and UMAP embedding was performed using runPCA and runUMAP functions. 
# Clustering was performed by FindNeighbors and FindClusters functions. 
# Cell type identification was performed using known cell type markers expressed in the brain including 
# excitatory (RORB, CUX2, SATB2), inhibitory neuron (GAD1, GAD2), astrocyte (SLC1A2, SLC1A3), 
# oligodendrocyte (MOBP, PLP1), immature oligodendrocyte (BCAS1), oligodendrocyte precursor cell (PDGFRA), 
# microglia (PTPRC), and endothelial cell markers (CLDN5, ID1) as well as using positive markers found 
# by FindAllMarkers function with 3000 most variable features in scaled data. 

# DEG analysis was performed by ‘FindMarkers’ function in Seurat v4.0 with all genes available in the assay. 
# The genes with adjusted p-value < 0.01 were taken and listed in Supplementary Table 5c. 
# The final visualization of various snRNAseq data was performed by ggplot2 (v3.3.5) and matplotlib (v3.5.0).

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218022
# https://pmc.ncbi.nlm.nih.gov/articles/PMC9961399/

library(future)
options(future.globals.maxSize = 20 * 1024^3) 

barcodes <- read.table(gzfile("GSE218022_barcodes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
genes <- read.table(gzfile("GSE218022_genes.tsv.gz"), header = FALSE, stringsAsFactors = FALSE)
matrix <- readMM(gzfile("GSE218022_matrix.mtx.gz"))
metadata <- read.csv(gzfile("GSE218022_metadata.csv.gz"), header = TRUE, row.names = 1)

head(barcodes, 2)
head(genes, 2)
# head(matrix, 2)
head(metadata, 2)
colnames(metadata, 2)

head(metadata, 2)

# Convert matrix to dgCMatrix (Seurat-compatible sparse format)
matrix <- as(matrix, "CsparseMatrix")

# Assign row and column names to match metadata
rownames(matrix) <- genes$V1  
colnames(matrix) <- barcodes$V1

# Check if metadata row names match colnames of matrix
if (!all(colnames(matrix) %in% rownames(metadata))) {
  stop("Error: Metadata row names do not match cell barcodes in the matrix.")
}

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = matrix, meta.data = metadata)

# Compute mitochondrial percentage if MT genes exist
if (any(grepl("^MT-", rownames(seurat_obj)))) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
} else {
  warning("No mitochondrial genes detected with the ^MT- pattern.")
}

# Set the default assay
DefaultAssay(seurat_obj) <- "RNA"

# Print summary
# print(seurat_obj)

num_cells <- nrow(GetAssayData(seurat_obj, assay = "RNA", slot = "counts"))  # Number of cells
num_genes <- ncol(GetAssayData(seurat_obj, assay = "RNA", slot = "counts"))  # Number of genes

cat("Number of cells:", num_cells, "\n")
cat("Number of genes:", num_genes, "\n")

# Convert all three columns to numeric
metadata$nFeature_RNA <- as.numeric(metadata$nFeature_RNA)
metadata$nCount_RNA <- as.numeric(metadata$nCount_RNA)
metadata$percent_mt <- as.numeric(metadata$percent_mt)

# Print summary statistics for each column
summary(metadata[, c("nFeature_RNA", "nCount_RNA", "percent_mt")])

# Check if metadata columns are present
required_features <- c("nFeature_RNA", "nCount_RNA", "percent_mt", "seurat_clusters", "cell_type")
missing_features <- setdiff(required_features, colnames(seurat_obj@meta.data))

if (length(missing_features) > 0) {
  stop(paste("Error: Missing required metadata columns:", paste(missing_features, collapse = ", ")))
}

df = metadata

# Ensure seurat_clusters is a factor
df$seurat_clusters <- as.factor(df$seurat_clusters)
# Define colors for each cluster
cluster_colors <- rainbow(length(levels(df$seurat_clusters)))

# Convert to numeric safely, handling non-numeric values
df$nFeature_RNA <- as.numeric(as.character(df$nFeature_RNA))
df$nCount_RNA <- as.numeric(as.character(df$nCount_RNA))
df$percent_mt <- as.numeric(as.character(df$percent_mt))

# Remove NA values to prevent plotting errors
df <- df[!is.na(df$nFeature_RNA) & !is.na(df$nCount_RNA) & !is.na(df$percent_mt), ]

# Check if there are empty clusters after filtering
valid_clusters <- levels(df$seurat_clusters)[table(df$seurat_clusters) > 0]
print(valid_clusters)


# Set a larger figure size for Jupyter Notebook
options(repr.plot.width = 16, repr.plot.height = 8)  # Increase width & height

# Ensure seurat_clusters is a factor
df$seurat_clusters <- as.factor(df$seurat_clusters)

# Convert to numeric safely
df$nFeature_RNA <- as.numeric(as.character(df$nFeature_RNA))

# Remove NA values
df <- df[!is.na(df$nFeature_RNA), ]

# Generate violin + jitter plot
ggplot(df, aes(x = seurat_clusters, y = nFeature_RNA, fill = seurat_clusters)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_jitter(width = 0.2, size = 1.8, color = "black", alpha = 0.5) +  # Jitter points
  theme_minimal(base_size = 20) +  # Larger base font size
  labs(title = "nFeature_RNA by Seurat Cluster", x = "Seurat Cluster", y = "nFeature_RNA") +
  theme(legend.position = "none",  # Remove legend
        plot.title = element_text(hjust = 0.5, size = 24),  # Center & enlarge title
        axis.text.x = element_text(size = 18),  # Increase x-axis text size
        axis.text.y = element_text(size = 18))  # Increase y-axis text size


# Set a larger figure size for Jupyter Notebook
options(repr.plot.width = 16, repr.plot.height = 8)  # Increase width & height

# Ensure seurat_clusters is a factor
df$seurat_clusters <- as.factor(df$seurat_clusters)

# Convert to numeric safely
df$nCount_RNA <- as.numeric(as.character(df$nCount_RNA))

# Remove NA values
df <- df[!is.na(df$nCount_RNA), ]

# Generate violin + jitter plot for nCount_RNA with y-limit
ggplot(df, aes(x = seurat_clusters, y = nCount_RNA, fill = seurat_clusters)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black") +  # Violin plot with transparency and black outline
  geom_jitter(width = 0.2, size = 1.8, color = "black", alpha = 0.5) +  # Jitter points
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  ylim(0, 100000) +  # Set y-axis limit
  theme_minimal(base_size = 20) +  # Larger base font size
  labs(title = "nCount_RNA by Seurat Cluster", 
       x = "Seurat Cluster", 
       y = "nCount_RNA") +
  theme(legend.position = "none",  # Remove legend
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),  # Center & bold title
        axis.text.x = element_text(size = 18),  # Increase x-axis text size
        axis.text.y = element_text(size = 18))  # Increase y-axis text size

# Set a larger figure size for Jupyter Notebook
options(repr.plot.width = 16, repr.plot.height = 8)  # Increase width & height

# Ensure seurat_clusters is a factor
df$seurat_clusters <- as.factor(df$seurat_clusters)

# Convert to numeric safely
df$percent_mt <- as.numeric(as.character(df$percent_mt))

# Remove NA values
df <- df[!is.na(df$percent_mt), ]

# Generate violin + jitter plot for percent_mt with y-limit
ggplot(df, aes(x = seurat_clusters, y = percent_mt, fill = seurat_clusters)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = "black") +  # Violin plot with transparency and black outline
  geom_jitter(width = 0.2, size = 1.8, color = "black", alpha = 0.5) +  # Jitter points
  scale_fill_brewer(palette = "Set2") +  # Nice color palette
  ylim(0, 15) +  # Set y-axis limit to 100%
  theme_minimal(base_size = 20) +  # Larger base font size
  labs(title = "percent_mt by Seurat Cluster", 
       x = "Seurat Cluster", 
       y = "percent_mt") +
  theme(legend.position = "none",  # Remove legend
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),  # Center & bold title
        axis.text.x = element_text(size = 18),  # Increase x-axis text size
        axis.text.y = element_text(size = 18))  # Increase y-axis text size


# Generate scatter plot
ggplot(df, aes(x = nCount_RNA, y = nFeature_RNA, color = seurat_clusters)) +
  geom_point(alpha = 0.6, size = 2) +  # Scatter plot points with transparency
  scale_color_brewer(palette = "Dark2") +  # Color scheme for clusters
  theme_minimal(base_size = 18) +  # Bigger font sizes
  labs(title = "Scatter Plot: nCount_RNA vs. nFeature_RNA", 
       x = "nCount_RNA", 
       y = "nFeature_RNA", 
       color = "Seurat Cluster") +  # Color legend
  xlim(0, 30000) +  # Custom x-axis limits
  ylim(0, 8000) +  # Custom y-axis limits
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),  # Center & bold title
        axis.text.x = element_text(size = 16),  # Larger x-axis labels
        axis.text.y = element_text(size = 16))  # Larger y-axis labels

# seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 &
#                                  nCount_RNA > 200 & nCount_RNA < 10000 &
#                                  percent_mt < 10)

# num_cells <- ncol(GetAssayData(seurat_obj, assay = "RNA", slot = "counts"))  # Number of cells
# num_genes <- nrow(GetAssayData(seurat_obj, assay = "RNA", slot = "counts"))  # Number of genes

# cat("Number of cells:", num_cells, "\n")
# cat("Number of genes:", num_genes, "\n")
# head(metadata)

# head(GetAssayData(seurat_obj, assay = "RNA", slot = "counts"), 3)

seurat_obj@meta.data$cell_type <- as.factor(seurat_obj@meta.data$cell_type_detail)
unique(seurat_obj@meta.data$cell_type)
length(unique(seurat_obj@meta.data$cell_type))

check_seurat_data <- function(seurat_obj) {
  found_data <- FALSE  # Flag to track if any data is found
  
  for (assay in Assays(seurat_obj)) {
    slots <- slotNames(seurat_obj[[assay]])
    
    if ("data" %in% slots) {
      message("✅ Normalized data found in: ", assay)
      found_data <- TRUE
    }
    
    if ("scale.data" %in% slots) {
      message("✅ Standardized (scaled) data found in: ", assay)
      found_data <- TRUE
    }
  }
  
  if (!found_data) {
    message("❌ No normalized or standardized data found in any assay.")
  }
}

check_seurat_data(seurat_obj)

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale the data
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

# Alternatively, we can use SCTransform. 

# Normalize & standardize data using SCTransform (with top 3000 variable features)
# seurat_obj <- SCTransform(seurat_obj, variable.features.n = 3000, verbose = FALSE)

print("variable features") 

variable_features_df <- HVFInfo(seurat_obj)
variable_features_df$gene <- rownames(variable_features_df)  # Add gene names as a column

# Extract top 10 variable genes
top10_genes <- head(VariableFeatures(seurat_obj), 10)
print(top10_genes)
      
# Generate base variable feature plot
plot <- VariableFeaturePlot(seurat_obj)

# Add labels for the top 10 most variable features
plot + 
  geom_text_repel(data = subset(variable_features_df, gene %in% top10_genes), 
                  aes(x = mean, y = variance.standardized, label = gene), 
                  size = 5, color = "red", max.overlaps = Inf) +  # Use ggrepel to avoid overlap
  theme_minimal(base_size = 16) +
  labs(title = "Top 10 Highly Variable Features in Seurat Object") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))

# save.image(file = "estephany1.RData")

head(metadata,3)
metadata[1,1]
unique(seurat_obj@meta.data$cell_type_detail)
length(unique(metadata$cell_type_detail))



# Perform PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# Plot PCA using DimPlot
DimPlot(seurat_obj, reduction = "pca") +
  theme_minimal(base_size = 14) +
  labs(title = "PCA Plot of Seurat Object") +
  xlim(-20, 20) +  
  ylim(-20, 20) +  
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))

DimHeatmap(seurat_obj, dims = 1:5, cells = 500, balanced = TRUE)
ElbowPlot(seurat_obj, ndims = 50)  

# Find neighbors and clusters
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)

# Run UMAP for visualization
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# Visualize the clustering
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

# save.image(file = "estephany3.RData")
# load("estephany3.RData")

options(repr.plot.width = 8, repr.plot.height = 6)
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

print("UMAP visualization of pre-annotated cell clusters")
DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type_detail")

length(unique(seurat_obj@meta.data$cell_type_detail))

options(repr.plot.width = 12, repr.plot.height = 10)
print("excitatory neurons") 
FeaturePlot(seurat_obj, features = c("RORB", "CUX2", "SATB2"))  

options(repr.plot.width = 12, repr.plot.height = 6)
print("inhibitory neurons")
FeaturePlot(seurat_obj, features = c("GAD1", "GAD2"))  

options(repr.plot.width = 12, repr.plot.height = 6)
print("astrocytes") 
FeaturePlot(seurat_obj, features = c("SLC1A2", "SLC1A3"))  

options(repr.plot.width = 12, repr.plot.height = 6)
print("oligodendrocyte") 
FeaturePlot(seurat_obj, features = c("MOBP", "PLP1"))  

options(repr.plot.width = 6, repr.plot.height = 6)
print("immature oligodendrocyte") 
FeaturePlot(seurat_obj, features = c("BCAS1"))  

options(repr.plot.width = 6, repr.plot.height = 6)
print("oligodendrocyte precursor cell") 
FeaturePlot(seurat_obj, features = c("PDGFRA"))  

options(repr.plot.width = 6, repr.plot.height = 6)
print("microglia") 
FeaturePlot(seurat_obj, features = c("PTPRC"))  

options(repr.plot.width = 12, repr.plot.height = 6)
print("endothelial cells") 
FeaturePlot(seurat_obj, features = c("CLDN5", "ID1"))  

DotPlot(seurat_obj, features = list(
  Excitatory = c("RORB", "CUX2", "SATB2"),
  Inhibitory = c("GAD1", "GAD2"),
  Astrocyte = c("SLC1A2", "SLC1A3"),
  Oligodendrocyte = c("MOBP", "PLP1"),
  Immature_Oligo = c("BCAS1"),
  OPC = c("PDGFRA"),
  Microglia = c("PTPRC"),
  Endothelial = c("CLDN5", "ID1")
)) + RotatedAxis()

# Find markers for each cluster
cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster_markers, file = "cluster_markers.csv", row.names = TRUE)


