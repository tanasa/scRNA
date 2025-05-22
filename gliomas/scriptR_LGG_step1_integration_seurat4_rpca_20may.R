
# Load libraries
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

# Create output directories
plot_dir <- "./seurat4_rpca_plots"
rds_dir <- "./seurat4_rpca_rds"
dir.create(plot_dir, showWarnings = FALSE)
dir.create(rds_dir, showWarnings = FALSE)

# Define samples
base_dir <- "/home/tanasa/CNS/LGG/LGG_matrix_redoing"
samples <- list(
  "GSM5518630_LGG-04-1" = "LGG-04-1",
  "GSM5518631_LGG-04-2" = "LGG-04-2",
  "GSM5518632_LGG-04-3" = "LGG-04-3",
  "GSM5518638_LGG-03"   = "LGG-03"
)

seurat_list <- list()

# Step 1: Process each sample individually

for (prefix in names(samples)) {
  sample_id <- samples[[prefix]]
  cat("Processing sample:", sample_id, "\n")

  # File paths
  matrix_path   <- file.path(base_dir, paste0(prefix, "_matrix.mtx.gz"))
  features_path <- file.path(base_dir, paste0(prefix, "_features.tsv.gz"))
  barcodes_path <- file.path(base_dir, paste0(prefix, "_barcodes.tsv.gz"))

  # Read raw data
  counts   <- readMM(matrix_path)
  features <- read.table(features_path, header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.table(barcodes_path, header = FALSE, stringsAsFactors = FALSE)

  # Assign names
  gene_names <- make.unique(features$V2)
  rownames(counts) <- gene_names
  colnames(counts) <- paste(sample_id, barcodes$V1, sep = "_")  # Ensure unique cell names

  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_id)
  seurat_obj$sample <- sample_id
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

  # QC plot
  qc_plot <- VlnPlot(seurat_obj, 
                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                     pt.size = 0.1, ncol = 3) +
                     ggtitle(paste("QC -", sample_id))

  ggsave(filename = file.path(plot_dir, paste0(sample_id, "_QC.pdf")), 
                              plot = qc_plot, width = 10, height = 4)

  # Filtering
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)

  # Normalize and preprocess
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

  # head(seurat_obj@meta.data) 
  # UMAP plot
  umap_plot <- DimPlot(seurat_obj, 
                       reduction = "umap", 
                       label = TRUE, 
                       group.by = "seurat_clusters") +
                       ggtitle(paste("UMAP -", sample_id))
  
  ggsave(filename = file.path(plot_dir, paste0(sample_id, "_UMAP.pdf")), plot = umap_plot, width = 6, height = 5)

  # Save object
  saveRDS(seurat_obj, file = file.path(rds_dir, paste0(sample_id, ".rds")))
  seurat_list[[sample_id]] <- seurat_obj

  cat("Finished sample:", sample_id, "\n\n")
}

# Step 2: Integration

cat("Performing integration...\n")

for (i in names(seurat_list)) {
  seurat_list[[i]] <- NormalizeData(seurat_list[[i]])
  seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]], selection.method = "vst", nfeatures = 2000)
}

anchors <- FindIntegrationAnchors(object.list = seurat_list, reduction = "rpca", dims = 1:20)
seurat_integrated <- IntegrateData(anchorset = anchors, dims = 1:20)

# Step 3: Downstream on integrated object
DefaultAssay(seurat_integrated) <- "integrated"
seurat_integrated <- ScaleData(seurat_integrated)
seurat_integrated <- RunPCA(seurat_integrated)
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:20)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.5)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:20)

# Integrated UMAP
umap_integrated <- DimPlot(seurat_integrated, 
                           reduction = "umap", 
                           group.by = "sample", 
                           label = FALSE, 
                           pt.size = 0.6) +
                           ggtitle("UMAP - LGG - Integrated Samples")
ggsave(filename = file.path(plot_dir, "LGG_Integrated_UMAP.pdf"), 
                            plot = umap_integrated, 
                            width = 7, height = 6)

# Integrated UMAP with cluster labels
umap_integrated2 <- DimPlot(seurat_integrated, 
                           reduction = "umap", 
                           group.by = "seurat_clusters", 
                           label = TRUE, 
                           pt.size = 0.6) +
                   ggtitle("UMAP - LGG - Integrated Clusters")

# Save UMAP with cluster labels
ggsave(filename = file.path(plot_dir, "LGG_Integrated_UMAP_Clusters.pdf"), 
       plot = umap_integrated2, 
       width = 7, height = 6)


# Arrange side by side using patchwork
combined_umap <- umap_integrated + umap_integrated2

# Save the combined plot
ggsave(filename = file.path(plot_dir, "LGG_Combined_UMAP.pdf"), 
       plot = combined_umap, 
       width = 14, height = 6)  # wider for side-by-side

# Save integrated object
saveRDS(seurat_integrated, file = file.path(rds_dir, "LGG_seurat4_integrated.rpca.rds"))

# Save the session
save.image(file = "seurat4_integrated.rpca.RData")

cat("Integration complete. All outputs saved.\n")
