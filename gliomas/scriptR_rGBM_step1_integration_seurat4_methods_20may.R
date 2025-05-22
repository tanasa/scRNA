
# Load libraries
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(future)
options(future.globals.maxSize = 1000 * 1024^3)  # 5 GB
plan(multisession, workers = 10)

# Create output directories
plot_dir <- "./seurat4_methods_plots"
rds_dir <- "./seurat4_methods_rds"
dir.create(plot_dir, showWarnings = FALSE)
dir.create(rds_dir, showWarnings = FALSE)

# Define samples
base_dir <- "/home/tanasa/CNS/rGBM"
samples <- list(
"GSM5518596_rGBM-01-A" = "rGBM-01-A",
"GSM5518597_rGBM-01-B" = "rGBM-01-B",
"GSM5518598_rGBM-01-C" = "rGBM-01-C",
"GSM5518599_rGBM-01-D" = "rGBM-01-D",
"GSM5518612_rGBM-02-2" = "rGBM-02-2",
"GSM5518613_rGBM-02-3" = "rGBM-02-3",
"GSM5518614_rGBM-02-4" = "rGBM-02-4",
"GSM5518615_rGBM-02-5" = "rGBM-02-5",
"GSM5518616_rGBM-03-1" = "rGBM-03-1",
"GSM5518617_rGBM-03-2" = "rGBM-03-2",
"GSM5518618_rGBM-03-3" = "rGBM-03-3",
"GSM5518619_rGBM-04-1" = "rGBM-04-1",
"GSM5518620_rGBM-04-2" = "rGBM-04-2",
"GSM5518621_rGBM-04-3" = "rGBM-04-3",
"GSM5518622_rGBM-04-4" = "rGBM-04-4",
"GSM5518626_rGBM-05-1" = "rGBM-05-1",
"GSM5518627_rGBM-05-2" = "rGBM-05-2",
"GSM5518628_rGBM-05-3" = "rGBM-05-3"
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

# Step 2: Integration with multiple methods
reduction_methods <- c("rpca", "jpca", "rlsi", "cca")

for (method in reduction_methods) {
  cat("\n=============================\n")
  cat("Integrating using:", method, "\n")
  cat("=============================\n")

  # Re-normalize and find HVGs
  for (i in names(seurat_list)) {
    seurat_list[[i]] <- NormalizeData(seurat_list[[i]], verbose = FALSE)
    seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }

  # Integration Anchors
  anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                    reduction = method,
                                    dims = 1:20)

  # Integrate Data
  seurat_integrated <- IntegrateData(anchorset = anchors, dims = 1:20)

  # Downstream processing
  DefaultAssay(seurat_integrated) <- "integrated"
  seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE)
  seurat_integrated <- RunPCA(seurat_integrated, verbose = FALSE)
  seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:20, verbose = FALSE)
  seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.5, verbose = FALSE)
  seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:20, verbose = FALSE)

  # Plotting
  umap_by_sample <- DimPlot(seurat_integrated, reduction = "umap", group.by = "sample", label = FALSE, pt.size = 0.6) +
                    ggtitle(paste("UMAP -", toupper(method), "- Integrated Samples"))

  umap_by_cluster <- DimPlot(seurat_integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 0.6) +
                     ggtitle(paste("UMAP -", toupper(method), "- Integrated Clusters"))

  # Save plots
  ggsave(filename = file.path(plot_dir, paste0("LGG_", method, "_Integrated_UMAP.pdf")),
         plot = umap_by_sample, width = 7, height = 6)

  ggsave(filename = file.path(plot_dir, paste0("LGG_", method, "_Integrated_UMAP_Clusters.pdf")),
         plot = umap_by_cluster, width = 7, height = 6)

  ggsave(filename = file.path(plot_dir, paste0("LGG_", method, "_Combined_UMAP.pdf")),
         plot = umap_by_sample + umap_by_cluster, width = 14, height = 6)

  # Save integrated object
  saveRDS(seurat_integrated, file = file.path(rds_dir, paste0("LGG_seurat4_integrated_", method, ".rds")))
}

cat("All integration methods completed.\n")

# Save the session
save.image(file = "LGG_seurat4_integrated_methods.all.RData")

cat("Integration complete. All outputs saved.\n")
