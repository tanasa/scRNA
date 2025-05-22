# Load libraries
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(batchelor)
library(future)
library(future.apply)

options(future.globals.maxSize = 128 * 1024^3)  # 128 GB
plan("multicore")  # Use "multisession" on Windows

# Define samples
base_dir <- "/home/tanasa/CNS/LGG/LGG_matrix_redoing"
samples <- list(
  "GSM5518630_LGG-04-1" = "LGG-04-1",
  "GSM5518631_LGG-04-2" = "LGG-04-2",
  "GSM5518632_LGG-04-3" = "LGG-04-3",
  "GSM5518638_LGG-03"   = "LGG-03"
)

# Output directories
plot_dir <- "./seurat5_plots"
rds_dir <- "./seurat5_rds"
dir.create(plot_dir, showWarnings = FALSE)
dir.create(rds_dir, showWarnings = FALSE)

# Step 1: Individual preprocessing
seurat_list <- future_lapply(names(samples), function(prefix) {
  sample_id <- samples[[prefix]]
  cat("Processing sample:", sample_id, "\n")

  # Load raw data
  matrix_path   <- file.path(base_dir, paste0(prefix, "_matrix.mtx.gz"))
  features_path <- file.path(base_dir, paste0(prefix, "_features.tsv.gz"))
  barcodes_path <- file.path(base_dir, paste0(prefix, "_barcodes.tsv.gz"))

  counts   <- readMM(matrix_path)
  features <- read.table(features_path, header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.table(barcodes_path, header = FALSE, stringsAsFactors = FALSE)

  gene_names <- make.unique(features$V2)
  rownames(counts) <- gene_names
  colnames(counts) <- paste(sample_id, barcodes$V1, sep = "_")

  seurat_obj <- CreateSeuratObject(counts = counts, project = sample_id)
  seurat_obj$sample <- sample_id
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

  qc_plot <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3) +
             ggtitle(paste("QC -", sample_id))
  ggsave(file.path(plot_dir, paste0(sample_id, "_QC.pdf")), plot = qc_plot, width = 10, height = 4)

  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

  umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, group.by = "seurat_clusters") +
               ggtitle(paste("UMAP -", sample_id))
  ggsave(file.path(plot_dir, paste0(sample_id, "_UMAP.pdf")), plot = umap_plot, width = 6, height = 5)

  saveRDS(seurat_obj, file = file.path(rds_dir, paste0(sample_id, ".rds")))
  cat("Finished sample:", sample_id, "\n\n")
  return(seurat_obj)
})
names(seurat_list) <- samples

cat("Preparing merged object for integration...\n")
obj_merged <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
obj_merged <- NormalizeData(obj_merged)
obj_merged <- FindVariableFeatures(obj_merged, selection.method = "vst", nfeatures = 2000)
obj_merged <- ScaleData(obj_merged)
obj_merged <- RunPCA(obj_merged)

integration_methods <- list(
  RPCA = list(method = RPCAIntegration, new_reduction = "integrated.rpca", cluster_name = "rpca_clusters", umap_name = "umap.rpca"),
  CCA = list(method = CCAIntegration, new_reduction = "integrated.cca", cluster_name = "cca_clusters", umap_name = "umap.cca"),
  Harmony = list(method = HarmonyIntegration, new_reduction = "harmony", cluster_name = "harmony_clusters", umap_name = "umap.harmony")
)

for (method_name in names(integration_methods)) {
  cat("Running", method_name, "integration...\n")
  obj <- obj_merged
  method_info <- integration_methods[[method_name]]

  if (method_name %in% c("RPCA", "CCA")) {
    obj <- IntegrateLayers(
      object = obj,
      method = method_info$method,
      orig.reduction = "pca",
      new.reduction = method_info$new_reduction,
      verbose = TRUE
    )
  } else if (tolower(method_name) %in% c("harmony")) {
    obj <- IntegrateLayers(
      object = obj,
      method = method_info$method,
      new.reduction = method_info$new_reduction,
      verbose = TRUE
    )
  } else {
    stop(paste("Unsupported integration method:", method_name))
  }

  obj <- FindNeighbors(obj, reduction = method_info$new_reduction, dims = 1:30)
  obj <- FindClusters(obj, resolution = 0.5, cluster.name = method_info$cluster_name)
  obj <- RunUMAP(obj, reduction = method_info$new_reduction, dims = 1:30, reduction.name = method_info$umap_name)

  method <- method_name
  umap_by_sample <- DimPlot(obj, reduction = method_info$umap_name, group.by = "sample", label = FALSE, pt.size = 0.6) +
                    ggtitle(paste("UMAP -", toupper(method), "- Integrated Samples"))

  umap_by_cluster <- DimPlot(obj, reduction = method_info$umap_name, group.by = method_info$cluster_name, label = TRUE, pt.size = 0.6) +
                     ggtitle(paste("UMAP -", toupper(method), "- Integrated Clusters"))

  ggsave(file.path(plot_dir, paste0("LGG_", method, "_Integrated_UMAP.pdf")), plot = umap_by_sample, width = 7, height = 6)
  ggsave(file.path(plot_dir, paste0("LGG_", method, "_Integrated_UMAP_Clusters.pdf")), plot = umap_by_cluster, width = 7, height = 6)
  ggsave(file.path(plot_dir, paste0("LGG_", method, "_Combined_UMAP.pdf")), plot = umap_by_sample + umap_by_cluster, width = 14, height = 6)

  saveRDS(obj, file = file.path(rds_dir, paste0("seurat5_integrated_", method_name, ".rds")))
  write.csv(obj@meta.data, file = file.path(rds_dir, paste0("metadata_", method_name, ".csv")))
  save.image(file = paste0("seurat5_integrated_", method_name, ".RData"))
  cat("Finished", method_name, "integration.\n\n")
}

save.image(file = "seurat5_integrated_methods.all.RData")
cat("Integration complete. All outputs saved.\n")
