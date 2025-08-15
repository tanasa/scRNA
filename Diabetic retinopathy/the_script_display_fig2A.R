#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################

#set.seed(423)

#library(Cairo)
#options(bitmapType='cairo')
#library(sf)
#library(rgl)
#library(patchwork)

#library(monocle3)
#library(plotly)

#library(ggplot2)
#library(dplyr)
#library(edgeR)

#library(Seurat)
#library(SeuratData)
#library(SeuratWrappers)

#library(conos)
#library(liger)
#library(harmony)

#library(cowplot)
#library(Matrix)
#library(gridExtra)

#library(stringr)
#library(magrittr)
#library(dplyr)
#library(purrr)

#library(future)
#library(future.apply)

#library(data.table)
#library(vroom)
#library(tidyverse)

#plan("multiprocess", workers = 8)
#options(future.globals.maxSize = 64000 * 1024^2)

##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################

Idents(samples.combined.list.combined) <- "integrated_snn_res.0.4"
Idents(samples.combined.list.combined) <- "seurat_clusters"

##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
################################################################################## SUI would like to plot the data with different RESOLUTIONS

RESOLUTION_0 = 0

RESOLUTION_02 = 0.2

RESOLUTION_04 = 0.4

##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################

cluster0 = subset(x = samples.combined.list.combined, subset = (seurat_clusters == "0") ) 
cluster3 = subset(x = samples.combined.list.combined, subset = (seurat_clusters == "3") ) 
cluster8 = subset(x = samples.combined.list.combined, subset = (seurat_clusters == "8") )
cluster11 = subset(x = samples.combined.list.combined, subset = (seurat_clusters == "11") ) 
cluster12 = subset(x = samples.combined.list.combined, subset = (seurat_clusters == "12") )
cluster15 = subset(x = samples.combined.list.combined, subset = (seurat_clusters == "15") )
cluster17 = subset(x = samples.combined.list.combined, subset = (seurat_clusters == "17") )
cluster28 = subset(x = samples.combined.list.combined, subset = (seurat_clusters == "28") )
cluster30 = subset(x = samples.combined.list.combined, subset = (seurat_clusters == "30") ) 
cluster44 = subset(x = samples.combined.list.combined, subset = (seurat_clusters == "44") )

##################################################################################
################################################################################## CLUSTER_ID = 0
##################################################################################
##################################################################################
####################################################################################################################################################################
####################################################################################################################################################################

CLUSTER_ID = 0

# Visualization using TSNE :
p1 <- DimPlot(cluster0, reduction = "tsne", group.by = "orig.ident")
plot_grid(p1)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.TSNE", "png", sep="."), width = 20, height = 20, units = "cm")

# Visualization using UMAP :
p2 <- DimPlot(cluster0, reduction = "umap", group.by = "orig.ident")
plot_grid(p2)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.UMAP", "png", sep="."), width = 20, height = 20, units = "cm")

##################################################################################
################################################################################## using less DIMENSIONS for speed purposes : 1:10
##################################################################################
################################################################################## RESOLUTION_0
CLUSTER_ID = 0
CLUSTER = cluster0
DefaultAssay(CLUSTER) <- "integrated"
CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_0)
CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

TSNEPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

UMAPPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_02
CLUSTER_ID = 0
CLUSTER = cluster0
DefaultAssay(CLUSTER) <- "integrated"
CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_02)
CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

TSNEPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

UMAPPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_04
CLUSTER_ID = 0
CLUSTER = cluster0
DefaultAssay(CLUSTER) <- "integrated"
CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_04)
CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

TSNEPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

UMAPPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## 
################################################################################## CLUSTER_ID = 3
##################################################################################
##################################################################################
##################################################################################
##################################################################################
####################################################################################################################################################################
####################################################################################################################################################################

CLUSTER_ID = 3

# Visualization using TSNE :
p1 <- DimPlot(cluster3, reduction = "tsne", group.by = "orig.ident")
plot_grid(p1)
ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.TSNE", "png", sep="."), width = 20, height = 20, units = "cm")

# Visualization using UMAP :
p2 <- DimPlot(cluster3, reduction = "umap", group.by = "orig.ident")
plot_grid(p2)
ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.UMAP", "png", sep="."), width = 20, height = 20, units = "cm")

##################################################################################
################################################################################## using less DIMENSIONS for speed purposes : 1:10
##################################################################################
################################################################################## RESOLUTION_0
CLUSTER_ID = 3
CLUSTER = cluster3
DefaultAssay(CLUSTER) <- "integrated"
CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_0)
CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

TSNEPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

UMAPPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_02
CLUSTER_ID = 3
CLUSTER = cluster3
DefaultAssay(CLUSTER) <- "integrated"
CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_02)
CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

TSNEPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

UMAPPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_04
CLUSTER_ID = 3
CLUSTER = cluster3
DefaultAssay(CLUSTER) <- "integrated"
CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_04)
CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

TSNEPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

UMAPPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
##################################################################################
################################################################################## using less DIMENSIONS for speed purposes : 1:10
##################################################################################
##################################################################################
##################################################################################
####################################################################################################################################################################
####################################################################################################################################################################
# CLUSTER_ID = 8

# Visualization using TSNE :
# p1 <- DimPlot(cluster8, reduction = "tsne", group.by = "orig.ident")
# plot_grid(p1)
# ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.TSNE", "png", sep="."), width = 20, height = 20, units = "cm")

# Visualization using UMAP :
# p2 <- DimPlot(cluster8, reduction = "umap", group.by = "orig.ident")
# plot_grid(p2)
# ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.UMAP", "png", sep="."), width = 20, height = 20, units = "cm")

##################################################################################
################################################################################## using less DIMENSIONS for speed purposes : 1:10
##################################################################################
################################################################################## RESOLUTION_0
# CLUSTER_ID = 8
# CLUSTER = cluster8
# DefaultAssay(CLUSTER) <- "integrated"
# CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

# CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_0)
# CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
# CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

# TSNEPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

# UMAPPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_02
# CLUSTER_ID = 8
# CLUSTER = cluster8
# DefaultAssay(CLUSTER) <- "integrated"
# CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

# CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_02)
# CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
# CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

# TSNEPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

# UMAPPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_04
# CLUSTER_ID = 8
# CLUSTER = cluster8
# DefaultAssay(CLUSTER) <- "integrated"
# CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

# CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_04)
# CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
# CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

# TSNEPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

# UMAPPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
##################################################################################
################################################################################## using less DIMENSIONS for speed purposes : 1:10
##################################################################################
##################################################################################
##################################################################################
####################################################################################################################################################################
####################################################################################################################################################################

# CLUSTER_ID = 11

# Visualization using TSNE :
# p1 <- DimPlot(cluster11, reduction = "tsne", group.by = "orig.ident")
# plot_grid(p1)
# ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.TSNE", "png", sep="."), width = 20, height = 20, units = "cm")

# Visualization using UMAP :
# p2 <- DimPlot(cluster11, reduction = "umap", group.by = "orig.ident")
# plot_grid(p2)
# ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.UMAP", "png", sep="."), width = 20, height = 20, units = "cm")

##################################################################################
################################################################################## using less DIMENSIONS for speed purposes : 1:10
##################################################################################
################################################################################## RESOLUTION_0
# CLUSTER_ID = 11
# CLUSTER = cluster11
# DefaultAssay(CLUSTER) <- "integrated"
# CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

# CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_0)
# CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
# CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

# TSNEPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

# UMAPPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_02
# CLUSTER_ID = 11
# CLUSTER = cluster11
# DefaultAssay(CLUSTER) <- "integrated"
# CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

# CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_02)
# CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
# CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

# TSNEPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

# UMAPPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_04
# CLUSTER_ID = 11
# CLUSTER = cluster11
# DefaultAssay(CLUSTER) <- "integrated"
# CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

# CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_04)
# CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
# CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

# TSNEPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

# UMAPPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
##################################################################################
################################################################################## using less DIMENSIONS for speed purposes : 1:10
##################################################################################
##################################################################################
##################################################################################
####################################################################################################################################################################
####################################################################################################################################################################
CLUSTER_ID = 12

# Visualization using TSNE :
p1 <- DimPlot(cluster12, reduction = "tsne", group.by = "orig.ident")
plot_grid(p1)
ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.TSNE", "png", sep="."), width = 20, height = 20, units = "cm")

# Visualization using UMAP :
p2 <- DimPlot(cluster12, reduction = "umap", group.by = "orig.ident")
plot_grid(p2)
ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.UMAP", "png", sep="."), width = 20, height = 20, units = "cm")

##################################################################################
################################################################################## using less DIMENSIONS for speed purposes : 1:10
##################################################################################
################################################################################## RESOLUTION_0
CLUSTER_ID = 12
CLUSTER = cluster12
DefaultAssay(CLUSTER) <- "integrated"
CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_0)
CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

TSNEPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

UMAPPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_02
CLUSTER_ID = 12
CLUSTER = cluster12
DefaultAssay(CLUSTER) <- "integrated"
CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_02)
CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

TSNEPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

UMAPPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_04
CLUSTER_ID = 12
CLUSTER = cluster12
DefaultAssay(CLUSTER) <- "integrated"
CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_04)
CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

TSNEPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

UMAPPlot(CLUSTER, label=TRUE)
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
##################################################################################
################################################################################## using less DIMENSIONS for speed purposes : 1:10
##################################################################################
##################################################################################
##################################################################################
####################################################################################################################################################################
####################################################################################################################################################################
# CLUSTER_ID = 15

# Visualization using TSNE :
# p1 <- DimPlot(cluster15, reduction = "tsne", group.by = "orig.ident")
# plot_grid(p1)
# ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.TSNE", "png", sep="."), width = 20, height = 20, units = "cm")

# Visualization using UMAP :
# p2 <- DimPlot(cluster15, reduction = "umap", group.by = "orig.ident")
# plot_grid(p2)
# ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.UMAP", "png", sep="."), width = 20, height = 20, units = "cm")

##################################################################################
################################################################################## 
##################################################################################
################################################################################## RESOLUTION_0
# CLUSTER_ID = 15
# CLUSTER = cluster15
# DefaultAssay(CLUSTER) <- "integrated"
# CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

# CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_0)
# CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
# CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

# TSNEPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

# UMAPPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_02
# CLUSTER_ID = 15
# CLUSTER = cluster15
# DefaultAssay(CLUSTER) <- "integrated"
# CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

# CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_02)
# CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
# CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

# TSNEPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

# UMAPPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_04
# CLUSTER_ID = 15
# CLUSTER = cluster15
# DefaultAssay(CLUSTER) <- "integrated"
# CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

# CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_04)
# CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
# CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

# TSNEPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

# UMAPPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
##################################################################################
################################################################################## using less DIMENSIONS for speed purposes : 1:10
##################################################################################
##################################################################################
##################################################################################
####################################################################################################################################################################
####################################################################################################################################################################
# CLUSTER_ID = 17

# Visualization using TSNE :
# p1 <- DimPlot(cluster17, reduction = "tsne", group.by = "orig.ident")
# plot_grid(p1)
# ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.TSNE", "png", sep="."), width = 20, height = 20, units = "cm")

# Visualization using UMAP :
# p2 <- DimPlot(cluster17, reduction = "umap", group.by = "orig.ident")
# plot_grid(p2)
# ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.UMAP", "png", sep="."), width = 20, height = 20, units = "cm")

##################################################################################
################################################################################## using less DIMENSIONS for speed purposes : 1:10
##################################################################################
################################################################################## RESOLUTION_0
# CLUSTER_ID = 17
# CLUSTER = cluster17
# DefaultAssay(CLUSTER) <- "integrated"
# CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

# CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_0)
# CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
# CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

# TSNEPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

# UMAPPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_02
#CLUSTER_ID = 17
#CLUSTER = cluster17
#DefaultAssay(CLUSTER) <- "integrated"
#CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

#CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_02)
#CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
#CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

#TSNEPlot(CLUSTER, label=TRUE)
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
#DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

#UMAPPlot(CLUSTER, label=TRUE)
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
#DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_04
#CLUSTER_ID = 17
#CLUSTER = cluster17
#DefaultAssay(CLUSTER) <- "integrated"
#CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

#CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_04)
#CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
#CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

#TSNEPlot(CLUSTER, label=TRUE)
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
#DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

#UMAPPlot(CLUSTER, label=TRUE)
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
#DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
##################################################################################
################################################################################## using less DIMENSIONS for speed purposes : 1:10
##################################################################################
##################################################################################
##################################################################################
####################################################################################################################################################################
####################################################################################################################################################################

# CLUSTER_ID = 28
# Visualization using TSNE :
# p1 <- DimPlot(cluster28, reduction = "tsne", group.by = "orig.ident")
# plot_grid(p1)
# ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.TSNE", "png", sep="."), width = 20, height = 20, units = "cm")

# Visualization using UMAP :
# p2 <- DimPlot(cluster28, reduction = "umap", group.by = "orig.ident")
# plot_grid(p2)
# ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.UMAP", "png", sep="."), width = 20, height = 20, units = "cm")

##################################################################################
################################################################################## using less DIMENSIONS for speed purposes : 1:10
##################################################################################
################################################################################## RESOLUTION_0

# CLUSTER_ID = 28
# CLUSTER = cluster28
# DefaultAssay(CLUSTER) <- "integrated"
# CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

# CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_0)
# CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
# CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

# TSNEPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

# UMAPPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_02

# CLUSTER_ID = 28
# CLUSTER = cluster28
# DefaultAssay(CLUSTER) <- "integrated"
# CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

# CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_02)
# CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
# CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

# TSNEPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

# UMAPPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_04

#CLUSTER_ID = 28
#CLUSTER = cluster28
#DefaultAssay(CLUSTER) <- "integrated"
#CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

#CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_04)
#CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
#CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

#TSNEPlot(CLUSTER, label=TRUE)
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
#DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

#UMAPPlot(CLUSTER, label=TRUE)
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
#DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
##################################################################################
################################################################################## using less DIMENSIONS for speed purposes : 1:10
##################################################################################
##################################################################################
##################################################################################
####################################################################################################################################################################
####################################################################################################################################################################

#CLUSTER_ID = 30

# Visualization using TSNE :
#p1 <- DimPlot(cluster30, reduction = "tsne", group.by = "orig.ident")
#plot_grid(p1)
#ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.TSNE", "png", sep="."), width = 20, height = 20, units = "cm")

# Visualization using UMAP :
#p2 <- DimPlot(cluster30, reduction = "umap", group.by = "orig.ident")
#plot_grid(p2)
#ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.UMAP", "png", sep="."), width = 20, height = 20, units = "cm")

##################################################################################
################################################################################## using less DIMENSIONS for speed purposes : 1:10
##################################################################################
################################################################################## RESOLUTION_0
#CLUSTER_ID = 30
#CLUSTER = cluster30
#DefaultAssay(CLUSTER) <- "integrated"
#CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

#CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_0)
#CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
#CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

#TSNEPlot(CLUSTER, label=TRUE)
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
#DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

#UMAPPlot(CLUSTER, label=TRUE)
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
#DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_02
#CLUSTER_ID = 30
#CLUSTER = cluster30
#DefaultAssay(CLUSTER) <- "integrated"
#CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

#CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_02)
#CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
#CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

#TSNEPlot(CLUSTER, label=TRUE)
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
#DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

#UMAPPlot(CLUSTER, label=TRUE)
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
#DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_04
#CLUSTER_ID = 30
#CLUSTER = cluster30
#DefaultAssay(CLUSTER) <- "integrated"
#CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

#CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_04)
#CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
#CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

#TSNEPlot(CLUSTER, label=TRUE)
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
#DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

#UMAPPlot(CLUSTER, label=TRUE)
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
#DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
#ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
##################################################################################
################################################################################## using less DIMENSIONS for speed purposes : 1:10
##################################################################################
##################################################################################
##################################################################################
####################################################################################################################################################################
####################################################################################################################################################################

# CLUSTER_ID = 44

# Visualization using TSNE :
# p1 <- DimPlot(cluster44, reduction = "tsne", group.by = "orig.ident")
# plot_grid(p1)
# ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.TSNE", "png", sep="."), width = 20, height = 20, units = "cm")

# Visualization using UMAP :
# p2 <- DimPlot(cluster44, reduction = "umap", group.by = "orig.ident")
# plot_grid(p2)
# ggsave(paste("figure2A.sui.Dim.Plot", "CLUSTER", CLUSTER_ID, "assay.INTEGRATED.display.UMAP", "png", sep="."), width = 20, height = 20, units = "cm")

##################################################################################
################################################################################## 
##################################################################################
################################################################################## RESOLUTION_0

# CLUSTER_ID = 44
# CLUSTER = cluster44
# DefaultAssay(CLUSTER) <- "integrated"
# CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

# CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_0)
# CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
# CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

# TSNEPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

# UMAPPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_0, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_02

# CLUSTER_ID = 44
# CLUSTER = cluster44
# DefaultAssay(CLUSTER) <- "integrated"
# CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

# CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_02)
# CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
# CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

# TSNEPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

# UMAPPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_02, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##################################################################################
################################################################################## RESOLUTION_04
# CLUSTER_ID = 44
# CLUSTER = cluster44
# DefaultAssay(CLUSTER) <- "integrated"
# CLUSTER <- FindNeighbors(CLUSTER, dims = 1:10)

# CLUSTER <- FindClusters(CLUSTER, resolution = RESOLUTION_04)
# CLUSTER <- RunTSNE(CLUSTER, dims = 1:10)
# CLUSTER <- RunUMAP(CLUSTER, dims = 1:10)

# TSNEPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "tsne", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "TSNE", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

# UMAPPlot(CLUSTER, label=TRUE)
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "png", sep="."), width = 40, height = 40, units = "cm")
# DimPlot(CLUSTER, reduction = "umap", group.by = "orig.ident")
# ggsave(paste("figure2A.sui", "CLUSTER", CLUSTER_ID, "UMAP", "for.RESOLUTION", RESOLUTION_04, "v2.png", sep="."), width = 40, height = 40, units = "cm")

##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################
##########################################################################################################################################################################################################################################################################