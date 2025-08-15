###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

library(Cairo)
options(bitmapType='cairo')

library(conos)
library(liger)
library(harmony)

library(dplyr)
library(edgeR)

library(Seurat)
library(SeuratData)
library(SeuratWrappers)

library(cowplot)
library(Matrix)

library(gridExtra)
library(ggplot2)

library(future)
library(future.apply)

plan("multiprocess", workers = 8)
options(future.globals.maxSize = 64000 * 1024^2)

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

# Identify cell types that are present in both datasets
# Obtain cell type markers that are conserved in both control and stimulated cells
# Compare the datasets to find cell-type specific responses to stimulation

###########################################################################################################
########################################################################################################### Setup the Seurat objects

# InstallData("ifnb")
# data("ifnb")

# str(ifnb)
# Formal class 'Seurat' [package "Seurat"] with 12 slots
#  ..@ assays      :List of 1
#  .. ..$ RNA:Formal class 'Assay' [package "Seurat"] with 8 slots
#  .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#  .. .. .. .. .. ..@ i       : int [1:9787436] 20 27 37 64 65 83 87 131 139 175 ...
#  .. .. .. .. .. ..@ p       : int [1:14000] 0 877 1590 2440 3549 4183 4740 5720 6301 7181 ...
#  .. .. .. .. .. ..@ Dim     : int [1:2] 14053 13999
#  .. .. .. .. .. ..@ Dimnames:List of 2
#  .. .. .. .. .. .. ..$ : chr [1:14053] "AL627309.1" "RP11-206L10.2" "LINC00115" "NOC2L" ...
#  .. .. .. .. .. .. ..$ : chr [1:13999] "AAACATACATTTCC.1" "AAACATACCAGAAA.1" "AAACATACCTCGCT.1" "AAACATACCTGGTA.1" ...
#  .. .. .. .. .. ..@ x       : num [1:9787436] 1 1 1 1 1 2 1 1 1 1 ...
#  .. .. .. .. .. ..@ factors : list()
#  .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#  .. .. .. .. .. ..@ i       : int [1:9787436] 20 27 37 64 65 83 87 131 139 175 ...
#  .. .. .. .. .. ..@ p       : int [1:14000] 0 877 1590 2440 3549 4183 4740 5720 6301 7181 ...
#  .. .. .. .. .. ..@ Dim     : int [1:2] 14053 13999
#  .. .. .. .. .. ..@ Dimnames:List of 2
#  .. .. .. .. .. .. ..$ : chr [1:14053] "AL627309.1" "RP11-206L10.2" "LINC00115" "NOC2L" ...
#  .. .. .. .. .. .. ..$ : chr [1:13999] "AAACATACATTTCC.1" "AAACATACCAGAAA.1" "AAACATACCTCGCT.1" "AAACATACCTGGTA.1" ...
#  .. .. .. .. .. ..@ x       : num [1:9787436] 1 1 1 1 1 2 1 1 1 1 ...
#  .. .. .. .. .. ..@ factors : list()
#  .. .. .. ..@ scale.data   : num[0 , 0 ] 
#  .. .. .. ..@ key          : chr "rna_"
#  .. .. .. ..@ var.features : logi(0) 
#  .. .. .. ..@ meta.features:'data.frame':	14053 obs. of  0 variables
#  .. .. .. ..@ misc         : symbol NULL
#  .. .. .. ..@ NA           : NULL
#  ..@ meta.data   :'data.frame':	13999 obs. of  5 variables:
#  .. ..$ orig.ident        : chr [1:13999] "IMMUNE_CTRL" "IMMUNE_CTRL" "IMMUNE_CTRL" "IMMUNE_CTRL" ...
#  .. ..$ nCount_RNA        : num [1:13999] 3017 2481 3420 3156 1868 ...
#  .. ..$ nFeature_RNA      : int [1:13999] 877 713 850 1109 634 557 980 581 880 669 ...
#  .. ..$ stim              : chr [1:13999] "CTRL" "CTRL" "CTRL" "CTRL" ...
#  .. ..$ seurat_annotations: Factor w/ 13 levels "CD14 Mono","CD4 Naive T",..: 1 1 1 12 3 1 7 2 6 1 ...
#  ..@ active.assay: chr "RNA"
#  ..@ active.ident: Factor w/ 2 levels "IMMUNE_CTRL",..: 1 1 1 1 1 1 1 1 1 1 ...
#  .. ..- attr(*, "names")= chr [1:13999] "AAACATACATTTCC.1" "AAACATACCAGAAA.1" "AAACATACCTCGCT.1" "AAACATACCTGGTA.1" ...
#  ..@ graphs      : list()
#  ..@ neighbors   : list()
#  ..@ reductions  : list()
#  ..@ project.name: chr "ifnb"
#  ..@ misc        : list()
#  ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
#  .. ..$ : int [1:3] 3 0 0
#  ..@ commands    : list()
#  ..@ tools       : list()

# write.table(as.data.frame(ifnb@meta.data), 
#            file=paste(NAME, "meta.data.START.txt", sep="."),
#            quote=FALSE, sep="\t",   
#            row.names = TRUE, col.names = TRUE)    

# colnames(as.data.frame(ifnb[["RNA"]]@counts))
# rownames(as.data.frame(ifnb[["RNA"]]@counts))

# write.table(colnames(as.data.frame(ifnb[["RNA"]]@counts)), 
#            file=paste(NAME, "meta.data.START.COLNAMES.txt", sep="."),
#            quote=FALSE, sep="\t",   
#            row.names = TRUE, col.names = TRUE)  


# write.table(rownames(as.data.frame(ifnb[["RNA"]]@counts)), 
#            file=paste(NAME, "meta.data.START.ROWNAMES.txt", sep="."),
#            quote=FALSE, sep="\t",   
#            row.names = TRUE, col.names = TRUE)  

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

# ifnb.list <- SplitObject(ifnb, split.by = "stim")

# $CTRL
# An object of class Seurat 
# 14053 features across 6548 samples within 1 assay 
# Active assay: RNA (14053 features)

# $STIM
# An object of class Seurat 
# 14053 features across 7451 samples within 1 assay 
# Active assay: RNA (14053 features)

################################################## the matrices look the same :
# ifnb.list$CTRL[["RNA"]]@counts
# ifnb.list$CTRL[["RNA"]]@data

################################################## the matrices look the same :
# ifnb.list$STIM[["RNA"]]@counts
# ifnb.list$STIM[["RNA"]]@data

########################################################################################################### SPLIT
########################################################################################################### OBJECT
# info regarding SPLIT OBJECT : 
# assign the test object a three level attribute : an example from the MANUAL
########################################################################################################### an 
########################################################################################################### example

# groups <- sample(c("group1", "group2", "group3"), size = 80, replace = TRUE)
# names(groups) <- colnames(samples_small)

# head(samples_small@meta.data)
#                   orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8
# ATGCCAGAACGACT SeuratProject         70           47               0
# CATGGCCTGTGCAT SeuratProject         85           52               0
# GAACCTGATGAACC SeuratProject         87           50               1
# TGACTGGATTCTCA SeuratProject        127           56               0
# AGTCAGACTGCACA SeuratProject        173           53               0
# TCTGATACACGTGT SeuratProject         70           48               0
#               letter.idents groups RNA_snn_res.1
# ATGCCAGAACGACT             A     g2             0
# CATGGCCTGTGCAT             A     g1             0
# GAACCTGATGAACC             B     g2             0
# TGACTGGATTCTCA             A     g2             0
# AGTCAGACTGCACA             A     g2             0
# TCTGATACACGTGT             A     g1             0

# samples_small <- AddMetaData(object = samples_small, metadata = groups, col.name = "group")

# head(samples_small@meta.data)
#                   orig.ident nCount_RNA nFeature_RNA RNA_snn_res.0.8
# ATGCCAGAACGACT SeuratProject         70           47               0
# CATGGCCTGTGCAT SeuratProject         85           52               0
# GAACCTGATGAACC SeuratProject         87           50               1
# TGACTGGATTCTCA SeuratProject        127           56               0
# AGTCAGACTGCACA SeuratProject        173           53               0
# TCTGATACACGTGT SeuratProject         70           48               0
#                letter.idents groups RNA_snn_res.1  group
# ATGCCAGAACGACT             A     g2             0 group2
# CATGGCCTGTGCAT             A     g1             0 group2
# GAACCTGATGAACC             B     g2             0 group3

# Splits object based on a single attribute into a list of subsetted
#     objects, one for each level of the attribute. For example, useful
#     for taking an object that contains cells from many patients, and
#     subdividing it into patient-specific objects.
# SplitObject(object, split.by = "ident")
     
# obj.list <- SplitObject(samples_small, split.by = "group")

# Splits object based on a single attribute into a list of subsetted
#     objects, one for each level of the attribute. For example, useful
#     for taking an object that contains cells from many patients, and
#     subdividing it into patient-specific objects.
# SplitObject(object, split.by = "ident")

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
# from a tutorial on ZINB-WAVE : 

################################################## the matrices look the same :
# ifnb.list$CTRL[["RNA"]]@counts
# ifnb.list$CTRL[["RNA"]]@data

################################################## the matrices look the same :
# ifnb.list$STIM[["RNA"]]@counts
# ifnb.list$STIM[["RNA"]]@data

###########################################################################################################
########################################################################################################### OUTPUT the folder CTRL

# library("Matrix")

# data_mtx_CTRL <- GetAssayData(ifnb.list$CTRL[["RNA"]], slot="counts")

# write.table(rownames(data_mtx_CTRL), file="features.data_mtx_CTRL.tsv", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
# write.table(colnames(data_mtx_CTRL), file="barcodes.data_mtx_CTRL.tsv", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

# Matrix::writeMM(data_mtx_CTRL, file="data_mtx_CTRL.mtx")

########################################################################################################### OUTPUT the folder STIM
###########################################################################################################

# library("Matrix")

# data_mtx_STIM <- GetAssayData(ifnb.list$STIM[["RNA"]], slot="counts")

# write.table(rownames(data_mtx_STIM), file="features.data_mtx_STIM.tsv", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
# write.table(colnames(data_mtx_STIM), file="barcodes.data_mtx_STIM.tsv", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

# Matrix::writeMM(data_mtx_STIM, file="data_mtx_STIM.mtx")

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
############################################################################################### to VERIFY the PIPELINE after we read again :
###############################################################################################

# DIR_CTRL="/media/bogdan/e7a05079-e4ce-4cfa-b0cd-bf2d7a54b46d/PIPELINE_SEURAT3/data_CTRL_again4/"

# DIR_STIM="/media/bogdan/e7a05079-e4ce-4cfa-b0cd-bf2d7a54b46d/PIPELINE_SEURAT3/data_STIM_again4/"

###############################################################################################
###############################################################################################

# ctrl.data <- Read10X(data.dir = DIR_CTRL)
# stim.data <- Read10X(data.dir = DIR_STIM)

# colnames(ctrl.data) <- paste(CTRL, colnames(ctrl.data), sep="-")
# colnames(stim.data) <- paste(STIM, colnames(stim.data), sep="-")

# read10X(mtx="matrix.mtx", genes="features.tsv", barcodes="barcodes.tsv")
# Error in read10X(mtx = "matrix.mtx", genes = "features.tsv", barcodes = "barcodes.tsv") 
# it does not work ...
# read10X(mtx="data_mtx_CTRL.mtx", genes="features.data_mtx_CTRL.tsv", barcodes="barcodes.data_mtx_CTRL.tsv")
 
# str(data_mtx_CTRL)
#Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#  ..@ i       : int [1:4626015] 20 27 37 64 65 83 87 131 139 175 ...
#  ..@ p       : int [1:6549] 0 877 1590 2440 3549 4183 4740 5720 6301 7181 ...
#  ..@ Dim     : int [1:2] 14053 6548
#  ..@ Dimnames:List of 2
#  .. ..$ : chr [1:14053] "AL627309.1" "RP11-206L10.2" "LINC00115" "NOC2L" ...
#  .. ..$ : chr [1:6548] "AAACATACATTTCC.1" "AAACATACCAGAAA.1" "AAACATACCTCGCT.1" "AAACATACCTGGTA.1" ...
#  ..@ x       : num [1:4626015] 1 1 1 1 1 2 1 1 1 1 ...
#  ..@ factors : list()

# str(data_mtx_STIM)
# Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#  ..@ i       : int [1:5161421] 7 10 20 25 27 79 87 113 176 194 ...
#  ..@ p       : int [1:7452] 0 588 1380 1964 2696 3242 4558 5140 5702 6210 ...
#  ..@ Dim     : int [1:2] 14053 7451
#  ..@ Dimnames:List of 2
#  .. ..$ : chr [1:14053] "AL627309.1" "RP11-206L10.2" "LINC00115" "NOC2L" ...
#  .. ..$ : chr [1:7451] "AAACATACCAAGCT.1" "AAACATACCCCTAC.1" "AAACATACCCGTAA.1" "AAACATACCCTCGT.1" ...
#  ..@ x       : num [1:5161421] 15 1 1 1 1 1 1 1 1 1 ...
#  ..@ factors : list()

####################################################################################################################### 
####################################################################################################################### 
####################################################################################################################### 
####################################################################################################################### 
####################################################################################################################### 
####################################################################################################################### 
####################################################################################################################### 
####################################################################################################################### NAMES

CTRL1="WT.1M"
STIM1="DIS.1M"

CTRL2="WT.2M"
STIM2="DIS.2M"

CTRL3="WT.3M"
STIM3="DIS.3M"

CTRL4="WT.4M"
STIM4="DIS.4M"

CTRL5="WT.5M"
STIM5="DIS.5M"

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

DIR_CTRL1 = "/labs/scRNAseq/cell_aggregate_samples_1M_WT/"
DIR_STIM1 = "/labs/scRNAseq/cell_aggregate_samples_1M_DIS/"

DIR_CTRL2 = "/labs/scRNAseq/cell_aggregate_samples_2M_WT"
DIR_STIM2 = "/labs/scRNAseq/cell_aggregate_samples_2M_DIS"

DIR_CTRL3 = "/labs/scRNAseq/cell_aggregate_samples_3M_WT"
DIR_STIM3 = "/labs/scRNAseq/cell_aggregate_samples_3M_DIS"

DIR_CTRL4 = "/labs/scRNAseq/cell_aggregate_samples_4M_WT"
DIR_STIM4 = "/labs/scRNAseq/cell_aggregate_samples_4M_DIS"

DIR_CTRL5 = "/labs/scRNAseq/cell_aggregate_samples_5M_WT"
DIR_STIM5 = "/labs/scRNAseq/cell_aggregate_samples_5M_DIS"

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

NAME <- paste("integrated", sep=".")

LIST.KNOWN.MARKERS = read.delim("a.LIST.MARKERS.RETINA.mouse.txt", header=T, sep="\t", stringsAsFactors=F) ### to read the LIST at the start..

### in the folders above :

###############################################################################################
############################################################################################### THRESHOLDS : use 600 !

RESOLUTION = 0.6 

THRESHOLD_min_counts <- 500 
THRESHOLD_max_counts <- 20000 

THRESHOLD_min_genes <- 600 
THRESHOLD_max_genes <- 2500 

THRESHOLD_mito <- 20

###############################################################################################
################################################################################################
############################################################################################### MERGE VIGNETTE

# https://satijalab.org/seurat/v3.1/merge_vignette.html

# ctrl.data <- read.table("control_expression_matrix.txt.gz", sep = "\t")
# stim.data <- read.table("stimulated_expression_matrix.txt.gz", sep = "\t")

###############################################################################################
###############################################################################################

ctrl1.data <- Read10X(data.dir = DIR_CTRL1)
stim1.data <- Read10X(data.dir = DIR_STIM1)

colnames(ctrl1.data) <- paste(CTRL1, colnames(ctrl1.data), sep="-")
colnames(stim1.data) <- paste(STIM1, colnames(stim1.data), sep="-")

###############################################################################################
###############################################################################################

ctrl2.data <- Read10X(data.dir = DIR_CTRL2)
stim2.data <- Read10X(data.dir = DIR_STIM2)

colnames(ctrl2.data) <- paste(CTRL2, colnames(ctrl2.data), sep="-")
colnames(stim2.data) <- paste(STIM2, colnames(stim2.data), sep="-")

###############################################################################################
###############################################################################################

ctrl3.data <- Read10X(data.dir = DIR_CTRL3)
stim3.data <- Read10X(data.dir = DIR_STIM3)

colnames(ctrl3.data) <- paste(CTRL3, colnames(ctrl3.data), sep="-")
colnames(stim3.data) <- paste(STIM3, colnames(stim3.data), sep="-")

###############################################################################################
###############################################################################################

ctrl4.data <- Read10X(data.dir = DIR_CTRL4)
stim4.data <- Read10X(data.dir = DIR_STIM4)

colnames(ctrl4.data) <- paste(CTRL4, colnames(ctrl4.data), sep="-")
colnames(stim4.data) <- paste(STIM4, colnames(stim4.data), sep="-")

###############################################################################################
###############################################################################################

ctrl5.data <- Read10X(data.dir = DIR_CTRL5)
stim5.data <- Read10X(data.dir = DIR_STIM5)

colnames(ctrl5.data) <- paste(CTRL5, colnames(ctrl5.data), sep="-")
colnames(stim5.data) <- paste(STIM5, colnames(stim5.data), sep="-")

####################################################################################################################### 
####################################################################################################################### 
####################################################################################################################### 
####################################################################################################################### 

##############################################################################################################################################################################################################################################  
##############################################################################################################################################################################################################################################  
####################################################################################################################### 
####################################################################################################################### making SEURAT OBJECTS : 1
####################################################################################################################### 

ctrl1 <- CreateSeuratObject(ctrl1.data, project = CTRL1)

head(colnames(ctrl1[["RNA"]]@counts))
head(ctrl1@meta.data)
# GetAssayData(ctrl1)[1:10, 1:10]
# as.data.frame(ctrl1[["RNA"]]@counts[c("Rbpms"),])

write.table(as.data.frame(ctrl1[["RNA"]]@counts[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.COUNTS", "Rbpms", "in.assay.RNA.sample", CTRL1, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

write.table(as.data.frame(ctrl1[["RNA"]]@data[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.DATA", "Rbpms", "in.assay.RNA.sample", CTRL1, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

##################################################################################################################
##################################################################################################################

stim1 <- CreateSeuratObject(stim1.data, project = STIM1)

head(colnames(stim1[["RNA"]]@counts))
head(stim1@meta.data)
# GetAssayData(stim1)[1:10, 1:10]
# as.data.frame(stim1[["RNA"]]@counts[c("Rbpms"),])

write.table(as.data.frame(stim1[["RNA"]]@counts[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.COUNTS", "Rbpms", "in.assay.RNA.sample", STIM1, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

write.table(as.data.frame(stim1[["RNA"]]@data[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.DATA", "Rbpms", "in.assay.RNA.sample", STIM1, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

##############################################################################################################################################################################################################################################  
##############################################################################################################################################################################################################################################  
####################################################################################################################### 
####################################################################################################################### making SEURAT OBJECTS : 2
#######################################################################################################################  

ctrl2 <- CreateSeuratObject(ctrl2.data, project = CTRL2)

head(colnames(ctrl2[["RNA"]]@counts))
head(ctrl2@meta.data)
# GetAssayData(ctrl2)[1:10, 1:10]
# as.data.frame(ctrl2[["RNA"]]@counts[c("Rbpms"),])

write.table(as.data.frame(ctrl2[["RNA"]]@counts[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.COUNTS", "Rbpms", "in.assay.RNA.sample", CTRL2, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

write.table(as.data.frame(ctrl2[["RNA"]]@data[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.DATA", "Rbpms", "in.assay.RNA.sample", CTRL2, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

##################################################################################################################
##################################################################################################################

stim2 <- CreateSeuratObject(stim2.data, project = STIM2)

head(colnames(stim2[["RNA"]]@counts))
head(stim2@meta.data)
# GetAssayData(stim2)[1:10, 1:10]
# as.data.frame(stim2[["RNA"]]@counts[c("Rbpms"),])

write.table(as.data.frame(stim2[["RNA"]]@counts[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.COUNTS", "Rbpms", "in.assay.RNA.sample", STIM2, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

write.table(as.data.frame(stim2[["RNA"]]@data[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.DATA", "Rbpms", "in.assay.RNA.sample", STIM2, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

##############################################################################################################################################################################################################################################  
##############################################################################################################################################################################################################################################  
####################################################################################################################### 
####################################################################################################################### making SEURAT OBJECTS : 3
#######################################################################################################################  

ctrl3 <- CreateSeuratObject(ctrl3.data, project = CTRL3)

head(colnames(ctrl3[["RNA"]]@counts))
head(ctrl3@meta.data)
# GetAssayData(ctrl3)[1:10, 1:10]
# as.data.frame(ctrl3[["RNA"]]@counts[c("Rbpms"),])

write.table(as.data.frame(ctrl3[["RNA"]]@counts[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.COUNTS", "Rbpms", "in.assay.RNA.sample", CTRL3, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

write.table(as.data.frame(ctrl3[["RNA"]]@data[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.DATA", "Rbpms", "in.assay.RNA.sample", CTRL3, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

##################################################################################################################
##################################################################################################################

stim3 <- CreateSeuratObject(stim3.data, project = STIM3)

head(colnames(stim3[["RNA"]]@counts))
head(stim3@meta.data)
# GetAssayData(stim3)[1:10, 1:10]
# as.data.frame(stim3[["RNA"]]@counts[c("Rbpms"),])

write.table(as.data.frame(stim3[["RNA"]]@counts[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.COUNTS", "Rbpms", "in.assay.RNA.sample", STIM3, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

write.table(as.data.frame(stim3[["RNA"]]@data[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.DATA", "Rbpms", "in.assay.RNA.sample", STIM3, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

##############################################################################################################################################################################################################################################  
##############################################################################################################################################################################################################################################  
####################################################################################################################### 
####################################################################################################################### making SEURAT OBJECTS : 4
#######################################################################################################################  

ctrl4 <- CreateSeuratObject(ctrl4.data, project = CTRL4)

head(colnames(ctrl4[["RNA"]]@counts))
head(ctrl4@meta.data)
# GetAssayData(ctrl4)[1:10, 1:10]
# as.data.frame(ctrl4[["RNA"]]@counts[c("Rbpms"),])

write.table(as.data.frame(ctrl4[["RNA"]]@counts[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.COUNTS", "Rbpms", "in.assay.RNA.sample", CTRL4, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

write.table(as.data.frame(ctrl4[["RNA"]]@data[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.DATA", "Rbpms", "in.assay.RNA.sample", CTRL4, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

##################################################################################################################
##################################################################################################################

stim4 <- CreateSeuratObject(stim4.data, project = STIM4)

head(colnames(stim4[["RNA"]]@counts))
head(stim4@meta.data)
# GetAssayData(stim4)[1:10, 1:10]
# as.data.frame(stim4[["RNA"]]@counts[c("Rbpms"),])

write.table(as.data.frame(stim4[["RNA"]]@counts[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.COUNTS", "Rbpms", "in.assay.RNA.sample", STIM4, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

write.table(as.data.frame(stim4[["RNA"]]@data[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.DATA", "Rbpms", "in.assay.RNA.sample", STIM4, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

##############################################################################################################################################################################################################################################  
##############################################################################################################################################################################################################################################  
####################################################################################################################### 
####################################################################################################################### making SEURAT OBJECTS : 5
#######################################################################################################################  

ctrl5 <- CreateSeuratObject(ctrl5.data, project = CTRL5)

head(colnames(ctrl5[["RNA"]]@counts))
head(ctrl5@meta.data)
# GetAssayData(ctrl5)[1:10, 1:10]
# as.data.frame(ctrl5[["RNA"]]@counts[c("Rbpms"),])

write.table(as.data.frame(ctrl5[["RNA"]]@counts[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.COUNTS", "Rbpms", "in.assay.RNA.sample", CTRL5, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

write.table(as.data.frame(ctrl5[["RNA"]]@data[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.DATA", "Rbpms", "in.assay.RNA.sample", CTRL5, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

##################################################################################################################
##################################################################################################################

stim5 <- CreateSeuratObject(stim5.data, project = STIM5)

head(colnames(stim5[["RNA"]]@counts))
head(stim5@meta.data)
# GetAssayData(stim5)[1:10, 1:10]
# as.data.frame(stim5[["RNA"]]@counts[c("Rbpms"),])

write.table(as.data.frame(stim5[["RNA"]]@counts[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.COUNTS", "Rbpms", "in.assay.RNA.sample", STIM5, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

write.table(as.data.frame(stim5[["RNA"]]@data[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.DATA", "Rbpms", "in.assay.RNA.sample", STIM5, "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

####################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################
################################################################################################################## it adds ID to specific CELLS
##################################################################################################################
# samples.combined <- merge(x=ctrl, y = stim, add.cell.ids = c(CTRL, STIM), project = "TUTORIAL")
####################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################
##################################################################################################################
##################################################################################################################

samples.combined1 <- merge(x=ctrl1, y = stim1)
samples.combined1

head(colnames(samples.combined1))
tail(colnames(samples.combined1))

head(samples.combined1@meta.data)
tail(samples.combined1@meta.data)
 
table(samples.combined1$orig.ident)

head(colnames(samples.combined1[["RNA"]]@counts))
tail(colnames(samples.combined1[["RNA"]]@counts))

head(colnames(samples.combined1[["RNA"]]@data))
tail(colnames(samples.combined1[["RNA"]]@data))

##############################################################################################################################################################################################################################################  
##############################################################################################################################################################################################################################################  
################################################################################################################## we can use all the samples or, 
################################################################################################################## if we wanna reduce the number of samples, we do get the following 

# samples.combined2 <- merge(x=samples.combined1, y = c(ctrl2, stim2, ctrl3, stim3, ctrl4, stim4, ctrl5, stim5))

samples.combined2 <- merge(x=samples.combined1, y = c(ctrl3, stim3, ctrl5, stim5))
samples.combined2

# samples.combined2
# An object of class Seurat 
# 25399 features across 119074 samples within 1 assay 
# Active assay: RNA (25399 features, 0 variable features)

head(colnames(samples.combined2))
tail(colnames(samples.combined2))

head(samples.combined2@meta.data)
tail(samples.combined2@meta.data)
 
table(samples.combined2$orig.ident)

head(colnames(samples.combined2[["RNA"]]@counts))
tail(colnames(samples.combined2[["RNA"]]@counts))

head(colnames(samples.combined2[["RNA"]]@data))
tail(colnames(samples.combined2[["RNA"]]@data))

#  STZ.1M STZ.2M STZ.3M STZ.4M STZ.5M  WT.1M  WT.2M  WT.3M  WT.4M  WT.5M 
#  9500  13749  15088  14882  13797   9206   8946  13996  10701   9209 

# to print RBPMS for verifications :

write.table(as.data.frame(samples.combined2[["RNA"]]@counts[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.COUNTS", "Rbpms", "in.assay.RNA.sample", "COMBINED", "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

write.table(as.data.frame(samples.combined2[["RNA"]]@data[c("Rbpms"),]), 
                        file=paste(NAME, "figure2.DATA", "Rbpms", "in.assay.RNA.sample", "COMBINED", "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

##############################################################################################################################################################################################################################################  
##############################################################################################################################################################################################################################################  
##############################################################################################################################################################################################################################################  
##############################################################################################################################################################################################################################################  

dim(samples.combined2[["RNA"]]@counts)
dim(samples.combined2[["RNA"]]@data) 

# dim(samples.combined2[["RNA"]]@data)
# [1]  25399 119074
# dim(samples.combined2[["RNA"]]@counts)
# [1]  25399 119074

dim(samples.combined2@meta.data) 
# 70796     3
# 119074      3

write.table(as.data.frame(samples.combined2@meta.data), 
                        file=paste(NAME, "figure2.DATA", "the.META.DATA", "in.assay.RNA.sample", "COMBINED", "txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

##############################################################################################################################################################################################################################################  
##############################################################################################################################################################################################################################################  
##############################################################################################################################################################################################################################################  
##############################################################################################################################################################################################################################################  
 
### to print the INTEGRATED MATRIX

# write.table(as.data.frame(samples.combined2[["RNA"]]@counts), 
#                        file=paste(NAME, "figure2.DATA", "the.COUNTS.DATA", "in.assay.RNA.sample", "COMBINED", "txt", sep="."), 
#                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

# write.table(as.data.frame(samples.combined2[["RNA"]]@data), 
#                        file=paste(NAME, "figure2.DATA", "the.DATA.DATA", "in.assay.RNA.sample", "COMBINED", "txt", sep="."), 
#                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

##############################################################################################################################################################################################################################################  
##############################################################################################################################################################################################################################################  
##############################################################################################################################################################################################################################################  
##############################################################################################################################################################################################################################################  
 
### to save the memory 

rm(ctrl1)
rm(stim1)

rm(ctrl2)
rm(stim2)

rm(ctrl3)
rm(stim3)

rm(ctrl4)
rm(stim4)

rm(ctrl5)
rm(stim5)

### to continue the analysis using the variable : samples.combined

samples.combined = samples.combined2

rm(samples.combined2)
rm(samples.combined1)

##############################################################################################################################################################################################################################################  
##############################################################################################################################################################################################################################################  
############################################################################################################################################################################################################################################## 
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
################################################################################################################## compute the QC on the MITO GENES 

# rownames(samples.combined[["RNA"]]@data)
# samples.combined[["percent.mt"]] <- PercentageFeatureSet(samples.combined, pattern = "^MT-")        #### 've change to MT

samples.combined[["percent.mt"]] <- PercentageFeatureSet(samples.combined, pattern = "^mt-|^MT-|^Mt") #### 've change to MT

# from the other script :
# mito.genes <- grep(pattern = "^mt-|^MT-|^Mt", x = rownames(x=pbmc[["RNA"]]@data), value = TRUE)
# mito.genes

head(samples.combined@meta.data)
tail(samples.combined@meta.data)

VlnPlot(samples.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(paste(NAME, "figure1.samples.combined", "VIOLIN.PLOT.genes.counts.mito.plots", "png", sep="."), width=38, height=18, units = "cm") 

# group.by: Group (color) cells in different ways (for example, orig.ident) ###### by DEFAULT
# split.by: A variable to split the violin plots by, see ‘FetchData’ for more details

# rownames(samples.combined[["RNA"]]@counts)
# rownames(samples.combined[["RNA"]]@data)

###########################################################################################################
########################################################################################################### DISPLAY

plot1 <- FeatureScatter(samples.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(samples.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

CombinePlots(plots = list(plot1, plot2))
ggsave(paste(NAME, "figure1.samples.combined", "FEATURE.SCATTER.genes.counts.mito.plots", "png", sep="."), width=38, height=18, units = "cm") 

###########################################################################################################
###########################################################################################################
########################################################################################################### shall we FILTER the INTEGRATED DATA at this STEP : 
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
# THRESHOLD_genes <- 500 
# THRESHOLD_mito <- 20

# samples.combined <- subset(samples.combined, percent.mt <= THRESHOLD_mito)
# samples.combined

# THRESHOLD_min_genes <- 500 
# THRESHOLD_max_genes <- 20000 
# samples.combined <- subset(samples.combined, nCount_RNA >= THRESHOLD_min_counts & nCount_RNA <= THRESHOLD_max_counts)

# samples.combined <- subset(samples.combined, nFeature_RNA > THRESHOLD_min_genes & nFeature_RNA < THRESHOLD_max_genes)
# samples.combined

samples.combined <- subset(samples.combined, nFeature_RNA > THRESHOLD_min_genes & percent.mt < THRESHOLD_mito)
samples.combined

# THRESHOLD_min_counts <- 500 
# THRESHOLD_max_counts <- 20000 

###########################################################################################################
###########################################################################################################
###########################################################################################################
########################################################################################################### to verify if the pipeline gives the same results
###########################################################################################################
########################################################################################################### using the same parameters as at the begining

# from tutorial : samples.combined.list : ifnb = samples.combined

samples.combined.list <- SplitObject(samples.combined, split.by = "orig.ident")

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

for (i in 1:length(samples.combined.list)) 
{
    samples.combined.list[[i]] <- NormalizeData(samples.combined.list[[i]]) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
}

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

# pbmcsca.panel = samples.combined.list

samples.combined.list.con <- Conos$new(samples.combined.list)

samples.combined.list.con$buildGraph(k = 15, 
                                     k.self = 5, 
                                     space = "PCA", 
                                     ncomps = 30, 
                                     n.odgenes = 2000, 
                                     matching.method = "mNN", 
                                     metric = "angular", 
                                     score.component.variance = TRUE, verbose = TRUE)

samples.combined.list.con$findCommunities()

samples.combined.list.con$embedGraph()

samples.combined.list.CONOS <- as.Seurat(samples.combined.list.con)

DimPlot(samples.combined.list.CONOS, reduction = "largeVis", group.by = c("orig.ident"), ncol = 1, label=TRUE)
ggsave(paste(NAME, "figure1.Dim.Plot", "CONOS", "orig.ident", "png", sep="."), width = 40, height = 40, units = "cm")

DimPlot(samples.combined.list.CONOS, reduction = "largeVis", group.by = c("leiden"), ncol = 1, label=TRUE)
ggsave(paste(NAME, "figure1.Dim.Plot", "CONOS", "leiden", "png", sep="."), width = 40, height = 40, units = "cm")

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

str(samples.combined.list.CONOS)

#Formal class 'Seurat' [package "Seurat"] with 12 slots
#  ..@ assays      :List of 1
#  .. ..$ RNA:Formal class 'Assay' [package "Seurat"] with 8 slots
#  .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#  .. .. .. .. .. ..@ i       : int [1:6614279] 195 252 320 321 347 391 428 493 517 615 ...
#  .. .. .. .. .. ..@ p       : int [1:7292] 0 536 1081 1694 2386 2953 3454 4096 5531 6339 ...
#  .. .. .. .. .. ..@ Dim     : int [1:2] 25399 7291
#  .. .. .. .. .. ..@ Dimnames:List of 2
#  .. .. .. .. .. .. ..$ : chr [1:25399] "AABR07000089.1" "Vom2r6" "Vom2r5" "Raet1e" ...
#  .. .. .. .. .. .. ..$ : chr [1:7291] "WT.1M-AAACCTGAGATCCCAT" "WT.1M-AAACCTGAGGGCTCTC" "WT.1M-AAACCTGCAGTCACTA" "WT.1M-AAACCTGCATATACCG" ...
#  .. .. .. .. .. ..@ x       : num [1:6614279] 1 1 1 1 1 2 1 1 1 1 ...
#  .. .. .. .. .. ..@ factors : list()
#  .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#  .. .. .. .. .. ..@ i       : int [1:6614279] 195 252 320 321 347 391 428 493 517 615 ...
#  .. .. .. .. .. ..@ p       : int [1:7292] 0 536 1081 1694 2386 2953 3454 4096 5531 6339 ...
#  .. .. .. .. .. ..@ Dim     : int [1:2] 25399 7291
#  .. .. .. .. .. ..@ Dimnames:List of 2
#  .. .. .. .. .. .. ..$ : chr [1:25399] "AABR07000089.1" "Vom2r6" "Vom2r5" "Raet1e" ...
#  .. .. .. .. .. .. ..$ : chr [1:7291] "WT.1M-AAACCTGAGATCCCAT" "WT.1M-AAACCTGAGGGCTCTC" "WT.1M-AAACCTGCAGTCACTA" "WT.1M-AAACCTGCATATACCG" ...
#  .. .. .. .. .. ..@ x       : num [1:6614279] 2.62 2.62 2.62 2.62 2.62 ...
#  .. .. .. .. .. ..@ factors : list()
#  .. .. .. ..@ scale.data   : num[0 , 0 ] 
#  .. .. .. ..@ key          : chr "rna_"
#  .. .. .. ..@ assay.orig   : NULL
#  .. .. .. ..@ var.features : logi(0) 
#  .. .. .. ..@ meta.features:'data.frame':	25399 obs. of  0 variables
#  .. .. .. ..@ misc         : NULL
#  ..@ meta.data   :'data.frame':	7291 obs. of  5 variables:
#  .. ..$ orig.ident  : chr [1:7291] "WT.1M" "WT.1M" "WT.1M" "WT.1M" ...
#  .. ..$ nCount_RNA  : num [1:7291] 786 806 882 1016 833 ...
#  .. ..$ nFeature_RNA: int [1:7291] 536 545 613 692 567 501 642 1435 808 2127 ...
#  .. ..$ percent.mt  : num [1:7291] 8.14 5.71 6.01 5.91 5.52 ...
#  .. ..$ leiden      : Factor w/ 19 levels "1","2","3","4",..: 5 2 5 2 3 1 1 4 5 10 ...
#  ..@ active.assay: chr "RNA"
#  ..@ active.ident: Factor w/ 19 levels "1","2","3","4",..: 5 2 5 2 3 1 1 4 5 10 ...
#  .. ..- attr(*, "names")= chr [1:7291] "WT.1M-AAACCTGAGATCCCAT" "WT.1M-AAACCTGAGGGCTCTC" "WT.1M-AAACCTGCAGTCACTA" "WT.1M-AAACCTGCATATACCG" ...
#  ..@ graphs      :List of 1
#  .. ..$ RNA_mnn:Formal class 'Graph' [package "Seurat"] with 7 slots
#  .. .. .. ..@ assay.used: NULL
#  .. .. .. ..@ i         : int [1:97260] 101 574 578 970 2454 2793 7150 7231 1097 1263 ...
#  .. .. .. ..@ p         : int [1:7292] 0 8 17 31 47 55 76 83 106 114 ...
#  .. .. .. ..@ Dim       : int [1:2] 7291 7291
#  .. .. .. ..@ Dimnames  :List of 2
#  .. .. .. .. ..$ : chr [1:7291] "WT.1M-AAACCTGAGATCCCAT" "WT.1M-AAACCTGAGGGCTCTC" "WT.1M-AAACCTGCAGTCACTA" "WT.1M-AAACCTGCATATACCG" ...
#  .. .. .. .. ..$ : chr [1:7291] "WT.1M-AAACCTGAGATCCCAT" "WT.1M-AAACCTGAGGGCTCTC" "WT.1M-AAACCTGCAGTCACTA" "WT.1M-AAACCTGCATATACCG" ...
#  .. .. .. ..@ x         : num [1:97260] 1 1 1 1 1 1 1 1 1 1 ...
#  .. .. .. ..@ factors   : list()
#  ..@ neighbors   : list()
#  ..@ reductions  :List of 1
#  .. ..$ largeVis:Formal class 'DimReduc' [package "Seurat"] with 9 slots
#  .. .. .. ..@ cell.embeddings           : num [1:7291, 1:2] 2.03 1.65 2.67 2.62 1.28 ...
#  .. .. .. .. ..- attr(*, "dimnames")=List of 2
#  .. .. .. .. .. ..$ : chr [1:7291] "WT.1M-AAACCTGAGATCCCAT" "WT.1M-AAACCTGAGGGCTCTC" "WT.1M-AAACCTGCAGTCACTA" "WT.1M-AAACCTGCATATACCG" ...
#  .. .. .. .. .. ..$ : chr [1:2] "LARGEVIS_1" "LARGEVIS_2"
#  .. .. .. ..@ feature.loadings          : num[0 , 0 ] 
#  .. .. .. ..@ feature.loadings.projected: num[0 , 0 ] 
#  .. .. .. ..@ assay.used                : chr "RNA"
#  .. .. .. ..@ global                    : logi FALSE
#  .. .. .. ..@ stdev                     : num(0) 
#  .. .. .. ..@ key                       : chr "LARGEVIS_"
#  .. .. .. ..@ jackstraw                 :Formal class 'JackStrawData' [package "Seurat"] with 4 slots
#  .. .. .. .. .. ..@ empirical.p.values     : num[0 , 0 ] 
#  .. .. .. .. .. ..@ fake.reduction.scores  : num[0 , 0 ] 
#  .. .. .. .. .. ..@ empirical.p.values.full: num[0 , 0 ] 
#  .. .. .. .. .. ..@ overall.p.values       : num[0 , 0 ] 
#  .. .. .. ..@ misc                      : list()
#  ..@ project.name: chr "SeuratProject"
#  ..@ misc        :List of 2
#  .. ..$ conos.pairs        :List of 1
#  .. .. ..$ PCA:List of 1
#  .. .. .. ..$ WT.1M.vs.STZ.1M:List of 2
#  .. .. .. .. ..$ CPC: num [1:2000, 1:60] 5.46e-04 2.28e-05 5.55e-03 3.64e-02 5.68e-02 ...
#  .. .. .. .. .. ..- attr(*, "dimnames")=List of 2
#  .. .. .. .. .. .. ..$ : chr [1:2000] "AABR07001734.1" "AABR07003304.2" "AABR07006724.1" "AABR07007032.1" ...
#  .. .. .. .. .. .. ..$ : NULL
#  .. .. .. .. ..$ nv :List of 2
#  .. .. .. .. .. ..$ WT.1M : num [1:30] 0.172 0.0515 0.0363 0.0182 0.0103 ...
#  .. .. .. .. .. ..$ STZ.1M: num [1:30] 0.1924 0.0344 0.0238 0.0189 0.0133 ...
#  .. ..$ conos.leiden.result:List of 6
#  .. .. ..$ membership: Factor w/ 19 levels "0","1","2","3",..: 2 2 2 2 2 2 3 2 3 2 ...
#  .. .. .. ..- attr(*, "names")= chr [1:7291] "WT.1M-AATCGGTAGCTAGTCT" "STZ.1M-AAACCTGAGCCACTAT" "WT.1M-ACGAGGAGTCCGAACC" "WT.1M-ACTTGTTGTAAGTTCC" ...
#  .. .. ..$ dendrogram: NULL
#  .. .. ..$ algorithm : chr "leiden"
#  .. .. ..$ resolution: num 1
#  .. .. ..$ n.iter    : num 2
#  .. .. ..$ names     : chr [1:7291] "WT.1M-AATCGGTAGCTAGTCT" "STZ.1M-AAACCTGAGCCACTAT" "WT.1M-ACGAGGAGTCCGAACC" "WT.1M-ACTTGTTGTAAGTTCC" ...
#  ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
#  .. ..$ : int [1:3] 3 1 3
#  ..@ commands    : list()
#  ..@ tools       : list()

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
# samples.combined.list.CONOS@meta.data

write.table(as.data.frame(samples.combined.list.CONOS@meta.data), 
                        file=paste(NAME, "figure1.META.DATA", "CONOS.txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

samples.combined.list.combined = samples.combined.list.CONOS 

########################################################################################################### 
########################################################################################################### 
########################################################################################################### 
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

# samples.combined.list <- lapply(X = samples.combined.list, FUN = function(x) {
#     x <- NormalizeData(x)
#     x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })

########################################################################################################### Perform integration
###########################################################################################################
###########################################################################################################
###########################################################################################################

# We identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input, 
# and use these anchors to integrate the two datasets together with IntegrateData.

# samples.combined.list.anchors <- FindIntegrationAnchors(object.list = samples.combined.list, dims = 1:20)

# samples.combined.list.combined <- IntegrateData(anchorset = samples.combined.list.anchors, dims = 1:20)

# str(samples.combined.list.combined)
# Formal class 'Seurat' [package "Seurat"] with 12 slots
#  ..@ assays      :List of 2
#  .. ..$ RNA       :Formal class 'Assay' [package "Seurat"] with 8 slots
#  .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#  .. .. .. .. .. ..@ i       : int [1:8212162] 195 252 320 321 347 391 428 493 517 615 ...
#  .. .. .. .. .. ..@ p       : int [1:10355] 0 536 958 1503 1909 2522 3214 3781 4282 4924 ...
#  .. .. .. .. .. ..@ Dim     : int [1:2] 25399 10354
#  .. .. .. .. .. ..@ Dimnames:List of 2
#  .. .. .. .. .. .. ..$ : chr [1:25399] "AABR07000089.1" "Vom2r6" "Vom2r5" "Raet1e" ...
#  .. .. .. .. .. .. ..$ : chr [1:10354] "WT.1M-AAACCTGAGATCCCAT" "WT.1M-AAACCTGAGGAGTAGA" "WT.1M-AAACCTGAGGGCTCTC" "WT.1M-AAACCTGCAAAGGAAG" ...
#  .. .. .. .. .. ..@ x       : num [1:8212162] 1 1 1 1 1 2 1 1 1 1 ...
#  .. .. .. .. .. ..@ factors : list()
#  .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#  .. .. .. .. .. ..@ i       : int [1:8212162] 195 252 320 321 347 391 428 493 517 615 ...
#  .. .. .. .. .. ..@ p       : int [1:10355] 0 536 958 1503 1909 2522 3214 3781 4282 4924 ...
#  .. .. .. .. .. ..@ Dim     : int [1:2] 25399 10354
#  .. .. .. .. .. ..@ Dimnames:List of 2
#  .. .. .. .. .. .. ..$ : chr [1:25399] "AABR07000089.1" "Vom2r6" "Vom2r5" "Raet1e" ...
#  .. .. .. .. .. .. ..$ : chr [1:10354] "WT.1M-AAACCTGAGATCCCAT" "WT.1M-AAACCTGAGGAGTAGA" "WT.1M-AAACCTGAGGGCTCTC" "WT.1M-AAACCTGCAAAGGAAG" ...
#  .. .. .. .. .. ..@ x       : num [1:8212162] 2.62 2.62 2.62 2.62 2.62 ...
#  .. .. .. .. .. ..@ factors : list()
#  .. .. .. ..@ scale.data   : num[0 , 0 ] 
#  .. .. .. ..@ key          : chr "rna_"
#  .. .. .. ..@ assay.orig   : NULL
#  .. .. .. ..@ var.features : logi(0) 
#  .. .. .. ..@ meta.features:'data.frame':	25399 obs. of  0 variables
#  .. .. .. ..@ misc         : NULL
#  .. ..$ integrated:Formal class 'Assay' [package "Seurat"] with 8 slots
#  .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#  .. .. .. .. .. ..@ i       : int(0) 
#  .. .. .. .. .. ..@ p       : int 0
#  .. .. .. .. .. ..@ Dim     : int [1:2] 0 0
#  .. .. .. .. .. ..@ Dimnames:List of 2
#  .. .. .. .. .. .. ..$ : NULL
#  .. .. .. .. .. .. ..$ : NULL
#  .. .. .. .. .. ..@ x       : num(0) 
#  .. .. .. .. .. ..@ factors : list()
# .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#  .. .. .. .. .. ..@ i       : int [1:3639944] 1 5 132 249 283 332 391 419 451 477 ...
#  .. .. .. .. .. ..@ p       : int [1:10355] 0 34 82 127 165 215 268 307 345 384 ...
#  .. .. .. .. .. ..@ Dim     : int [1:2] 2000 10354
#  .. .. .. .. .. ..@ Dimnames:List of 2
#  .. .. .. .. .. .. ..$ : chr [1:2000] "Arr3" "Cst3" "Pde6h" "Cartpt" ...
#  .. .. .. .. .. .. ..$ : chr [1:10354] "WT.1M-AAACCTGAGATCCCAT" "WT.1M-AAACCTGAGGAGTAGA" "WT.1M-AAACCTGAGGGCTCTC" "WT.1M-AAACCTGCAAAGGAAG" ...
#  .. .. .. .. .. ..@ x       : num [1:3639944] 3.95 3.28 2.62 2.62 2.62 ...
#  .. .. .. .. .. ..@ factors : list()
#  .. .. .. ..@ scale.data   : logi [1, 1] NA
#  .. .. .. ..@ key          : chr "integrated_"
#  .. .. .. ..@ assay.orig   : NULL
#  .. .. .. ..@ var.features : chr [1:2000] "Arr3" "Cst3" "Pde6h" "Cartpt" ...
#  .. .. .. ..@ meta.features:'data.frame':	2000 obs. of  0 variables
#  .. .. .. ..@ misc         : NULL
#  ..@ meta.data   :'data.frame':	10354 obs. of  3 variables:
#  .. ..$ orig.ident  : chr [1:10354] "WT.1M" "WT.1M" "WT.1M" "WT.1M" ...
#  .. ..$ nCount_RNA  : num [1:10354] 786 660 806 598 882 ...
#  .. ..$ nFeature_RNA: int [1:10354] 536 422 545 406 613 692 567 501 642 1435 ...
#  ..@ active.assay: chr "integrated"
# etc etc etc

#####################################################################################################################################
#####################################################################################################################################

# write.table(as.data.frame(samples.combined.list.combined@meta.data), 
#                         file=paste(NAME, "META.DATA", "in.assay.RNA.at.the.start.SCRIPT.txt", sep="."), 
#                         sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

# head(samples.combined.list.combined@meta.data)
#                       orig.ident nCount_RNA nFeature_RNA percent.mt
# WT.1M-AAACCTGAGATCCCAT      WT.1M        786          536   8.142494
# WT.1M-AAACCTGAGGGCTCTC      WT.1M        806          545   5.707196

#####################################################################################################################################
#####################################################################################################################################
###########################################################################################################
########################################################################################################### Perform an INTEGRATED ANALYSIS 
###########################################################################################################
########################################################################################################### it changes the DEFAULT ASSAY to "INTEGRATED"

# Now we can run a single integrated analysis on all cells !
# however, we will change the group.by = "stim" with group.by = "orig.ident" 

# DefaultAssay(samples.combined.list.combined) <- "integrated"

# Run the standard workflow for visualization and clustering  

# samples.combined.list.combined <- ScaleData(samples.combined.list.combined, verbose = FALSE)

# samples.combined.list.combined <- RunPCA(samples.combined.list.combined, npcs = 30, verbose = FALSE)

# TSNE and UMAP 
# samples.combined.list.combined <- RunTSNE(samples.combined.list.combined, reduction = "pca", dims = 1:20)
# samples.combined.list.combined <- RunUMAP(samples.combined.list.combined, reduction = "pca", dims = 1:20)

# Clustering
# samples.combined.list.combined <- FindNeighbors(samples.combined.list.combined, reduction = "pca", dims = 1:20)
# samples.combined.list.combined <- FindClusters(samples.combined.list.combined, resolution = RESOLUTION)

##################################################################################
# Visualization using TSNE :
# p1 <- DimPlot(samples.combined.list.combined, reduction = "tsne", group.by = "orig.ident", label = TRUE)
# p2 <- DimPlot(samples.combined.list.combined, reduction = "tsne", label = TRUE)
# plot_grid(p1, p2)
# ggsave(paste(NAME, "figure2.Dim.Plot", "assay.INTEGRATED.display.TSNE", "png", sep="."), width = 40, height = 20, units = "cm")

##################################################################################
# Visualization using UMAP :
# p1 <- DimPlot(samples.combined.list.combined, reduction = "umap", group.by = "orig.ident", label = TRUE)
# p2 <- DimPlot(samples.combined.list.combined, reduction = "umap", label = TRUE)
# plot_grid(p1, p2)
# ggsave(paste(NAME, "figure2.Dim.Plot", "assay.INTEGRATED.display.UMAP", "png", sep="."), width = 40, height = 20, units = "cm")

###################################################################################
# To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.

# DimPlot(samples.combined.list.combined, reduction = "tsne", split.by = "orig.ident", label = TRUE)
# ggsave(paste(NAME, "figure2.Dim.Plot", "assay.INTEGRATED.display2.TSNE", "png", sep="."), width = 40, height = 20, units = "cm")

# DimPlot(samples.combined.list.combined, reduction = "umap", split.by = "orig.ident", label = TRUE)
# ggsave(paste(NAME, "figure2.Dim.Plot", "assay.INTEGRATED.display2.UMAP", "png", sep="."), width = 40, height = 20, units = "cm")

###########################################################################################################
########################################################################################################### to print META DATA in "INTEGRATED"

# write.table(as.data.frame(samples.combined.list.combined@meta.data), 
#                        file=paste(NAME, "figure2.META.DATA", "in.assay.COMBINED.txt", sep="."), 
#                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

###########################################################################################################
########################################################################################################### Identify CONSERVED CELL TYPE MARKERS
###########################################################################################################
###########################################################################################################

# To identify canonical cell type marker genes that are conserved across conditions, we provide the FindConservedMarkers function. 
# This function performs differential gene expression testing for each dataset/group and combines the p-values using meta-analysis methods from the MetaDE R package. 
# For example, we can calculated the genes that are conserved markers irrespective of stimulation condition in cluster 7 (NK cells).

###########################################################################################################
########################################################################################################### it changes the DEFAULT ASSAY to "RNA"
###########################################################################################################
###########################################################################################################

# about the CLUSTERS : in LEIDEN SLOT 
# orig.ident	nCount_RNA	nFeature_RNA	percent.mt	leiden

DefaultAssay(samples.combined.list.combined) <- "RNA"

###########################################################################################################
###########################################################################################################
########################################################################################################### to print META DATA in "RNA" 
########################################################################################################### the files are the same as 've printed before ...

write.table(as.data.frame(samples.combined.list.combined@meta.data), 
                        file=paste(NAME, "figure2.META.DATA", "in.assay.RNA.txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
########################################################################################################### finding DIFF MARKERS
########################################################################################################### between a CLUSTER 
########################################################################################################### and all other CLUSTERS
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

samples.combined.list.combined.markers.clusters <- FindAllMarkers(object = samples.combined.list.combined, only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)

###################################################################################################################################
################################################################################################################################### THE MARKERS
###################################################################################################################################
### printing ALL MARKERS :

samples.combined.list.combined.markers.clusters.all <- samples.combined.list.combined.markers.clusters %>% group_by(cluster)                      

samples.combined.list.combined.markers.clusters.all <- as.data.frame(samples.combined.list.combined.markers.clusters.all)

write.table(samples.combined.list.combined.markers.clusters.all, 
            file=paste(NAME, "figure5.samples.combined.here.CELL.CLUSTER.SPECIFIC.MARKERS.between.CLUSTERS.indep.of.STIM.all.txt",  sep="."), 
            sep="\t", quote=F, row.names=T, col.names=T)

####################################################################################################################################
#################################################################################################################################### printing the TOP 12
### printing the TOP 12 : if tare at least 12 MARKERS in each CLUSTER, if not we print as many as we have ..,
####################################################################################################################################

# samples.combined.list.combined.markers.clusters.top <- samples.combined.list.combined.markers.clusters %>% group_by(cluster) %>% top_n(12, avg_logFC)                     
# samples.combined.list.combined.markers.clusters.top <- as.data.frame(samples.combined.list.combined.markers.clusters.top)

####################################################################################################################################
#################################################################################################################################### printing the TOP 3

samples.combined.list.combined.markers.clusters.top <- samples.combined.list.combined.markers.clusters %>% group_by(cluster) %>% top_n(3, avg_logFC)
samples.combined.list.combined.markers.clusters.top <- as.data.frame(samples.combined.list.combined.markers.clusters.top)

write.table(samples.combined.list.combined.markers.clusters.top, 
            file=paste(NAME, "figure5.samples.combined.here.CELL.CLUSTER.SPECIFIC.MARKERS.between.CLUSTERS.indep.of.STIM.top3.avg_logFC.txt",  sep="."), 
            sep="\t", quote=F, row.names=T, col.names=T)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# samples.combined.list.combined.markers.clusters.top.list.markers <- split(samples.combined.list.combined.markers.clusters.top, 
#                                                                   samples.combined.list.combined.markers.clusters.top$cluster)

samples.combined.list.combined.markers.clusters.all.list.markers <- split(samples.combined.list.combined.markers.clusters.all, 
                                                                          samples.combined.list.combined.markers.clusters.all$cluster)

####################################################################################################################################
####################################################################################################################################

# start CLUSTERING from 1 in CONOS 

####################################################################################################################################
####################################################################################################################################

for (i in 1:length(samples.combined.list.combined.markers.clusters.all.list.markers)){
    
     # print(i) 
     # LENGTH_LIST = dim(samples.combined.list.combined.markers.clusters.all.list.markers[[i]])[1]  
     # printing the MINIMUM between LENGTH_LIST, and 12 (in such a way that the lists shorter than 3 are still printed)  

if (dim(samples.combined.list.combined.markers.clusters.all.list.markers[[i]])[1] > 0){
    
    for (j in 1:min(dim(samples.combined.list.combined.markers.clusters.all.list.markers[[i]])[1], 3) )    ##### printing top 3 MARKERS in each CLUSTER or the MIN 
    {
         gene <- samples.combined.list.combined.markers.clusters.all.list.markers[[i]]$gene[j] 
         print(gene)  

        # VIOLIN PLOT 
        # png(paste("figure", "cluster", i, "top.gene", j, "plot.VIOLIN.for", gene, "png", sep="."),  width=20, height=20, units="cm", res=400)
          VlnPlot(object = samples.combined.list.combined, features = gene)
          ggsave(paste(NAME, "figure5.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.VIOLIN.for", gene, "png", sep="."), 
                 width=40, height=20, units="cm")
        # dev.off()

        # FEATURE PLOT
        # png(paste("figure", "cluster", i, "top.gene", j, "plot.FEATURE.for", gene, "png", sep="."),  width=10, height=10, units="cm", res=400)
        #  FeaturePlot(object = samples.combined.list.combined, features = gene , cols = c("lightgrey", "blue"), reduction = "tsne")        
        #  ggsave(paste(NAME, "figure5.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.FEATURE.TSNE.for", gene, "png", sep="."))
        # dev.off()

        # FEATURE PLOT
        # png(paste("figure", "cluster", i, "top.gene", j, "plot.FEATURE.for", gene, "png", sep="."),  width=10, height=10, units="cm", res=400)
        #  FeaturePlot(object = samples.combined.list.combined, features = gene , cols = c("lightgrey", "blue"), reduction = "umap")  
        #  ggsave(paste(NAME, "figure5.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.FEATURE.UMAP.for", gene, "png", sep="."))      
        # dev.off()

        # FEATURE PLOT
        # png(paste("figure", "cluster", i, "top.gene", j, "plot.FEATURE.for", gene, "png", sep="."),  width=10, height=10, units="cm", res=400)
          FeaturePlot(object = samples.combined.list.combined, features = gene , cols = c("lightgrey", "blue"), reduction = "largeVis")  
          ggsave(paste(NAME, "figure5.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.FEATURE.largeVis.for", gene, "png", sep="."), 
                 width=20, height=20, units="cm")      
        # dev.off()

        # DOT PLOT
        # png(paste("figure", "cluster", i, "top.gene", j, "plot.DOT.for", gene, "png", sep="."),  width=30, height=30, units="cm", res=400)
          DotPlot(object = samples.combined.list.combined, features = gene)
          ggsave(paste(NAME, "figure5.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.DOT", gene, "png", sep="."), 
                 width=20, height=20, units="cm") 
        # dev.off()

        # RIDGE PLOT 
        # png(paste("figure", "cluster", i, "top.gene", j, "plot.RIDGE.for", gene, "png", sep="."),  width=30, height=30, units="cm", res=400)
          RidgePlot(object = samples.combined.list.combined, features = gene)
          ggsave(paste(NAME, "figure5.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.RIDGE", gene, "png", sep="."), 
                 width=20, height=20, units="cm") 
        # dev.off()
     }
  }
}

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# THE NUMBER OF CELLS in each CLUSTER :
table(Idents(samples.combined.list.combined))
write.table(table(Idents(samples.combined.list.combined)), 
                             file=paste(NAME, "figure5.samples.combined", "NUMBER.CELLS.per.CLUSTER.considering.combined.samples.txt", sep="."), 
                             sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)

# THE NUMBER OF CLUSTERS is :
dim(table(Idents(samples.combined.list.combined)))

write.table(dim(table(Idents(samples.combined.list.combined))), 
                             file=paste(NAME, "figure5.samples.combined", "NUMBER.CLUSTERS.considering.combined.samples.txt", sep="."), 
                             sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
#################################################################################################################################### 
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
#################################################################################################################################### 
#################################################################################################################################### 

# DefaultAssay(samples.combined.list.combined) <- "RNA"
# start CLUSTERING from 1 in CONOS 

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

TABLE_CELLS_CLUSTERS_SAMPLES = table(samples.combined.list.combined@meta.data$orig.ident, 
                                     samples.combined.list.combined@meta.data$leiden)
         
TABLE_CELLS_CLUSTERS_SAMPLES.df = as.data.frame(TABLE_CELLS_CLUSTERS_SAMPLES)

# head(TABLE_CELLS_CLUSTERS_SAMPLES.df)
# THE NUMBER OF CELLS in each CLUSTER per SAMPLE :

write.table(TABLE_CELLS_CLUSTERS_SAMPLES, 
            file=paste(NAME, "figure5.samples.combined", "NUMBER.CELLS.per.CLUSTER.considering.combined.samples.PER.SAMPLE.v1.txt", sep="."), 
                       sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)

write.table(TABLE_CELLS_CLUSTERS_SAMPLES.df, 
            file=paste(NAME, "figure5.samples.combined", "NUMBER.CELLS.per.CLUSTER.considering.combined.samples.PER.SAMPLE.v2.txt", sep="."), 
                       sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)

# table(samples.combined.list.combined@meta.data$orig.ident)

####################################################################################################################################
####################################################################################################################################
#################################################################################################################################### the numbers of cells per CLUSTER
####################################################################################################################################

TABLE.CTRL1 = TABLE_CELLS_CLUSTERS_SAMPLES.df[TABLE_CELLS_CLUSTERS_SAMPLES.df$Var1 == CTRL1, ]
TABLE.STIM1 = TABLE_CELLS_CLUSTERS_SAMPLES.df[TABLE_CELLS_CLUSTERS_SAMPLES.df$Var1 == STIM1, ]

TABLE.CTRL2 = TABLE_CELLS_CLUSTERS_SAMPLES.df[TABLE_CELLS_CLUSTERS_SAMPLES.df$Var1 == CTRL2, ]
TABLE.STIM2 = TABLE_CELLS_CLUSTERS_SAMPLES.df[TABLE_CELLS_CLUSTERS_SAMPLES.df$Var1 == STIM2, ]

TABLE.CTRL3 = TABLE_CELLS_CLUSTERS_SAMPLES.df[TABLE_CELLS_CLUSTERS_SAMPLES.df$Var1 == CTRL3, ]
TABLE.STIM3 = TABLE_CELLS_CLUSTERS_SAMPLES.df[TABLE_CELLS_CLUSTERS_SAMPLES.df$Var1 == STIM3, ]

TABLE.CTRL4 = TABLE_CELLS_CLUSTERS_SAMPLES.df[TABLE_CELLS_CLUSTERS_SAMPLES.df$Var1 == CTRL4, ]
TABLE.STIM4 = TABLE_CELLS_CLUSTERS_SAMPLES.df[TABLE_CELLS_CLUSTERS_SAMPLES.df$Var1 == STIM4, ]

TABLE.CTRL5 = TABLE_CELLS_CLUSTERS_SAMPLES.df[TABLE_CELLS_CLUSTERS_SAMPLES.df$Var1 == CTRL5, ]
TABLE.STIM5 = TABLE_CELLS_CLUSTERS_SAMPLES.df[TABLE_CELLS_CLUSTERS_SAMPLES.df$Var1 == STIM5, ]

####################################################################################################################################
#################################################################################################################################### computing the length of the LISTS

dim(TABLE.CTRL1)[1]
dim(TABLE.STIM1)[1]

dim(TABLE.CTRL2)[1]
dim(TABLE.STIM2)[1]

dim(TABLE.CTRL3)[1]
dim(TABLE.STIM3)[1]

dim(TABLE.CTRL4)[1]
dim(TABLE.STIM4)[1]

dim(TABLE.CTRL5)[1]
dim(TABLE.STIM5)[1]

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
#################################################################################################################################### to compute CONSERVED MARKERS 
####################################################################################################################################

# start CLUSTERING from 1 in CONOS 
# so we have changed some of the indexes 
# dim_array_clusters =  dim(table(Idents(samples.combined.list.combined))) - 1 

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

LIST.CLUSTERS.and.CONSERVED.MARKERS <- list()
dim_array_clusters =  dim(table(Idents(samples.combined.list.combined))) 

for (i in 1:dim_array_clusters) {
if (  (TABLE.CTRL1[TABLE.CTRL1$Var2 == i,]$Freq >=3) && (TABLE.STIM1[TABLE.STIM1$Var2 == i,]$Freq >=3) && 
      (TABLE.CTRL3[TABLE.CTRL3$Var2 == i,]$Freq >=3) && (TABLE.STIM3[TABLE.STIM3$Var2 == i,]$Freq >=3) && 
      (TABLE.CTRL5[TABLE.CTRL5$Var2 == i,]$Freq >=3) && (TABLE.STIM5[TABLE.STIM5$Var2 == i,]$Freq >=3) ){
  
   LIST.CLUSTERS.and.CONSERVED.MARKERS[[i]] <- FindConservedMarkers(samples.combined.list.combined, 
                                                                    ident.1 = i, 
                                                                    grouping.var = "orig.ident", print.bar = FALSE, only.pos = FALSE)

    x <- as.data.frame(as.matrix(LIST.CLUSTERS.and.CONSERVED.MARKERS[[i]]))  
    x$gene <- row.names(x)

    write.table(x, file=paste(NAME, "figure6.samples.combined.CONSERVED.MARKERS.CLUSTER", i, "LIST.txt", sep="."), 
                   sep="\t", quote=F, row.names=T, col.names=T)

    ### printing the top CONSERVED MARKERS in a CLUSTER

    for (j in 1:min(dim(x)[1], 3))         #################################### printing top 3 MARKERS in each CLUSTER or the MIN 
    {
        gene <- x$gene[j] 
        print(gene)   

        ### VLN PLOT
        ### VlnPlot(object = samples.combined.list.combined, features = x$gene[1:6])
        VlnPlot(object = samples.combined.list.combined, features = gene, split.by = "orig.ident",  split.plot = TRUE)  
        ggsave(paste(NAME, "figure6.samples.combined.CONSERVED.MARKERS.CLUSTER", i, "top.gene", j, gene, "VIOLIN.PLOT.png", sep="."), width=100, height=20, units = "cm")

        ### FEATURE PLOT
        ### FeaturePlot(object = samples.combined.list.combined, features = x$gene[1:6] , cols = c("lightgrey", "blue"), reduction = "tsne")   
        # FeaturePlot(object = samples.combined.list.combined, features = gene, cols = c("lightgrey", "blue"), reduction = "tsne", split.by = "orig.ident")        
        # ggsave(paste(NAME, "figure6.samples.combined.CONSERVED.MARKERS.CLUSTER", i, "top.gene", j, gene, "FEATURE.TSNE.png", sep="."), width=60, height=40, units = "cm")

        ### FEATURE PLOT 
        ### FeaturePlot(object = samples.combined.list.combined, features = x$gene[1:6] , cols = c("lightgrey", "blue"), reduction = "umap")    
        # FeaturePlot(object = samples.combined.list.combined, features = gene, cols = c("lightgrey", "blue"), reduction = "umap", split.by = "orig.ident")        
        # ggsave(paste(NAME, "figure6.samples.combined.CONSERVED.MARKERS.CLUSTER", i, "top.gene", j, gene, "FEATURE.UMAP.png", sep="."), width=60, height=40, units = "cm")

        ### FEATURE PLOT 
        ### FeaturePlot(object = samples.combined.list.combined, features = x$gene[1:6] , cols = c("lightgrey", "blue"), reduction = "umap")    
        FeaturePlot(object = samples.combined.list.combined, features = gene, cols = c("lightgrey", "blue"), reduction = "largeVis", split.by = "orig.ident")        
        ggsave(paste(NAME, "figure6.samples.combined.CONSERVED.MARKERS.CLUSTER", i, "top.gene", j, gene, "FEATURE.largeVis.png", sep="."), width=80, height=30, units = "cm")

        ### DOT PLOT
        ### DotPlot(object = samples.combined.list.combined, features = x$gene[1:6])
        # DotPlot(object = samples.combined.list.combined, features = gene,  split.by = "orig.ident")
        # ggsave(paste(NAME, "figure6.samples.combined.CONSERVED.MARKERS.CLUSTER", i, "top.gene", j, gene, "DOT.PLOT.png", sep="."), width=20, height=24, units = "cm") 

        ### RIDGE PLOT 
        ### RidgePlot(object = samples.combined.list.combined, features = x$gene[1:6])
        RidgePlot(object = samples.combined.list.combined, features = gene, group.by = "orig.ident")
        ggsave(paste(NAME, "figure6.samples.combined.CONSERVED.MARKERS.CLUSTER", i, "top.gene", j, gene, "RIDGE.PLOT.png", sep="."), width=40, height=40, units = "cm") 
     }
  }
}

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
############################################################################################################### from SEURAT TUTORIAL

# samples.combined.list.combined <- RenameIdents(samples.combined.list.combined, "0" = "CD14 Mono", "1" = "CD4 Naive T", "2" = "CD4 Memory T", 
#                                                 "3" = "CD16 Mono", "4" = "B", "5" = "CD8 T", "6" = "T activated", 
#                                                 "7" = "NK", "8" = "DC", "9" = "B Activated", "10" = "Mk", 
#                                                 "11" = "pDC", "12" = "Eryth", "13" = "Mono/Mk Doublets")

# DimPlot(samples.combined.list.combined, label = TRUE)
# ggsave(paste(NAME, "Dim.Plot", "after.Rename.Idents", "png", sep="."), width = 40, height = 20, units = "cm")

# The DotPlot function with the split.by parameter can be useful for viewing conserved cell type markers across conditions,
# showing both the expression level and the percentage of cells in a cluster expressing any given gene. 

# Idents(samples.combined.list.combined) <- factor(Idents(samples.combined.list.combined), levels = c("Mono/Mk Doublets", "pDC", 
#                                                          "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated", 
#                                                          "CD4 Naive T", "CD4 Memory T"))

# markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5")

# DotPlot(samples.combined.list.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, split.by = "orig.ident") + RotatedAxis()
# ggsave(paste(NAME, "Dot.Plot", "CD3D", "CREM", "CD8A", "FCGR3A", "png", sep="."), width = 40, height = 20, units = "cm")

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
######################################################################################################## writing again the META DATA

# write.table(as.data.frame(samples.combined.list.combined@meta.data), 
#                        file=paste(NAME, "META.DATA", "in.assay.RNA.after.relabeling.CLUSTERS.txt", sep="."), 
#                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

# but why the RENAMED IDENTITIES were not added to the CLUSTERS and to the META DATA ?

# Formal class 'Seurat' [package "Seurat"] with 12 slots
#  ..@ assays      :List of 2
#  .. ..$ RNA       :Formal class 'Assay' [package "Seurat"] with 8 slots
#  .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#  .. .. .. .. .. ..@ i       : int [1:9787436] 20 27 37 64 65 83 87 131 139 175 ...
#  .. .. .. .. .. ..@ p       : int [1:14000] 0 877 1590 2440 3549 4183 4740 5720 6301 7181 ...
#  .. .. .. .. .. ..@ Dim     : int [1:2] 14053 13999
#  .. .. .. .. .. ..@ Dimnames:List of 2
#  .. .. .. .. .. .. ..$ : chr [1:14053] "AL627309.1" "RP11-206L10.2" "LINC00115" "NOC2L" ...
#  .. .. .. .. .. .. ..$ : chr [1:13999] "AAACATACATTTCC.1" "AAACATACCAGAAA.1" "AAACATACCTCGCT.1" "AAACATACCTGGTA.1" ...
#  .. .. .. .. .. ..@ x       : num [1:9787436] 1 1 1 1 1 2 1 1 1 1 ...
#  .. .. .. .. .. ..@ factors : list()
#  .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#  .. .. .. .. .. ..@ i       : int [1:9787436] 20 27 37 64 65 83 87 131 139 175 ...
#  .. .. .. .. .. ..@ p       : int [1:14000] 0 877 1590 2440 3549 4183 4740 5720 6301 7181 ...
#  .. .. .. .. .. ..@ Dim     : int [1:2] 14053 13999
#  .. .. .. .. .. ..@ Dimnames:List of 2
#  .. .. .. .. .. .. ..$ : chr [1:14053] "AL627309.1" "RP11-206L10.2" "LINC00115" "NOC2L" ...
#  .. .. .. .. .. .. ..$ : chr [1:13999] "AAACATACATTTCC.1" "AAACATACCAGAAA.1" "AAACATACCTCGCT.1" "AAACATACCTGGTA.1" ...
#  .. .. .. .. .. ..@ x       : num [1:9787436] 1.46 1.46 1.46 1.46 1.46 ...
#  .. .. .. .. .. ..@ factors : list()
#  .. .. .. ..@ scale.data   : num[0 , 0 ] 
#  .. .. .. ..@ key          : chr "rna_"
#  .. .. .. ..@ assay.orig   : NULL
#  .. .. .. ..@ var.features : logi(0) 
#  .. .. .. ..@ meta.features:'data.frame':	14053 obs. of  0 variables
#  .. .. .. ..@ misc         : NULL
#  .. ..$ integrated:Formal class 'Assay' [package "Seurat"] with 8 slots
#  .. .. .. ..@ counts       :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#  .. .. .. .. .. ..@ i       : int(0) 
#  .. .. .. .. .. ..@ p       : int 0
#  .. .. .. .. .. ..@ Dim     : int [1:2] 0 0
#  .. .. .. .. .. ..@ Dimnames:List of 2
#  .. .. .. .. .. .. ..$ : NULL
#  .. .. .. .. .. .. ..$ : NULL
#  .. .. .. .. .. ..@ x       : num(0) 
#  .. .. .. .. .. ..@ factors : list()
#  .. .. .. ..@ data         :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#  .. .. .. .. .. ..@ i       : int [1:8216142] 0 1 2 3 4 5 6 7 8 9 ...
#  .. .. .. .. .. ..@ p       : int [1:14000] 0 1244 2359 3544 4848 5914 6928 8051 9085 10245 ...
#  .. .. .. .. .. ..@ Dim     : int [1:2] 2000 13999
#  .. .. .. .. .. ..@ Dimnames:List of 2
#  .. .. .. .. .. .. ..$ : chr [1:2000] "HBB" "HBA2" "HBA1" "CCL4" ...
#  .. .. .. .. .. .. ..$ : chr [1:13999] "AAACATACATTTCC.1" "AAACATACCAGAAA.1" "AAACATACCTCGCT.1" "AAACATACCTGGTA.1" ...
#  .. .. .. .. .. ..@ x       : num [1:8216142] -0.0592 0.0176 0.0127 0.3794 0.7167 ...
#  .. .. .. .. .. ..@ factors : list()
#  .. .. .. ..@ scale.data   : num [1:2000, 1:13999] -0.3907 -0.1129 -0.089 -0.1944 -0.0208 ...
#  .. .. .. .. ..- attr(*, "dimnames")=List of 2
#  .. .. .. .. .. ..$ : chr [1:2000] "HBB" "HBA2" "HBA1" "CCL4" ...
#  .. .. .. .. .. ..$ : chr [1:13999] "AAACATACATTTCC.1" "AAACATACCAGAAA.1" "AAACATACCTCGCT.1" "AAACATACCTGGTA.1" ...
#  .. .. .. ..@ key          : chr "integrated_"
#  .. .. .. ..@ assay.orig   : NULL
#  .. .. .. ..@ var.features : chr [1:2000] "HBB" "HBA2" "HBA1" "CCL4" ...
#  .. .. .. ..@ meta.features:'data.frame':	2000 obs. of  0 variables
#  .. .. .. ..@ misc         : NULL
#  ..@ meta.data   :'data.frame':	13999 obs. of  7 variables:
#  .. ..$ orig.ident            : chr [1:13999] "IMMUNE_CTRL" "IMMUNE_CTRL" "IMMUNE_CTRL" "IMMUNE_CTRL" ...
#  .. ..$ nCount_RNA            : num [1:13999] 3017 2481 3420 3156 1868 ...
#  .. ..$ nFeature_RNA          : int [1:13999] 877 713 850 1109 634 557 980 581 880 669 ...
#  .. ..$ stim                  : chr [1:13999] "CTRL" "CTRL" "CTRL" "CTRL" ...
#  .. ..$ seurat_annotations    : chr [1:13999] "CD14 Mono" "CD14 Mono" "CD14 Mono" "pDC" ...
#  .. ..$ integrated_snn_res.0.5: Factor w/ 13 levels "0","1","2","3",..: 1 1 1 12 3 1 8 2 6 1 ...
#  .. ..$ seurat_clusters       : Factor w/ 13 levels "0","1","2","3",..: 1 1 1 12 3 1 8 2 6 1 ...
#  ..@ active.assay: chr "RNA"
#  ..@ active.ident: Factor w/ 13 levels "pDC","Eryth",..: 5 5 5 1 13 5 10 12 9 5 ...
#  .. ..- attr(*, "names")= chr [1:13999] "AAACATACATTTCC.1" "AAACATACCAGAAA.1" "AAACATACCTCGCT.1" "AAACATACCTGGTA.1" ...
#  ..@ graphs      :List of 2
#  .. ..$ integrated_nn :Formal class 'Graph' [package "Seurat"] with 7 slots
#  .. .. .. ..@ assay.used: chr "integrated"
#  .. .. .. ..@ i         : int [1:279980] 0 68 87 134 184 327 786 887 1177 1402 ...
#  .. .. .. ..@ p         : int [1:14000] 0 53 74 87 101 110 132 149 153 157 ...
#  .. .. .. ..@ Dim       : int [1:2] 13999 13999
#  .. .. .. ..@ Dimnames  :List of 2
#  .. .. .. .. ..$ : chr [1:13999] "AAACATACATTTCC.1" "AAACATACCAGAAA.1" "AAACATACCTCGCT.1" "AAACATACCTGGTA.1" ...
#  .. .. .. .. ..$ : chr [1:13999] "AAACATACATTTCC.1" "AAACATACCAGAAA.1" "AAACATACCTCGCT.1" "AAACATACCTGGTA.1" ...
#  .. .. .. ..@ x         : num [1:279980] 1 1 1 1 1 1 1 1 1 1 ...
#  .. .. .. ..@ factors   : list()
#  .. ..$ integrated_snn:Formal class 'Graph' [package "Seurat"] with 7 slots
#  .. .. .. ..@ assay.used: chr "integrated"
#  .. .. .. ..@ i         : int [1:1081997] 0 68 87 134 182 184 327 379 589 621 ...
#  .. .. .. ..@ p         : int [1:14000] 0 112 168 239 303 404 508 565 656 706 ...
#  .. .. .. ..@ Dim       : int [1:2] 13999 13999
#  .. .. .. ..@ Dimnames  :List of 2
#  .. .. .. .. ..$ : chr [1:13999] "AAACATACATTTCC.1" "AAACATACCAGAAA.1" "AAACATACCTCGCT.1" "AAACATACCTGGTA.1" ...
#  .. .. .. .. ..$ : chr [1:13999] "AAACATACATTTCC.1" "AAACATACCAGAAA.1" "AAACATACCTCGCT.1" "AAACATACCTGGTA.1" ...
#  .. .. .. ..@ x         : num [1:1081997] 1 0.111 0.176 0.212 0.111 ...
#  .. .. .. ..@ factors   : list()
#  ..@ neighbors   : list()
#  ..@ reductions  :List of 2
#  .. ..$ pca :Formal class 'DimReduc' [package "Seurat"] with 9 slots
#  .. .. .. ..@ cell.embeddings           : num [1:13999, 1:30] 11.736 13.666 11.625 -0.908 -8.224 ...
#  .. .. .. .. ..- attr(*, "dimnames")=List of 2
#  .. .. .. .. .. ..$ : chr [1:13999] "AAACATACATTTCC.1" "AAACATACCAGAAA.1" "AAACATACCTCGCT.1" "AAACATACCTGGTA.1" ...
#  .. .. .. .. .. ..$ : chr [1:30] "PC_1" "PC_2" "PC_3" "PC_4" ...
#  .. .. .. ..@ feature.loadings          : num [1:2000, 1:30] -0.000166 0.000153 -0.001022 0.045661 0.059561 ...
#  .. .. .. .. ..- attr(*, "dimnames")=List of 2
#  .. .. .. .. .. ..$ : chr [1:2000] "HBB" "HBA2" "HBA1" "CCL4" ...
#  .. .. .. .. .. ..$ : chr [1:30] "PC_1" "PC_2" "PC_3" "PC_4" ...
#  .. .. .. ..@ feature.loadings.projected: num[0 , 0 ] 
#  .. .. .. ..@ assay.used                : chr "integrated"
#  .. .. .. ..@ global                    : logi FALSE
#  .. .. .. ..@ stdev                     : num [1:30] 9.46 3.87 3.49 3.47 2.85 ...
#  .. .. .. ..@ key                       : chr "PC_"
#  .. .. .. ..@ jackstraw                 :Formal class 'JackStrawData' [package "Seurat"] with 4 slots
#  .. .. .. .. .. ..@ empirical.p.values     : num[0 , 0 ] 
#  .. .. .. .. .. ..@ fake.reduction.scores  : num[0 , 0 ] 
#  .. .. .. .. .. ..@ empirical.p.values.full: num[0 , 0 ] 
#  .. .. .. .. .. ..@ overall.p.values       : num[0 , 0 ] 
#  .. .. .. ..@ misc                      :List of 1
#  .. .. .. .. ..$ total.variance: num 1417
#  .. ..$ umap:Formal class 'DimReduc' [package "Seurat"] with 9 slots
#  .. .. .. ..@ cell.embeddings           : num [1:13999, 1:2] -8.37 -7.95 -7.84 -3.62 7.74 ...
#  .. .. .. .. ..- attr(*, "scaled:center")= num [1:2] 0.368 -0.322
#  .. .. .. .. ..- attr(*, "dimnames")=List of 2
#  .. .. .. .. .. ..$ : chr [1:13999] "AAACATACATTTCC.1" "AAACATACCAGAAA.1" "AAACATACCTCGCT.1" "AAACATACCTGGTA.1" ...
#  .. .. .. .. .. ..$ : chr [1:2] "UMAP_1" "UMAP_2"
#  .. .. .. ..@ feature.loadings          : num[0 , 0 ] 
#  .. .. .. ..@ feature.loadings.projected: num[0 , 0 ] 
#  .. .. .. ..@ assay.used                : chr "integrated"
#  .. .. .. ..@ global                    : logi TRUE
#  .. .. .. ..@ stdev                     : num(0) 
#  .. .. .. ..@ key                       : chr "UMAP_"
#  .. .. .. ..@ jackstraw                 :Formal class 'JackStrawData' [package "Seurat"] with 4 slots
#  .. .. .. .. .. ..@ empirical.p.values     : num[0 , 0 ] 
#  .. .. .. .. .. ..@ fake.reduction.scores  : num[0 , 0 ] 
#  .. .. .. .. .. ..@ empirical.p.values.full: num[0 , 0 ] 
#  .. .. .. .. .. ..@ overall.p.values       : num[0 , 0 ] 
#  .. .. .. ..@ misc                      : list()
#  ..@ project.name: chr "SeuratProject"
#  ..@ misc        : list()
#  ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
#  .. ..$ : int [1:3] 3 1 2
#  ..@ commands    :List of 7
#  .. ..$ FindIntegrationAnchors      :Formal class 'SeuratCommand' [package "Seurat"] with 5 slots
#  .. .. .. ..@ name       : chr "FindIntegrationAnchors"
#  .. .. .. ..@ time.stamp : POSIXct[1:1], format: "2020-01-30 12:04:53"
#  .. .. .. ..@ assay.used : Named chr [1:2] "RNA" "RNA"
#  .. .. .. .. ..- attr(*, "names")= chr [1:2] "CTRL" "STIM"
#  .. .. .. ..@ call.string: chr "FindIntegrationAnchors(object.list = ifnb.list, dims = 1:20)"
#  .. .. .. ..@ params     :List of 14
#  .. .. .. .. ..$ assay               : Named chr [1:2] "RNA" "RNA"
#  .. .. .. .. .. ..- attr(*, "names")= chr [1:2] "CTRL" "STIM"
#  .. .. .. .. ..$ anchor.features     : chr [1:2000] "HBB" "HBA2" "HBA1" "CCL4" ...
#  .. .. .. .. ..$ scale               : logi TRUE
#  .. .. .. .. ..$ normalization.method: chr "LogNormalize"
#  .. .. .. .. ..$ reduction           : chr "cca"
#  .. .. .. .. ..$ l2.norm             : logi TRUE
#  .. .. .. .. ..$ dims                : int [1:20] 1 2 3 4 5 6 7 8 9 10 ...
#  .. .. .. .. ..$ k.anchor            : num 5
#  .. .. .. .. ..$ k.filter            : num 200
#  .. .. .. .. ..$ k.score             : num 30
#  .. .. .. .. ..$ max.features        : num 200
#  .. .. .. .. ..$ nn.method           : chr "rann"
#  .. .. .. .. ..$ eps                 : num 0
#  .. .. .. .. ..$ verbose             : logi TRUE
#  .. ..$ PairwiseIntegrateReference  :Formal class 'SeuratCommand' [package "Seurat"] with 5 slots
#  .. .. .. ..@ name       : chr "PairwiseIntegrateReference"
#  .. .. .. ..@ time.stamp : POSIXct[1:1], format: "2020-01-30 12:05:03"
#  .. .. .. ..@ assay.used : NULL
#  .. .. .. ..@ call.string: chr [1:6] "PairwiseIntegrateReference(anchorset = anchorset, new.assay.name = new.assay.name,
# " "    normalization.method = normalization.method, features = features, 
# " "    features.to.integrate = features.to.integrate, dims = dims, 
# " "    k.weight = k.weight, weight.reduction = weight.reduction, " ...
#  .. .. .. ..@ params     :List of 12
#  .. .. .. .. ..$ new.assay.name       : chr "integrated"
#  .. .. .. .. ..$ normalization.method : chr "LogNormalize"
#  .. .. .. .. ..$ features             : chr [1:2000] "HBB" "HBA2" "HBA1" "CCL4" ...
#  .. .. .. .. ..$ features.to.integrate: chr [1:2000] "HBB" "HBA2" "HBA1" "CCL4" ...
#  .. .. .. .. ..$ dims                 : int [1:20] 1 2 3 4 5 6 7 8 9 10 ...
#  .. .. .. .. ..$ k.weight             : num 100
#  .. .. .. .. ..$ sd.weight            : num 1
#  .. .. .. .. ..$ sample.tree          : num [1, 1:2] -2 -1
#  .. .. .. .. ..$ preserve.order       : logi FALSE
#  .. .. .. .. ..$ do.cpp               : logi TRUE
#  .. .. .. .. ..$ eps                  : num 0
#  .. .. .. .. ..$ verbose              : logi TRUE
#  .. ..$ ScaleData.integrated        :Formal class 'SeuratCommand' [package "Seurat"] with 5 slots
# .. .. .. ..@ name       : chr "ScaleData.integrated"
#  .. .. .. ..@ time.stamp : POSIXct[1:1], format: "2020-01-30 12:09:27"
#  .. .. .. ..@ assay.used : chr "integrated"
#  .. .. .. ..@ call.string: chr "ScaleData(samples.combined.list.combined, verbose = FALSE)"
#  .. .. .. ..@ params     :List of 10
#  .. .. .. .. ..$ features          : chr [1:2000] "HBB" "HBA2" "HBA1" "CCL4" ...
#  .. .. .. .. ..$ assay             : chr "integrated"
#  .. .. .. .. ..$ model.use         : chr "linear"
#  .. .. .. .. ..$ use.umi           : logi FALSE
#  .. .. .. .. ..$ do.scale          : logi TRUE
#  .. .. .. .. ..$ do.center         : logi TRUE
#  .. .. .. .. ..$ scale.max         : num 10
#  .. .. .. .. ..$ block.size        : num 1000
#  .. .. .. .. ..$ min.cells.to.block: num 3000
#  .. .. .. .. ..$ verbose           : logi FALSE
#  .. ..$ RunPCA.integrated           :Formal class 'SeuratCommand' [package "Seurat"] with 5 slots
#  .. .. .. ..@ name       : chr "RunPCA.integrated"
#  .. .. .. ..@ time.stamp : POSIXct[1:1], format: "2020-01-30 12:09:35"
#  .. .. .. ..@ assay.used : chr "integrated"
#  .. .. .. ..@ call.string: chr "RunPCA(samples.combined.list.combined, npcs = 30, verbose = FALSE)"
#  .. .. .. ..@ params     :List of 10
#  .. .. .. .. ..$ assay          : chr "integrated"
#  .. .. .. .. ..$ npcs           : num 30
#  .. .. .. .. ..$ rev.pca        : logi FALSE
# .. .. .. .. ..$ weight.by.var  : logi TRUE
#  .. .. .. .. ..$ verbose        : logi FALSE
#  .. .. .. .. ..$ ndims.print    : int [1:5] 1 2 3 4 5
#  .. .. .. .. ..$ nfeatures.print: num 30
#  .. .. .. .. ..$ reduction.name : chr "pca"
#  .. .. .. .. ..$ reduction.key  : chr "PC_"
#  .. .. .. .. ..$ seed.use       : num 42
#  .. ..$ RunUMAP.integrated.pca      :Formal class 'SeuratCommand' [package "Seurat"] with 5 slots
#  .. .. .. ..@ name       : chr "RunUMAP.integrated.pca"
#  .. .. .. ..@ time.stamp : POSIXct[1:1], format: "2020-01-30 12:09:55"
#  .. .. .. ..@ assay.used : chr "integrated"
#  .. .. .. ..@ call.string: chr "RunUMAP(samples.combined.list.combined, reduction = \"pca\", dims = 1:20)"
#  .. .. .. ..@ params     :List of 20
#  .. .. .. .. ..$ dims                : int [1:20] 1 2 3 4 5 6 7 8 9 10 ...
#  .. .. .. .. ..$ reduction           : chr "pca"
#  .. .. .. .. ..$ assay               : chr "integrated"
#  .. .. .. .. ..$ umap.method         : chr "uwot"
#  .. .. .. .. ..$ n.neighbors         : int 30
#  .. .. .. .. ..$ n.components        : int 2
#  .. .. .. .. ..$ metric              : chr "cosine"
#  .. .. .. .. ..$ learning.rate       : num 1
#  .. .. .. .. ..$ min.dist            : num 0.3
#  .. .. .. .. ..$ spread              : num 1
#  .. .. .. .. ..$ set.op.mix.ratio    : num 1
#  .. .. .. .. ..$ local.connectivity  : int 1
#  .. .. .. .. ..$ repulsion.strength  : num 1
#  .. .. .. .. ..$ negative.sample.rate: int 5
#  .. .. .. .. ..$ uwot.sgd            : logi FALSE
#  .. .. .. .. ..$ seed.use            : int 42
#  .. .. .. .. ..$ angular.rp.forest   : logi FALSE
#  .. .. .. .. ..$ verbose             : logi TRUE
#  .. .. .. .. ..$ reduction.name      : chr "umap"
#  .. .. .. .. ..$ reduction.key       : chr "UMAP_"
#  .. ..$ FindNeighbors.integrated.pca:Formal class 'SeuratCommand' [package "Seurat"] with 5 slots
#  .. .. .. ..@ name       : chr "FindNeighbors.integrated.pca"
#  .. .. .. ..@ time.stamp : POSIXct[1:1], format: "2020-01-30 12:09:58"
#  .. .. .. ..@ assay.used : chr "integrated"
#  .. .. .. ..@ call.string: chr "FindNeighbors(samples.combined.list.combined, reduction = \"pca\", dims = 1:20)"
#  .. .. .. ..@ params     :List of 13
#  .. .. .. .. ..$ reduction   : chr "pca"
#  .. .. .. .. ..$ dims        : int [1:20] 1 2 3 4 5 6 7 8 9 10 ...
#  .. .. .. .. ..$ assay       : chr "integrated"
#  .. .. .. .. ..$ k.param     : num 20
#  .. .. .. .. ..$ compute.SNN : logi TRUE
#  .. .. .. .. ..$ prune.SNN   : num 0.0667
#  .. .. .. .. ..$ nn.method   : chr "rann"
#  .. .. .. .. ..$ annoy.metric: chr "euclidean"
#  .. .. .. .. ..$ nn.eps      : num 0
#  .. .. .. .. ..$ verbose     : logi TRUE
#  .. .. .. .. ..$ force.recalc: logi FALSE
#  .. .. .. .. ..$ do.plot     : logi FALSE
#  .. .. .. .. ..$ graph.name  : chr [1:2] "integrated_nn" "integrated_snn"
#  .. ..$ FindClusters                :Formal class 'SeuratCommand' [package "Seurat"] with 5 slots
#  .. .. .. ..@ name       : chr "FindClusters"
#  .. .. .. ..@ time.stamp : POSIXct[1:1], format: "2020-01-30 12:10:00"
#  .. .. .. ..@ assay.used : chr "integrated"
#  .. .. .. ..@ call.string: chr "FindClusters(samples.combined.list.combined, resolution = 0.5)"
#  .. .. .. ..@ params     :List of 10
#  .. .. .. .. ..$ graph.name      : chr "integrated_snn"
#  .. .. .. .. ..$ modularity.fxn  : num 1
#  .. .. .. .. ..$ resolution      : num 0.5
#  .. .. .. .. ..$ method          : chr "matrix"
#  .. .. .. .. ..$ algorithm       : num 1
#  .. .. .. .. ..$ n.start         : num 10
#  .. .. .. .. ..$ n.iter          : num 10
#  .. .. .. .. ..$ random.seed     : num 0
#  .. .. .. .. ..$ group.singletons: logi TRUE
#  .. .. .. .. ..$ verbose         : logi TRUE
#  ..@ tools       :List of 1
#  .. ..$ Integration:Formal class 'IntegrationData' [package "Seurat"] with 7 slots
#  .. .. .. ..@ neighbors         : NULL
#  .. .. .. ..@ weights           : NULL
#  .. .. .. ..@ integration.matrix: NULL
#  .. .. .. ..@ anchors           :'data.frame':	12770 obs. of  5 variables:
#  .. .. .. .. ..$ cell1   : num [1:12770] 4860 6213 1558 99 477 ...
#  .. .. .. .. ..$ cell2   : num [1:12770] 1 1 2 6 6 8 8 9 9 11 ...
#  .. .. .. .. ..$ score   : num [1:12770] 0.281 0.562 1 0.469 0.312 ...
#  .. .. .. .. ..$ dataset1: int [1:12770] 1 1 1 1 1 1 1 1 1 1 ...
#  .. .. .. .. ..$ dataset2: int [1:12770] 2 2 2 2 2 2 2 2 2 2 ...
#  .. .. .. ..@ offsets           : NULL
#  .. .. .. ..@ objects.ncell     : NULL
#  .. .. .. ..@ sample.tree       : num [1, 1:2] -2 -1

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
#################################################################################################################################### DIFFERENTIAL EXPRESSED GENES
#################################################################################################################################### across CONDITIONS
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# head(samples.combined.list.combined@meta.data)

# theme_set(theme_cowplot())
# t.cells <- subset(samples.combined.list.combined, idents = "CD4 Naive T")
# Idents(t.cells) <- "orig.ident"
# avg.t.cells <- log1p(AverageExpression(t.cells, verbose = FALSE)$RNA)
# avg.t.cells$gene <- rownames(avg.t.cells)

# cd14.mono <- subset(samples.combined.list.combined, idents = "CD14 Mono")
# Idents(cd14.mono) <- "orig.ident"
# avg.cd14.mono <- log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA)
# avg.cd14.mono$gene <- rownames(avg.cd14.mono)

# genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")

# p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
# p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)

# p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
# p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)

# plot_grid(p1, p2)
# ggsave(paste(NAME, "GRID.Plot", "T_CELLS", "MONOCYTES", "png", sep="."), width = 40, height = 20, units = "cm")

# Because we are confident in having identified common cell types across condition, we can ask what genes change in different conditions for cells of the same type. 
# First, we create a column in the meta.data slot to hold both the cell type and stimulation information and switch the current ident to that column. 
# Then we use FindMarkers to find the genes that are different between stimulated and control B cells. Notice that many of the top genes that show up are the same 
# as the ones we plotted earlier as core interferon response genes. 
# Additionally, genes like CXCL10 which we saw were specific to monocyte and B cell interferon response show up as highly significant in this list as well.

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

head(samples.combined.list.combined@meta.data)

# we may check initially if the order is the same : using : all.equal()

all.equal( rownames(as.data.frame(samples.combined.list.combined@meta.data)), rownames(as.data.frame(Idents(samples.combined.list.combined))) ) # TRUE

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

samples.combined.list.combined@meta.data$celltype.stim <- paste0(Idents(samples.combined.list.combined), "_", samples.combined.list.combined@meta.data$orig.ident)

head(samples.combined.list.combined@meta.data)
tail(samples.combined.list.combined@meta.data)

# samples.combined.list.combined <- StashIdent(samples.combined.list.combined, save.name = "stim.experiment")

# With Seurat 3.X, stashing identity classes can be accomplished with the following:
# samples.combined.list.combined[["stim.experiment"]] <- Idents(object = samples.combined.list.combined)

# samples.combined.list.combined <- SetIdent(samples.combined.list.combined, id = "celltype.stim")   # POINTING to the NEW COLUMN that we have just made !!!!!!

samples.combined.list.combined$celltype <- Idents(samples.combined.list.combined)

Idents(samples.combined.list.combined) <- "celltype.stim"

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################

write.table(as.data.frame(samples.combined.list.combined@meta.data), 
                        file=paste(NAME, "META.DATA", "in.assay.RNA.at.the.midpoint.SCRIPT.txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

#####################################################################################################################################
##################################################################################################################################### in CONOS 
##################################################################################################################################### INDEXING STARTS at 1 
#####################################################################################################################################
# DefaultAssay(samples.combined.list.combined) <- "RNA"
# they do start CLUSTERING from 1 in CONOS 

# LIST.CLUSTERS.and.DIFFERENTIAL.MARKERS <- list()

# CLUSTERS = as.numeric(unique(samples.combined.list.combined@meta.data$celltype))

# for (i in 1:length(CLUSTERS)){
   
   # i=i-1

#   if ( ( TABLE.CTRL[TABLE.CTRL$Var2 == i,]$Freq >=3) && (TABLE.STIM[TABLE.STIM$Var2 == i,]$Freq >=3) )
#   {

#      IDENT1 = paste0(i, "_", CTRL)
#      IDENT2 = paste0(i, "_", STIM)

#      LIST.CLUSTERS.and.DIFFERENTIAL.MARKERS[[i]] <- FindMarkers(samples.combined.list.combined, ident.1 = IDENT1, 
#                                                                                                 ident.2 = IDENT2, 
#                                                                                                 print.bar = FALSE, only.pos = FALSE)

   ### printing all DIFFERENTIAL MARKERS in a CELL CLUSTER
#   x <- as.data.frame(as.matrix(LIST.CLUSTERS.and.DIFFERENTIAL.MARKERS[[i]]))  
#   x$gene <- row.names(x)

#   write.table(x, file=paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS.CLUSTER", i, "LIST.txt", sep="."), 
#                  sep="\t", quote=F, row.names=T, col.names=T)
    
   ######################################################################################################### genes 1-12 or MIN
   ### VlnPlot(object = samples.combined.list.combined, features = gene, split.by = "celltype.stim",  split.plot = FALSE)

#   for (j in 1:min(dim(x)[1], 12))           #################################### printing top 12 MARKERS in each CLUSTER or the MIN 
#   {
#        gene <- x$gene[j] 
#        print(gene)   

        ### VLN PLOT
        ### VlnPlot(object = samples.combined.list.combined, features = x$gene[1:6])
        ### VlnPlot(object = samples.combined.list.combined, features = gene, group.by = "orig.ident", split.by = "orig.ident",  split.plot = TRUE)  
        ### VlnPlot(object = samples.combined.list.combined, features = gene, group.by = "orig.ident",  split.plot = FALSE) 
        
#        VlnPlot(object = samples.combined.list.combined, features = gene, split.by = "orig.ident",  split.plot = FALSE)
#        ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS.CLUSTER", i, "top.gene", j,  gene, "VIOLIN.PLOT1.png", sep="."), 
#               width=60, height=20, units = "cm")
  
#        VlnPlot(object = samples.combined.list.combined, features = gene, split.by = "celltype.stim",  split.plot = FALSE) 
#        ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS.CLUSTER", i, "top.gene", j,  gene, "VIOLIN.PLOT2.png", sep="."), 
#               width=60, height=20, units = "cm")
         
#        VlnPlot(object = samples.combined.list.combined, features = gene, group.by = "celltype.stim",  split.plot = FALSE) 
#        ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS.CLUSTER", i, "top.gene", j,  gene, "VIOLIN.PLOT3.png", sep="."), 
#               width=60, height=20, units = "cm")

        ### FEATURE PLOT
        ### FeaturePlot(object = samples.combined.list.combined, features = x$gene[1:6] , cols = c("lightgrey", "blue"), reduction = "tsne")   
        # FeaturePlot(object = samples.combined.list.combined, features = gene, cols = c("lightgrey", "blue"), reduction = "tsne", split.by = "orig.ident")        
        # ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS.CLUSTER", i, "top.gene", j, gene, "FEATURE.TSNE.png", sep="."), width=60, height=40, units = "cm")

        ### FEATURE PLOT 
        ### FeaturePlot(object = samples.combined.list.combined, features = x$gene[1:6] , cols = c("lightgrey", "blue"), reduction = "umap")    
        # FeaturePlot(object = samples.combined.list.combined, features = gene, cols = c("lightgrey", "blue"), reduction = "umap", split.by = "orig.ident")        
        # ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS.CLUSTER", i, "top.gene", j, gene, "FEATURE.UMAP.png", sep="."), width=60, height=40, units = "cm")

#        FeaturePlot(object = samples.combined.list.combined, features = gene, cols = c("lightgrey", "blue"), reduction = "largeVis", split.by = "orig.ident")        
#        ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS.CLUSTER", i, "top.gene", j, gene, "FEATURE.CONOS.png", sep="."), 
#               width=60, height=40, units = "cm")

        ### DOT PLOT
        ### DotPlot(object = samples.combined.list.combined, features = x$gene[1:6])
#        DotPlot(object = samples.combined.list.combined, features = gene,  split.by = "orig.ident")
#        ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS.CLUSTER", i, "top.gene", j, gene, "DOT.PLOT.png", sep="."), 
#               width=20, height=20, units = "cm") 

        ### RIDGE PLOT 
        ### RidgePlot(object = samples.combined.list.combined, features = x$gene[1:6])
#        RidgePlot(object = samples.combined.list.combined, features = gene, group.by = "orig.ident")
#        ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS.CLUSTER", i, "top.gene", j, gene, "RIDGE.PLOT.png", sep="."), 
#               width=40, height=20, units = "cm") 
#     }
#  }
#}

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
##################################################################################################################################### displaying the KNOWN MARKERS
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################

# library(Hmisc)
# LIST.KNOWN.MARKERS$MARKER = capitalize(LIST.KNOWN.MARKERS$MARKER)

DefaultAssay(samples.combined.list.combined) <- "RNA"

# LIST.KNOWN.MARKERS = read.delim("a.LIST.MARKERS.RETINA.SUI.mouse.txt", header=T, sep="\t", stringsAsFactors=F)
# dim(LIST.KNOWN.MARKERS)
# head(LIST.KNOWN.MARKERS)

for (i in 1:dim(LIST.KNOWN.MARKERS)[1]){
    if (LIST.KNOWN.MARKERS$MARKER[i] %in% row.names(samples.combined.list.combined[["RNA"]]@data)) {

        RidgePlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i])
        ggsave(paste(NAME, "figure9", "RESOLUTION", RESOLUTION, "KNOWN.MARKERS", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "RIDGE.png", sep="."), 
                     width=60, height=100, units = "cm")
        
        # FeaturePlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i], cols = c("lightgrey", "blue"), reduction = "tsne")
        # ggsave(paste(NAME, "figure7", "RESOLUTION", RESOLUTION, "KNOWN.markers", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "TSNE.png", sep="."))

        # FeaturePlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i], cols = c("lightgrey", "blue"), reduction = "umap")
        # ggsave(paste(NAME, "figure7", "RESOLUTION", RESOLUTION, "KNOWN.markers", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "UMAP.png", sep="."))

        FeaturePlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i], cols = c("lightgrey", "blue"), reduction = "largeVis")
        ggsave(paste(NAME, "figure9", "RESOLUTION", RESOLUTION, "KNOWN.MARKERS", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "largeVis.png", sep="."), 
               width=20, height=20, units = "cm")

        VlnPlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i])
        ggsave(paste(NAME, "figure9", "RESOLUTION", RESOLUTION, "KNOWN.MARKERS", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "VLN.png", sep="."), 
               width=120, height=40, units = "cm")
        
    }
}


for (i in 1:dim(LIST.KNOWN.MARKERS)[1]){
    if (LIST.KNOWN.MARKERS$MARKER[i] %in% row.names(samples.combined.list.combined[["RNA"]]@data)) {

        # RidgePlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i], split.by = "orig.ident")
        # ggsave(paste(NAME, "figure7", "RESOLUTION", RESOLUTION, "KNOWN.markers", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "v2.RIDGE.png", sep="."))
        
        # FeaturePlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i], split.by = "orig.ident", cols = c("lightgrey", "blue"), reduction = "tsne")
        # ggsave(paste(NAME, "figure7", "RESOLUTION", RESOLUTION, "KNOWN.markers", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "v2.TSNE.png", sep="."))

        # FeaturePlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i], split.by = "orig.ident", cols = c("lightgrey", "blue"), reduction = "umap")
        # ggsave(paste(NAME, "figure7", "RESOLUTION", RESOLUTION, "KNOWN.markers", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "v2.UMAP.png", sep="."))

        FeaturePlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i], split.by = "orig.ident", cols = c("lightgrey", "blue"), reduction = "largeVis")
        ggsave(paste(NAME, "figure9", "RESOLUTION", RESOLUTION, "KNOWN.MARKERS", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "v2.largeVis.png", sep="."), 
               width=60, height=20, units = "cm")

        VlnPlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i], split.by = "orig.ident")
        ggsave(paste(NAME, "figure9", "RESOLUTION", RESOLUTION, "KNOWN.MARKERS", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "v2.VLN.png", sep="."), 
               width=120, height=20, units = "cm")
        
    }
}

#####################################################################################################
##################################################################################################### for HUMAN GENE SYMBOLS
### using the UPPER CASES of the HUGO NAMES

LIST.KNOWN.MARKERS$MARKER = toupper(LIST.KNOWN.MARKERS$MARKER)

for (i in 1:dim(LIST.KNOWN.MARKERS)[1]){
    if (LIST.KNOWN.MARKERS$MARKER[i] %in% row.names(samples.combined.list.combined[["RNA"]]@data)) {

        RidgePlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i])
        ggsave(paste(NAME, "figure9", "RESOLUTION", RESOLUTION, "KNOWN.MARKERS", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "RIDGE.png", sep="."), 
                     width=60, height=100, units = "cm")
        
        # FeaturePlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i], cols = c("lightgrey", "blue"), reduction = "tsne")
        # ggsave(paste(NAME, "figure7", "RESOLUTION", RESOLUTION, "KNOWN.markers", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "TSNE.png", sep="."))

        # FeaturePlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i], cols = c("lightgrey", "blue"), reduction = "umap")
        # ggsave(paste(NAME, "figure7", "RESOLUTION", RESOLUTION, "KNOWN.markers", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "UMAP.png", sep="."))

        FeaturePlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i], cols = c("lightgrey", "blue"), reduction = "largeVis")
        ggsave(paste(NAME, "figure9", "RESOLUTION", RESOLUTION, "KNOWN.MARKERS", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "largeVis.png", sep="."), 
               width=20, height=20, units = "cm")

        VlnPlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i])
        ggsave(paste(NAME, "figure9", "RESOLUTION", RESOLUTION, "KNOWN.MARKERS", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "VLN.png", sep="."), 
               width=120, height=40, units = "cm")
        
    }
}


for (i in 1:dim(LIST.KNOWN.MARKERS)[1]){
    if (LIST.KNOWN.MARKERS$MARKER[i] %in% row.names(samples.combined.list.combined[["RNA"]]@data)) {

        # RidgePlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i], split.by = "orig.ident")
        # ggsave(paste(NAME, "figure7", "RESOLUTION", RESOLUTION, "KNOWN.markers", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "RIDGE.png", sep="."))
        
        # FeaturePlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i], split.by = "orig.ident", cols = c("lightgrey", "blue"), reduction = "tsne")
        # ggsave(paste(NAME, "figure7", "RESOLUTION", RESOLUTION, "KNOWN.markers", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "TSNE.png", sep="."))

        # FeaturePlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i], split.by = "orig.ident", cols = c("lightgrey", "blue"), reduction = "umap")
        # ggsave(paste(NAME, "figure7", "RESOLUTION", RESOLUTION, "KNOWN.markers", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "UMAP.png", sep="."))

        FeaturePlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i], split.by = "orig.ident", cols = c("lightgrey", "blue"), reduction = "largeVis")
        ggsave(paste(NAME, "figure9", "RESOLUTION", RESOLUTION, "KNOWN.MARKERS", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "v2.largeVis.png", sep="."), 
               width=60, height=20, units = "cm")

        VlnPlot(object = samples.combined.list.combined, features = LIST.KNOWN.MARKERS$MARKER[i], split.by = "orig.ident")
        ggsave(paste(NAME, "figure9", "RESOLUTION", RESOLUTION, "KNOWN.MARKERS", "cell", LIST.KNOWN.MARKERS$CELL[i], "gene", LIST.KNOWN.MARKERS$MARKER[i], "v2.VLN.png", sep="."), 
               width=120, height=20, units = "cm")
        
    }
}

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
##################################################################################################################################### to display the CLUSTER_SPECIFIC MARKERS
#####################################################################################################################################
#####################################################################################################################################
##################################################################################################################################### CLUSTERING_SPECIFIC MARKERS
##################################################################################################################################### 
##################################################################################################################################### 
##################################################################################################################################### in CONOS
##################################################################################################################################### the INDEXING STARTS at 1 ...
# WT and STZ conditions :
# they do start CLUSTERING from 1 in CONOS 

#####################################################################################################################################
#####################################################################################################################################
##################################################################################################################################### the old code is below 

# for (i in 1:length(samples.combined.list.combined.markers.clusters.all.list.markers))     
# {
    
     # print(i) 
     # LENGTH_LIST = dim(samples.combined.list.combined.markers.clusters.all.list.markers[[i]])[1]  
     # printing the MINIMUM between LENGTH_LIST, and 12 (in such a way that the lists shorter than 12 are still printed)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#     print(i) 
#     if (dim(samples.combined.list.combined.markers.clusters.all.list.markers[[i]])[1] > 0) {  
#
#         for (j in 1:min(dim(samples.combined.list.combined.markers.clusters.all.list.markers[[i]])[1], 12) )    ##### printing top 12 MARKERS in each CLUSTER or the MIN 
#         {
#          gene <- samples.combined.list.combined.markers.clusters.all.list.markers[[i]]$gene[j] 
#          print(gene)  

        # VIOLIN PLOT 
        # png(paste("figure", "cluster", i, "top.gene", j, "plot.VIOLIN.for", gene, "png", sep="."),  width=20, height=20, units="cm", res=400)
#         VlnPlot(object = samples.combined.list.combined, features = gene, split.by = "orig.ident")
#         ggsave(paste(NAME, "figure8.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.VIOLIN.for", gene, "png", sep="."), 
#                 width=40, height=20, units="cm")
        # dev.off()

        # FEATURE PLOT
        # png(paste("figure", "cluster", i, "top.gene", j, "plot.FEATURE.for", gene, "png", sep="."),  width=10, height=10, units="cm", res=400)
        #  FeaturePlot(object = samples.combined.list.combined, features = gene , cols = c("lightgrey", "blue"), reduction = "tsne", split.by = "orig.ident", label = TRUE)        
        #  ggsave(paste(NAME, "figure8.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.FEATURE.TSNE.for", gene, "png", sep="."), 
        #         width=40, height=20, units="cm")
        # dev.off()

        # FEATURE PLOT
        # png(paste("figure", "cluster", i, "top.gene", j, "plot.FEATURE.for", gene, "png", sep="."),  width=10, height=10, units="cm", res=400)
        #  FeaturePlot(object = samples.combined.list.combined, features = gene , cols = c("lightgrey", "blue"), reduction = "umap", split.by = "orig.ident", label = TRUE)  
        #  ggsave(paste(NAME, "figure8.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.FEATURE.UMAP.for", gene, "png", sep="."), 
        #         width=40, height=20, units="cm")      
        # dev.off()

#         FeaturePlot(object = samples.combined.list.combined, features = gene , cols = c("lightgrey", "blue"), reduction = "largeVis", split.by = "orig.ident", label = TRUE)  
#         ggsave(paste(NAME, "figure8.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.FEATURE.largeVis.for", gene, "png", sep="."), 
#               width=40, height=20, units="cm")     

        # DOT PLOT
        # png(paste("figure", "cluster", i, "top.gene", j, "plot.DOT.for", gene, "png", sep="."),  width=30, height=30, units="cm", res=400)
#          DotPlot(object = samples.combined.list.combined, features = gene, split.by = "orig.ident")
#          ggsave(paste(NAME, "figure8.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.DOT", gene, "png", sep="."), 
#                 width=20, height=20, units="cm") 
        # dev.off()

        # RIDGE PLOT 
        # png(paste("figure", "cluster", i, "top.gene", j, "plot.RIDGE.for", gene, "png", sep="."),  width=30, height=30, units="cm", res=400)
#          RidgePlot(object = samples.combined.list.combined, features = gene)
#          ggsave(paste(NAME, "figure8.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.RIDGE", gene, "png", sep="."), 
#                 width=20, height=20, units="cm") 
        # dev.off()
#       }
#   }
#}

#####################################################################################################################################
#####################################################################################################################################
##################################################################################################################################### re-display of CLUSTER SPECIFIC MARKERS

for (i in 1:length(samples.combined.list.combined.markers.clusters.all.list.markers))     
{
    
     # print(i) 
     # LENGTH_LIST = dim(samples.combined.list.combined.markers.clusters.all.list.markers[[i]])[1]  
     # printing the MINIMUM between LENGTH_LIST, and 12 (in such a way that the lists shorter than 12 are still printed)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     print(i) 
     if (dim(samples.combined.list.combined.markers.clusters.all.list.markers[[i]])[1] > 0) {  

         for (j in 1:min(dim(samples.combined.list.combined.markers.clusters.all.list.markers[[i]])[1], 3) )    ##### printing top 3 MARKERS in each CLUSTER or the MIN 
         {
          gene <- samples.combined.list.combined.markers.clusters.all.list.markers[[i]]$gene[j] 
          print(gene)  

        # VIOLIN PLOT 
        # png(paste("figure", "cluster", i, "top.gene", j, "plot.VIOLIN.for", gene, "png", sep="."),  width=20, height=20, units="cm", res=400)
          VlnPlot(object = samples.combined.list.combined, features = gene, split.by = "orig.ident")
          ggsave(paste(NAME, "figure8.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.VIOLIN.for", gene, "png", sep="."), 
                 width=80, height=20, units="cm")
        # dev.off()

        # FEATURE PLOT
        # png(paste("figure", "cluster", i, "top.gene", j, "plot.FEATURE.for", gene, "png", sep="."),  width=10, height=10, units="cm", res=400)
        #  FeaturePlot(object = samples.combined.list.combined, features = gene , cols = c("lightgrey", "blue"), reduction = "tsne", split.by = "orig.ident", label = TRUE)        
        #  ggsave(paste(NAME, "figure8.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.FEATURE.TSNE.for", gene, "png", sep="."), 
        #         width=40, height=20, units="cm")
        # dev.off()

        # FEATURE PLOT
        # png(paste("figure", "cluster", i, "top.gene", j, "plot.FEATURE.for", gene, "png", sep="."),  width=10, height=10, units="cm", res=400)
        #  FeaturePlot(object = samples.combined.list.combined, features = gene , cols = c("lightgrey", "blue"), reduction = "umap", split.by = "orig.ident", label = TRUE)  
        #  ggsave(paste(NAME, "figure8.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.FEATURE.UMAP.for", gene, "png", sep="."), 
        #         width=40, height=20, units="cm")      
        # dev.off()

         FeaturePlot(object = samples.combined.list.combined, features = gene , cols = c("lightgrey", "blue"), reduction = "largeVis", split.by = "orig.ident", label = TRUE)  
         ggsave(paste(NAME, "figure8.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.FEATURE.largeVis.for", gene, "png", sep="."), 
               width=60, height=20, units="cm")     

        # DOT PLOT
        # png(paste("figure", "cluster", i, "top.gene", j, "plot.DOT.for", gene, "png", sep="."),  width=30, height=30, units="cm", res=400)
        #  DotPlot(object = samples.combined.list.combined, features = gene, split.by = "orig.ident")
        #  ggsave(paste(NAME, "figure8.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.DOT", gene, "png", sep="."), 
        #         width=20, height=20, units="cm") 
        # dev.off()

        # RIDGE PLOT 
        # png(paste("figure", "cluster", i, "top.gene", j, "plot.RIDGE.for", gene, "png", sep="."),  width=30, height=30, units="cm", res=400)
          RidgePlot(object = samples.combined.list.combined, features = gene)
          ggsave(paste(NAME, "figure8.samples.combined", "RESOLUTION", RESOLUTION, "CLUSTER", i,  "top.gene", j, "plot.RIDGE", gene, "png", sep="."), 
                 width=80, height=80, units="cm") 
        # dev.off()
       }
   }
}

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################  DIFFERENTIAL ANALYSIS

# STZ1 versus STZ3
# STZ3 versus STZ5
# STZ1 versus STZ5

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
##################################################################################################################################### in CONOS 
##################################################################################################################################### INDEXING STARTS at 1 
#####################################################################################################################################
# DefaultAssay(samples.combined.list.combined) <- "RNA"
# they do start CLUSTERING from 1 in CONOS 
##################################################################################################################################### compare STIM1 versus STIM3
#####################################################################################################################################
##################################################################################################################################### aka compare STZ1 versus STZ3

LIST.CLUSTERS.and.DIFFERENTIAL.MARKERS <- list()
CLUSTERS = as.numeric(unique(samples.combined.list.combined@meta.data$celltype))

for (i in 1:length(CLUSTERS)){

#i = i-1

if ( ( TABLE.STIM1[TABLE.STIM1$Var2 == i,]$Freq >=3) && (TABLE.STIM3[TABLE.STIM3$Var2 == i,]$Freq >=3) ){


      IDENT1 = paste0(i, "_", STIM1)
      IDENT2 = paste0(i, "_", STIM3)

      LIST.CLUSTERS.and.DIFFERENTIAL.MARKERS[[i]] <- FindMarkers(samples.combined.list.combined, ident.1 = IDENT1, 
                                                                                                 ident.2 = IDENT2, 
                                                                                                 print.bar = FALSE, only.pos = FALSE)

   ### printing all DIFFERENTIAL MARKERS in a CELL CLUSTER
   x <- as.data.frame(as.matrix(LIST.CLUSTERS.and.DIFFERENTIAL.MARKERS[[i]]))  
   x$gene <- row.names(x)

   write.table(x, file=paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS",  STIM1, STIM3, "CLUSTER", i, "LIST.txt", sep="."), 
                  sep="\t", quote=F, row.names=T, col.names=T)
    
   ######################################################################################################### genes 1-12 or MIN
   ### VlnPlot(object = samples.combined.list.combined, features = gene, split.by = "celltype.stim",  split.plot = FALSE)

   for (j in 1:min(dim(x)[1], 3))           #################################### printing top 12 MARKERS in each CLUSTER or the MIN 
   {
        gene <- x$gene[j] 
        print(gene)   

        ### VLN PLOT
        ### VlnPlot(object = samples.combined.list.combined, features = x$gene[1:6])
        ### VlnPlot(object = samples.combined.list.combined, features = gene, group.by = "orig.ident", split.by = "orig.ident",  split.plot = TRUE)  
        ### VlnPlot(object = samples.combined.list.combined, features = gene, group.by = "orig.ident",  split.plot = FALSE) 
        
        VlnPlot(object = samples.combined.list.combined, features = gene, split.by = "orig.ident",  split.plot = FALSE)
        ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM1, STIM3, "CLUSTER", i, "top.gene", j,  gene, "VIOLIN.PLOT1.png", sep="."), 
               width=120, height=20, units = "cm")
  
        # VlnPlot(object = samples.combined.list.combined, features = gene, split.by = "celltype.stim",  split.plot = FALSE) 
        # ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM1, STIM3, "CLUSTER", i, "top.gene", j,  gene, "VIOLIN.PLOT2.png", sep="."), 
        #       width=60, height=20, units = "cm")
         
        # VlnPlot(object = samples.combined.list.combined, features = gene, group.by = "celltype.stim",  split.plot = FALSE) 
        # ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM1, STIM3, "CLUSTER", i, "top.gene", j,  gene, "VIOLIN.PLOT3.png", sep="."), 
        #       width=60, height=20, units = "cm")

        ### FEATURE PLOT
        ### FeaturePlot(object = samples.combined.list.combined, features = x$gene[1:6] , cols = c("lightgrey", "blue"), reduction = "tsne")   
        # FeaturePlot(object = samples.combined.list.combined, features = gene, cols = c("lightgrey", "blue"), reduction = "tsne", split.by = "orig.ident")        
        # ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS.CLUSTER", i, "top.gene", j, gene, "FEATURE.TSNE.png", sep="."), width=60, height=40, units = "cm")

        ### FEATURE PLOT 
        ### FeaturePlot(object = samples.combined.list.combined, features = x$gene[1:6] , cols = c("lightgrey", "blue"), reduction = "umap")    
        # FeaturePlot(object = samples.combined.list.combined, features = gene, cols = c("lightgrey", "blue"), reduction = "umap", split.by = "orig.ident")        
        # ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS.CLUSTER", i, "top.gene", j, gene, "FEATURE.UMAP.png", sep="."), width=60, height=40, units = "cm")

        FeaturePlot(object = samples.combined.list.combined, features = gene, cols = c("lightgrey", "blue"), reduction = "largeVis", split.by = "orig.ident")        
        ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM1, STIM3, "CLUSTER", i, "top.gene", j, gene, "FEATURE.CONOS.png", sep="."), 
               width=80, height=30, units = "cm")

        ### DOT PLOT
        ### DotPlot(object = samples.combined.list.combined, features = x$gene[1:6])
        # DotPlot(object = samples.combined.list.combined, features = gene,  split.by = "orig.ident")
        # ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM1, STIM3, "CLUSTER", i, "top.gene", j, gene, "DOT.PLOT.png", sep="."), 
        #        width=20, height=20, units = "cm") 

        ### RIDGE PLOT 
        ### RidgePlot(object = samples.combined.list.combined, features = x$gene[1:6])
        RidgePlot(object = samples.combined.list.combined, features = gene, group.by = "orig.ident")
        ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM1, STIM3, "CLUSTER", i, "top.gene", j, gene, "RIDGE.PLOT.png", sep="."), 
               width=40, height=20, units = "cm") 
     }
  }
}

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
##################################################################################################################################### in CONOS 
##################################################################################################################################### INDEXING STARTS at 1 
#####################################################################################################################################
# DefaultAssay(samples.combined.list.combined) <- "RNA"
# they do start CLUSTERING from 1 in CONOS 
##################################################################################################################################### compare STIM1 versus STIM5
#####################################################################################################################################
##################################################################################################################################### aka compare STZ1 versus STZ5

LIST.CLUSTERS.and.DIFFERENTIAL.MARKERS <- list()
CLUSTERS = as.numeric(unique(samples.combined.list.combined@meta.data$celltype))

for (i in 1:length(CLUSTERS)){

#i = i-1

if ( ( TABLE.STIM1[TABLE.STIM1$Var2 == i,]$Freq >=3) && (TABLE.STIM5[TABLE.STIM5$Var2 == i,]$Freq >=3) ){


      IDENT1 = paste0(i, "_", STIM1)
      IDENT2 = paste0(i, "_", STIM5)

      LIST.CLUSTERS.and.DIFFERENTIAL.MARKERS[[i]] <- FindMarkers(samples.combined.list.combined, ident.1 = IDENT1, 
                                                                                                 ident.2 = IDENT2, 
                                                                                                 print.bar = FALSE, only.pos = FALSE)

   ### printing all DIFFERENTIAL MARKERS in a CELL CLUSTER
   x <- as.data.frame(as.matrix(LIST.CLUSTERS.and.DIFFERENTIAL.MARKERS[[i]]))  
   x$gene <- row.names(x)

   write.table(x, file=paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS",  STIM1, STIM5, "CLUSTER", i, "LIST.txt", sep="."), 
                  sep="\t", quote=F, row.names=T, col.names=T)
    
   ######################################################################################################### genes 1-12 or MIN
   ### VlnPlot(object = samples.combined.list.combined, features = gene, split.by = "celltype.stim",  split.plot = FALSE)

   for (j in 1:min(dim(x)[1], 3))           #################################### printing top 12 MARKERS in each CLUSTER or the MIN 
   {
        gene <- x$gene[j] 
        print(gene)   

        ### VLN PLOT
        ### VlnPlot(object = samples.combined.list.combined, features = x$gene[1:6])
        ### VlnPlot(object = samples.combined.list.combined, features = gene, group.by = "orig.ident", split.by = "orig.ident",  split.plot = TRUE)  
        ### VlnPlot(object = samples.combined.list.combined, features = gene, group.by = "orig.ident",  split.plot = FALSE) 
        
        VlnPlot(object = samples.combined.list.combined, features = gene, split.by = "orig.ident",  split.plot = FALSE)
        ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM1, STIM5, "CLUSTER", i, "top.gene", j,  gene, "VIOLIN.PLOT1.png", sep="."), 
               width=120, height=20, units = "cm")
  
        # VlnPlot(object = samples.combined.list.combined, features = gene, split.by = "celltype.stim",  split.plot = FALSE) 
        # ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM1, STIM5, "CLUSTER", i, "top.gene", j,  gene, "VIOLIN.PLOT2.png", sep="."), 
        #       width=60, height=20, units = "cm")
         
        # VlnPlot(object = samples.combined.list.combined, features = gene, group.by = "celltype.stim",  split.plot = FALSE) 
        # ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM1, STIM5, "CLUSTER", i, "top.gene", j,  gene, "VIOLIN.PLOT3.png", sep="."), 
        #       width=60, height=20, units = "cm")

        ### FEATURE PLOT
        ### FeaturePlot(object = samples.combined.list.combined, features = x$gene[1:6] , cols = c("lightgrey", "blue"), reduction = "tsne")   
        # FeaturePlot(object = samples.combined.list.combined, features = gene, cols = c("lightgrey", "blue"), reduction = "tsne", split.by = "orig.ident")        
        # ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS.CLUSTER", i, "top.gene", j, gene, "FEATURE.TSNE.png", sep="."), width=60, height=40, units = "cm")

        ### FEATURE PLOT 
        ### FeaturePlot(object = samples.combined.list.combined, features = x$gene[1:6] , cols = c("lightgrey", "blue"), reduction = "umap")    
        # FeaturePlot(object = samples.combined.list.combined, features = gene, cols = c("lightgrey", "blue"), reduction = "umap", split.by = "orig.ident")        
        # ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS.CLUSTER", i, "top.gene", j, gene, "FEATURE.UMAP.png", sep="."), width=60, height=40, units = "cm")

        FeaturePlot(object = samples.combined.list.combined, features = gene, cols = c("lightgrey", "blue"), reduction = "largeVis", split.by = "orig.ident")        
        ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM1, STIM5, "CLUSTER", i, "top.gene", j, gene, "FEATURE.CONOS.png", sep="."), 
               width=80, height=30, units = "cm")

        ### DOT PLOT
        ### DotPlot(object = samples.combined.list.combined, features = x$gene[1:6])
        # DotPlot(object = samples.combined.list.combined, features = gene,  split.by = "orig.ident")
        # ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM1, STIM5, "CLUSTER", i, "top.gene", j, gene, "DOT.PLOT.png", sep="."), 
        #        width=20, height=20, units = "cm") 

        ### RIDGE PLOT 
        ### RidgePlot(object = samples.combined.list.combined, features = x$gene[1:6])
        RidgePlot(object = samples.combined.list.combined, features = gene, group.by = "orig.ident")
        ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM1, STIM5, "CLUSTER", i, "top.gene", j, gene, "RIDGE.PLOT.png", sep="."), 
               width=40, height=20, units = "cm") 
     }
  }
}

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
##################################################################################################################################### in CONOS 
##################################################################################################################################### INDEXING STARTS at 1 
#####################################################################################################################################
# DefaultAssay(samples.combined.list.combined) <- "RNA"
# they do start CLUSTERING from 1 in CONOS 
##################################################################################################################################### compare STIM3 versus STIM5
#####################################################################################################################################
##################################################################################################################################### compare STZ3 versus STZ5

LIST.CLUSTERS.and.DIFFERENTIAL.MARKERS <- list()
CLUSTERS = as.numeric(unique(samples.combined.list.combined@meta.data$celltype))

for (i in 1:length(CLUSTERS)){

#i = i-1

if ( ( TABLE.STIM3[TABLE.STIM3$Var2 == i,]$Freq >=3) && (TABLE.STIM5[TABLE.STIM5$Var2 == i,]$Freq >=3) ){


      IDENT1 = paste0(i, "_", STIM3)
      IDENT2 = paste0(i, "_", STIM5)

      LIST.CLUSTERS.and.DIFFERENTIAL.MARKERS[[i]] <- FindMarkers(samples.combined.list.combined, ident.1 = IDENT1, 
                                                                                                 ident.2 = IDENT2, 
                                                                                                 print.bar = FALSE, only.pos = FALSE)

   ### printing all DIFFERENTIAL MARKERS in a CELL CLUSTER
   x <- as.data.frame(as.matrix(LIST.CLUSTERS.and.DIFFERENTIAL.MARKERS[[i]]))  
   x$gene <- row.names(x)

   write.table(x, file=paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS",  STIM3, STIM5, "CLUSTER", i, "LIST.txt", sep="."), 
                  sep="\t", quote=F, row.names=T, col.names=T)
    
   ######################################################################################################### genes 1-12 or MIN
   ### VlnPlot(object = samples.combined.list.combined, features = gene, split.by = "celltype.stim",  split.plot = FALSE)

   for (j in 1:min(dim(x)[1], 3))           #################################### printing top 12 MARKERS in each CLUSTER or the MIN 
   {
        gene <- x$gene[j] 
        print(gene)   

        ### VLN PLOT
        ### VlnPlot(object = samples.combined.list.combined, features = x$gene[1:6])
        ### VlnPlot(object = samples.combined.list.combined, features = gene, group.by = "orig.ident", split.by = "orig.ident",  split.plot = TRUE)  
        ### VlnPlot(object = samples.combined.list.combined, features = gene, group.by = "orig.ident",  split.plot = FALSE) 
        
        VlnPlot(object = samples.combined.list.combined, features = gene, split.by = "orig.ident",  split.plot = FALSE)
        ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM3, STIM5, "CLUSTER", i, "top.gene", j,  gene, "VIOLIN.PLOT1.png", sep="."), 
               width=120, height=20, units = "cm")
  
        # VlnPlot(object = samples.combined.list.combined, features = gene, split.by = "celltype.stim",  split.plot = FALSE) 
        # ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM3, STIM5, "CLUSTER", i, "top.gene", j,  gene, "VIOLIN.PLOT2.png", sep="."), 
        #       width=60, height=20, units = "cm")
         
        # VlnPlot(object = samples.combined.list.combined, features = gene, group.by = "celltype.stim",  split.plot = FALSE) 
        # ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM3, STIM5, "CLUSTER", i, "top.gene", j,  gene, "VIOLIN.PLOT3.png", sep="."), 
        #       width=60, height=20, units = "cm")

        ### FEATURE PLOT
        ### FeaturePlot(object = samples.combined.list.combined, features = x$gene[1:6] , cols = c("lightgrey", "blue"), reduction = "tsne")   
        # FeaturePlot(object = samples.combined.list.combined, features = gene, cols = c("lightgrey", "blue"), reduction = "tsne", split.by = "orig.ident")        
        # ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS.CLUSTER", i, "top.gene", j, gene, "FEATURE.TSNE.png", sep="."), width=60, height=40, units = "cm")

        ### FEATURE PLOT 
        ### FeaturePlot(object = samples.combined.list.combined, features = x$gene[1:6] , cols = c("lightgrey", "blue"), reduction = "umap")    
        # FeaturePlot(object = samples.combined.list.combined, features = gene, cols = c("lightgrey", "blue"), reduction = "umap", split.by = "orig.ident")        
        # ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS.CLUSTER", i, "top.gene", j, gene, "FEATURE.UMAP.png", sep="."), width=60, height=40, units = "cm")

        FeaturePlot(object = samples.combined.list.combined, features = gene, cols = c("lightgrey", "blue"), reduction = "largeVis", split.by = "orig.ident")        
        ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM3, STIM5, "CLUSTER", i, "top.gene", j, gene, "FEATURE.CONOS.png", sep="."), 
               width=80, height=30, units = "cm")

        ### DOT PLOT
        ### DotPlot(object = samples.combined.list.combined, features = x$gene[1:6])
        # DotPlot(object = samples.combined.list.combined, features = gene,  split.by = "orig.ident")
        # ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM3, STIM5, "CLUSTER", i, "top.gene", j, gene, "DOT.PLOT.png", sep="."), 
        #        width=20, height=20, units = "cm") 

        ### RIDGE PLOT 
        ### RidgePlot(object = samples.combined.list.combined, features = x$gene[1:6])
        RidgePlot(object = samples.combined.list.combined, features = gene, group.by = "orig.ident")
        ggsave(paste(NAME, "figure7.samples.combined.DIFFERENTIAL.MARKERS", STIM3, STIM5, "CLUSTER", i, "top.gene", j, gene, "RIDGE.PLOT.png", sep="."), 
               width=40, height=20, units = "cm") 
     }
  }
}

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################

write.table(as.data.frame(samples.combined.list.combined@meta.data), 
                        file=paste(NAME, "META.DATA", "in.assay.RNA.at.the.end.SCRIPT.txt", sep="."), 
                        sep="\t", quote=FALSE, col.names = TRUE, row.names = TRUE)

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################
