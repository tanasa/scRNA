################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

library("Seurat")
library("cowplot")
library("dplyr")

library("gmodels")
library("gplots")

library("edgeR")
library("gridExtra")

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

# a name for the sample :

NAME="RETINA.filter900genes"
FILE="GSE63472_P14Retina_merged_digital_expression.txt"

################################################################################################################################################################
################################################################################################################################################################ READING

# we are starting from a DATA TABLE, a MATRIX that has GENES * CELLS 

pbmc.data <- read.delim(FILE, header=T, sep="\t", stringsAsFactors=F)

row.names(pbmc.data) <- pbmc.data$gene
pbmc.data <- pbmc.data[, -1]
head(pbmc.data[1:5,1:5])

################################################################################################################################################################
################################################################################################################################################################ the OBJECT
################################################################################################################################################################ min genes 900

pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 10, min.genes = 900, project = NAME)

# str(pbmc)
# pbmc@raw.data
# pbmc@data
# pbmc@meta.data

# head(pbmc@raw.data)
# head(pbmc@meta.data)

# IN ORDER tO ACCESS THE MATRIX and to write it :

# pbmc_data <- as.data.frame(as.matrix(pbmc@raw.data))

# write.table(pbmc_data, 
#            file = paste(NAME, ".SAMPLE.printing.the.matrix.as.DF.txt", sep=""), 
#            sep="\t", 
#            quote=FALSE, 
#            row.names = TRUE, col.names = TRUE) 

################################################################################################################################################################
################################################################################################################################################################ MITO GENES

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ]) / Matrix::colSums(pbmc@raw.data)

pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")

VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
ggsave(paste(NAME, "..figure1.ViolinPlot.nUMI.nGene.png", sep=""), width=20, height=10, units="cm")

png(paste(NAME, "..figure1.GenePlot.nUMI.nGene.png", sep=""), width=1000, height=500) 
par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
dev.off()

################################################################################################################################################################
################################################################################################################################################################ FILTERING

pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
                                   low.thresholds = c(900, -Inf), 
                                   high.thresholds = c(Inf, 0.2))

################################################################################################################################################################
################################################################################################################################################################ NORMALIZATION

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                                     scale.factor = 10000)

# writing the NORMALIZED DATA : 

# pbmc_data_normalized <- as.data.frame(as.matrix(pbmc@data))

# write.table(pbmc_data_normalized, 
#             file = paste(NAME, ".SAMPLE.printing.the.matrix.as.DF.normalized.txt", sep=""),  
#            sep="\t", 
#            quote=FALSE, 
#            row.names = TRUE, col.names = TRUE) 

################################################################################################################################################################
################################################################################################################################################################ VARIABLE GENES

png(paste(NAME,"..figure2.VARIABLE.genes.png", sep=""), width=1000, height=1000)

pbmc <- FindVariableGenes(object = pbmc, 
                          mean.function = ExpMean, 
                          dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, 
                          x.high.cutoff = 3, 
                          y.cutoff = 0.5)

dev.off()

# Printing the VARIABLE GENES : 

number_variable_genes <- length(x = pbmc@var.genes)

write.table(number_variable_genes, 
            file=paste(NAME, "..figure2.VARIABLE.genes.count.txt", sep=""), 
            quote=FALSE, 
            sep="\t", row.names=FALSE, col.names=FALSE)

################################################################################################################################################################
################################################################################################################################################################ SCALING

pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

################################################################################################################################################################
################################################################################################################################################################ PCA

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, 
                              do.print = TRUE, 
                              pcs.print = 1:10, 
                              genes.print = 10)

################################################################################################################################################################
################################################################################################################################################################

png(paste(NAME, "..figure3.viz.PCA.png", sep=""), width=1000, height=1000)
VizPCA(object = pbmc, pcs.use = 1:12)
dev.off()

# Here displaying the PCA analysis results : 

png(paste(NAME, "..figure3.viz.PCA.here.PC1.PC2.png", sep=""), width=500, height=500)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
dev.off()

png(paste(NAME, "..figure3.viz.PCA.here.PC1.PC3.png", sep=""), width=500, height=500)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 3)
dev.off()

png(paste(NAME, "..figure3.viz.PCA.here.PC2.PC3.png", sep=""), width=500, height=500)
PCAPlot(object = pbmc, dim.1 = 2, dim.2 = 3)
dev.off()

################################################################################################################################################################
################################################################################################################################################################

pbmc <- ProjectPCA(object = pbmc, do.print = TRUE)

################################################################################################################################################################
################################################################################################################################################################ PCA heatmaps

png(paste(NAME, "..figure4.viz.HEATMAP.here.PC.1.cells.500.genes.50.png", sep=""), width=500, height=500)
PCHeatmap(object = pbmc, pc.use = 1, cells.use = 500, num.genes=50, do.balanced = TRUE, label.columns = FALSE)
dev.off()

png(paste(NAME, "..figure4.viz.HEATMAP.here.PC.2.cells.500.genes.50.png", sep=""), width=500, height=500)
PCHeatmap(object = pbmc, pc.use = 2, cells.use = 500, num.genes=50, do.balanced = TRUE, label.columns = FALSE)
dev.off()

png(paste(NAME, "..figure4.viz.HEATMAP.here.PC.3.cells.500.genes.50.png", sep=""), width=500, height=500)
PCHeatmap(object = pbmc, pc.use = 3, cells.use = 500, num.genes=50, do.balanced = TRUE, label.columns = FALSE)
dev.off()

png(paste(".SAMPLE.figure4.viz.HEATMAP.here.PC.1.to.PC.12.cells.500.genes.50.png", sep=""),  width=500, height=500)
PCHeatmap(object = pbmc, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()

################################################################################################################################################################
################################################################################################################################################################ PCA SIGNIFICANT COMPONENTS

pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = FALSE, num.cores=8, do.par=TRUE)

################################################################################################################################################################
################################################################################################################################################################ 

png(paste(NAME, "..figure5.viz.PCA.JackStraw.here.PC.1.to.PC.12.with.JackStrawPlot.png", sep=""), width=500, height=500)
JackStrawPlot(object = pbmc, PCs = 1:20)
dev.off()

png(paste(NAME, "..figure5.viz.PCA.PCelbow.here.PC.1.to.PC.12.with.PCElbowplot.png", sep=""), width=500, height=500)
PCElbowPlot(object = pbmc)
dev.off()

################################################################################################################################################################
################################################################################################################################################################ CLUSTERING :
################################################################################################################################################################ RESOLUTION 0.6
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

RESOLUTION = 0.6

# CLUSTERING the cells and PRINT CLUSTER PARAMETERS : 

pbmc <- FindClusters(object = pbmc, 
                     reduction.type = "pca", 
                     dims.use = 1:12, 
                     resolution = RESOLUTION, 
                     print.output = 0, save.SNN = TRUE, plot.SNN=TRUE)

PrintFindClustersParams(object = pbmc)

################################################################################################################################################################
################################################################################################################################################################ TSNE :

pbmc <- RunTSNE(object = pbmc, dims.use = 1:12)

TSNEPlot(object = pbmc, do.label=TRUE)
ggsave(paste(NAME, "..figure5.viz.TSNE.for.RESOLUTION.", RESOLUTION, ".png", sep=""))

################################################################################################################################################################ to add 
################################################################################################################################################################ CELLS COUNTS
################################################################################################################################################################ on the display
### about clusters : https://github.com/satijalab/seurat/issues/354
### if you'd like a count of the number of cells in each cluster, use the [`table`] function to get counts of each cluster ID stored in the [`ident` slot].

# Calculate number of cells per cluster from pbmc@ident
cell.num <- table(pbmc@ident)

# Add cell number per cluster to cluster labels
ClusterLabels = paste("Cluster", names(cell.num), paste0("(n = ", cell.num, ")"))

# Order legend labels in plot in the same order as 'ClusterLabels'
ClusterBreaks = names(cell.num)

# Plot tSNE with new legend labels for clusters
TSNEPlot(object = pbmc, do.return = T, 
                        do.label = T, 
                        no.legend = FALSE, 
                        pt.size = 0.6) +
                        scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
                        labs(x = "t-SNE 1", y = "t-SNE 2") 

ggsave(paste(NAME, "..figure5.viz.TSNE.for.RESOLUTION", RESOLUTION, "and.the.number.of.cells", "png", sep="."))

write.table(data.frame(ClusterLabels), 
            file=paste(NAME, "figure5.viz.TSNE.for.RESOLUTION", RESOLUTION, "and.the.list.of.CLUSTERS", "txt", sep="."), 
            sep="\t", quote=F, row.names=T, col.names=T)

################################################################################################################################################################
################################################################################################################################################################ MARKERS :
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)                   #### we may change only.POS = FALSE 

pbmc.markers.df <- as.data.frame(pbmc.markers)

dim(pbmc.markers.df)

write.table(pbmc.markers.df, 
            file=paste(NAME, "..figure6.viz.TSNE.for.RESOLUTION", RESOLUTION, "MARKERS.all.txt", sep="."),
            sep="\t",
            quote=FALSE, row.names=FALSE)

pbmc.markers.df.clusters <- pbmc.markers %>% group_by(cluster)      

list.markers <- split(pbmc.markers.df.clusters, pbmc.markers.df.clusters$cluster)

################################################################################################################################################################ making figures 
################################################################################################################################################################ for each MARKER :

for (i in 1 : length(list.markers) )                                   
{
    print(i) 

    for (j in 1:dim(list.markers[[i]]))                                                  
    {
        gene <- list.markers[[i]]$gene[j] 
        print(gene)  

        k = i-1 ### in order to MATCH the NUMBER on TSNE PLOT

        VlnPlot(object = pbmc, features.plot = gene, x.lab.rot = TRUE)
        ggsave(paste(NAME, "figure7", "RESOLUTION", RESOLUTION, "CLUSTER", k, "gene", j, "plot.VIOLIN.for", gene, "png", sep="."))
 
        FeaturePlot(object = pbmc, features.plot = gene , cols.use = c("grey", "blue"), reduction.use = "tsne")        
        ggsave(paste(NAME, "figure7", "RESOLUTION", RESOLUTION, "CLUSTER", k, "gene", j, "plot.FEATURE.for", gene, "png", sep="."))
       
     }
}

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

# at the end, shall we save the session :
# saveRDS(pbmc, file = paste("SAMPLE, ".analysis.with.SEURAT.rds"))

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################ the MARKERS
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################

#####################################################################################################  resume the display : markers RGC
#####################################################################################################  

SET_GENES <- c("THY1", "RBPMS", "POU4F1", "POU4F2")
gene_name <- "THY1.RBPMS.POU4F1.POU4F2"
cell_name <- "RGC"

for (i in SET_GENES){
    if (i %in% row.names(pbmc@data)) {
        # png(paste(NAME, "..figure8.markers", "cell", cell_name, "gene", i, "png", sep="."), width=400, height=400)
        p <- FeaturePlot(object = pbmc, features.plot = i, cols.use = c("lightgrey", "blue"), reduction.use = "tsne",  no.legend = FALSE, do.return=TRUE)
        ggsave(paste(NAME, ".figure9.markers", "cell", cell_name, "gene", i, "png", sep="."))
        # dev.off()
    }
}

#####################################################################################################  resume the display : marker RODS : RHO 
#####################################################################################################

SET_GENES <- c("RHO")
gene_name <- "RHO"
cell_name <- "RODS"

for (i in SET_GENES){
    if (i %in% row.names(pbmc@data)) {
        # png(paste(NAME, "..figure8.markers", "cell", cell_name, "gene", i, "png", sep="."), width=400, height=400)
        FeaturePlot(object = pbmc, features.plot = i, cols.use = c("lightgrey", "blue"), reduction.use = "tsne",  no.legend = FALSE, do.return=TRUE)
        ggsave(paste(NAME, ".figure9.markers", "cell", cell_name, "gene", i, "png", sep="."))
        # dev.off()
    }
}

#####################################################################################################  resume the display : marker CONES : ARR3
#####################################################################################################  

SET_GENES <- c("ARR3")
gene_name <- "ARR3"
cell_name <- "CONES"

for (i in SET_GENES){
    if (i %in% row.names(pbmc@data)) {
        # png(paste(NAME, "..figure8.markers", "cell", cell_name, "gene", i, "png", sep="."), width=400, height=400)
        FeaturePlot(object = pbmc, features.plot = i, cols.use = c("lightgrey", "blue"), reduction.use = "tsne",  no.legend = FALSE, do.return=TRUE)
        ggsave(paste(NAME, ".figure9.markers", "cell", cell_name, "gene", i, "png", sep="."))
        # dev.off()
    }
}

#####################################################################################################  resume the display : marker BP : VSX2 
#####################################################################################################  resume the display : marker BP : PRKCA 
#####################################################################################################  resume the display : marker BP : GRM6 

SET_GENES <- c("VSX2","PRKCA","GRM6")
gene_name <- "VSX2.PRKCA.GRM6"
cell_name <- "BIPOLAR_CELLS"

for (i in SET_GENES){
    if (i %in% row.names(pbmc@data)) {
       # png(paste(NAME, "..figure8.markers", "cell", cell_name, "gene", i, "png", sep="."), width=400, height=400)
       FeaturePlot(object = pbmc, features.plot = i, cols.use = c("lightgrey", "blue"), reduction.use = "tsne")
       ggsave(paste(NAME, ".figure9.markers", "cell", cell_name, "gene", i, "png", sep="."))
       # dev.off()
    }
}

#####################################################################################################  resume the display : marker MG: SOX9 
#####################################################################################################  resume the display : marker MG: GLUL
#####################################################################################################  resume the display : marker MG: CLU
#####################################################################################################  resume the display : marker MG: DKK3
#####################################################################################################  resume the display : marker MG: S100A16
### from the article : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2665263/ 
### Glutamine synthetase (Glul), 
### clusterin (Clu), 
### dickkopf homolog 3 (Dkk3) (Blackshaw et al., 2004), 
### and S100 calcium binding protein A16 (Seigel et al., 1996) were successfully detected  

SET_GENES <- c("SOX9","GLUL","CLU", "DKK3", "S100A16")
gene_name <- "SOX9.GLUL.CLU.DKK3.S100A16"
cell_name <- "MULLER_GLIA"

for (i in SET_GENES){
    if (i %in% row.names(pbmc@data)) {
        # png(paste(NAME, "..figure8.markers", "cell", cell_name, "gene", i, "png", sep="."), width=400, height=400)
        FeaturePlot(object = pbmc, features.plot = i, cols.use = c("lightgrey", "blue"), reduction.use = "tsne")
        ggsave(paste(NAME, ".figure9.markers", "cell", cell_name, "gene", i, "png", sep="."))
        # dev.off()
     }
}

#####################################################################################################  resume the display : marker AMACRINE : PAX6 
#####################################################################################################  resume the display : marker AMACRINE : GAD1
#####################################################################################################  resume the display : marker AMACRINE : GLT1
#####################################################################################################  resume the display : marker AMACRINE : TFAP2A

SET_GENES <- c("PAX6","GAD1","TFAP2A")
gene_name <- "PAX6.GAD1.TFAP2A"
cell_name <- "AMACRINE"

for (i in SET_GENES){
    if (i %in% row.names(pbmc@data)) {
       # png(paste(NAME, "..figure8.markers", "cell", cell_name, "gene", i, "png", sep="."), width=400, height=400)
       FeaturePlot(object = pbmc, features.plot = i, cols.use = c("lightgrey", "blue"), reduction.use = "tsne")
       ggsave(paste(NAME, ".figure9.markers", "cell", cell_name, "gene", i, "png", sep="."))
       # dev.off()
    }
}

#####################################################################################################  resume the display : marker ENDOTHELIAL CELLS : PECAM1
#####################################################################################################  

SET_GENES <- c("PECAM1")
gene_name <- "PECAM1"
cell_name <- "ENDOTHELIAL_CELLS"

for (i in SET_GENES){
    if (i %in% row.names(pbmc@data)) {
       # png(paste(NAME, "..figure8.markers", "cell", cell_name, "gene", i, "png", sep="."), width=400, height=400)
       FeaturePlot(object = pbmc, features.plot = i, cols.use = c("lightgrey", "blue"), reduction.use = "tsne")
       ggsave(paste(NAME, ".figure9.markers", "cell", cell_name, "gene", i, "png", sep="."))
       # dev.off()
    }
}
 
#####################################################################################################  resume the display : marker ASTROCYTES : GFAP
#####################################################################################################  

SET_GENES <- c("GFAP")
gene_name <- "GFAP"
cell_name <- "ASTROCYTES"

for (i in SET_GENES){
    if (i %in% row.names(pbmc@data)) {
       # png(paste(NAME, "..figure8.markers", "cell", cell_name, "gene", i, "png", sep="."), width=400, height=400)
       FeaturePlot(object = pbmc, features.plot = i, cols.use = c("lightgrey", "blue"), reduction.use = "tsne")
       ggsave(paste(NAME, ".figure9.markers", "cell", cell_name, "gene", i, "png", sep="."))
       # dev.off()
    }
}

#####################################################################################################  resume the display : marker PERYCYTES : NG2
#####################################################################################################  

SET_GENES <- c("CSPG4")
gene_name <- "CSPG4"
cell_name <- "PERICYTES"

for (i in SET_GENES){
    if (i %in% row.names(pbmc@data)) {
       # png(paste(NAME, "..figure8.markers", "cell", cell_name, "gene", i, "png", sep="."), width=400, height=400)
       FeaturePlot(object = pbmc, features.plot = i, cols.use = c("lightgrey", "blue"), reduction.use = "tsne")
       ggsave(paste(NAME, ".figure9.markers", "cell", cell_name, "gene", i, "png", sep="."))
       # dev.off()
    }
}

#####################################################################################################  resume the display : marker of CHAT
#####################################################################################################  

SET_GENES <- c("CHAT")
gene_name <- "CHAT"
cell_name <- "COLINERGIC"

for (i in SET_GENES){
    if (i %in% row.names(pbmc@data)) {
       # png(paste(NAME, "..figure8.markers", "cell", cell_name, "gene", i, "png", sep="."), width=400, height=400)
       FeaturePlot(object = pbmc, features.plot = i, cols.use = c("lightgrey", "blue"), reduction.use = "tsne")
       ggsave(paste(NAME, ".figure9.markers", "cell", cell_name, "gene", i, "png", sep="."))
       # dev.off()
    }
}

#####################################################################################################  resume the display : marker of MICROGLIA
#####################################################################################################  

SET_GENES <- c("CX3CR1")
gene_name <- "CX3CR1"
cell_name <- "MICROGLIA"

for (i in SET_GENES){
    if (i %in% row.names(pbmc@data)) {
       # png(paste(NAME, "..figure8.markers", "cell", cell_name, "gene", i, "png", sep="."), width=400, height=400)
       FeaturePlot(object = pbmc, features.plot = i, cols.use = c("lightgrey", "blue"), reduction.use = "tsne")
       ggsave(paste(NAME, ".figure9.markers", "cell", cell_name, "gene", i, "png", sep="."))
       # dev.off()
    }
}

################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################ 
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
