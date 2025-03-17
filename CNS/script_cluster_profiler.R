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

################################################################################
################################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Please provide the input file. Usage: Rscript script.R <input_file>")
}

fin <- args[1]
# fin <- "gs_Cluster_1_Markers.txt"      
fin_name <- basename(fin)
print(fin_name)
w <- read.table(fin, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Filter significant genes (p_val_adj < 0.05 and avg_log2FC > 0.5)
filtered <- w %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5)

dim(filtered)
head(filtered, 3)

# Extract gene names
gene_list <- filtered$gene
head(gene_list, 2)

# Convert gene names to Entrez IDs
gene_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
head(gene_ids, 2)
gene_list_merged <- merge(gene_ids, filtered, by.x = "SYMBOL", by.y = "gene")
gene_list2 <- setNames(gene_list_merged$avg_log2FC, gene_list_merged$ENTREZID)
gene_list2 <- sort(gene_list2, decreasing = TRUE)
head(gene_list2, 2)

################################################################################ GO
################################################################################ ORA
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
  print(dotplot(ego, showCategory = 20))
  
}, error = function(e) {
  cat("An error occurred in GO over-representation analysis:", conditionMessage(e), "\n")
})

print(result)

################################################################################ GO
################################################################################ Enrichment (GSEA)
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
  
  print(dotplot(ego2, showCategory = 20))
  
}, error = function(e) {
  cat("An error occurred in GO enrichment analysis:", conditionMessage(e), "\n")
})

print(result2)

################################################################################ KEGG
################################################################################ ORA
# KEGG Over-Representation Analysis

result <- tryCatch({

  kegg_enrich <- enrichKEGG(gene = gene_ids$ENTREZID,
                            organism = "hsa",  
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05)
    
  if (nrow(kegg_enrich@result) > 0) {
    print(dotplot(kegg_enrich, showCategory = 20))
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

print(result)

################################################################################ KEGG
################################################################################ GSEA
# KEGG Gene Set Enrichment Analysis

result2 <- tryCatch({
    
  kegg_gse <- gseKEGG(geneList   = gene_list2,
                      organism     = 'hsa',
                      minGSSize    = 120,
                      pvalueCutoff = 0.05,
                      verbose      = FALSE)
  
  if (nrow(kegg_gse@result) > 0) {
    print(dotplot(kegg_gse, showCategory = 20))
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

print(result2)

################################################################################ WIKI
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

print(result)

################################################################################ WIKI
################################################################################ GSEA
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

print(result2)

################################################################################ REACTOME
################################################################################ ORA
# Reactome Over-Representation Analysis

result <- tryCatch({

  reactome_gse <- enrichPathway(
    gene = gene_ids$ENTREZID,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  if (nrow(reactome_gse@result) > 0) {
    write.table(reactome_gse@result,
                file = paste0(fin_name, "_Reactome_ORA_Results.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    png(paste0(fin_name, "_Reactome_ORA_Plot.png"), width = 1000, height = 800)
    print(dotplot(reactome_gse, showCategory = 20))
    dev.off()
  } else {
    cat("No enriched Reactome terms found. Please adjust your parameters.\n")
  }
  
}, error = function(e) {
  cat("An error occurred during Reactome over-representation analysis:", conditionMessage(e), "\n")
})

print(result)

################################################################################ REACTOME
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

print(result2)

################################################################################ MSigDB
################################################################################ ORA
# MSigDB Over-Representation Analysis

msig_genesets <- msigdbr(species = "Homo sapiens", category = "C2")
gene_list <- gene_ids$ENTREZID

result <- tryCatch({

  msig_enrich <- enricher(
    gene = gene_list,
    TERM2GENE = msig_genesets[, c("gs_name", "entrez_gene")]
  )
  
  if (!is.null(msig_enrich) && nrow(msig_enrich@result) > 0) {
    
    write.table(msig_enrich@result,
                file = paste0(fin_name, "_MSigDB_ORA_Results.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    png(paste0(fin_name, "_MSigDB_ORA_Plot.png"), width = 1000, height = 800)
    print(dotplot(msig_enrich, showCategory = 20))
    dev.off()
  } else {
    cat("No significant MSigDB pathways found. Try adjusting p-value cutoff or using more genes.\n")
  }
  
}, error = function(e) {
  cat("An error occurred during MSigDB ORA analysis:", conditionMessage(e), "\n")
})

print(result)

################################################################################ MSIGDB
################################################################################ Enrichment (GSEA)
# MSigDB Gene Set Enrichment Analysis

C2_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
head(C2_t2g)

result2 <- tryCatch({

  em2 <- GSEA(gene_list2, TERM2GENE = C2_t2g)
  
  if (!is.null(em2) && nrow(em2@result) > 0) {
    write.table(em2@result,
                file = paste0(fin_name, "_MSigDB_GSEA_Results.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    png(paste0(fin_name, "_MSigDB_GSEA_Plot.png"), width = 1000, height = 800)
    print(dotplot(em2, showCategory = 20))
    dev.off()
  } else {
    cat("No significant MSigDB GSEA pathways found. Try adjusting parameters or using more genes.\n")
  }
  
  em2
}, error = function(e) {
  cat("An error occurred during MSigDB GSEA analysis:", conditionMessage(e), "\n")
  NULL
})

print(result2)

################################################################################
################################################################################