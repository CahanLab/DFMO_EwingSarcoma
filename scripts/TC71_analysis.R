library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(fgsea)
library(singleCellNet)

##### perform Deseq2 ######
expDatGood = singleCellNet::utils_loadObject("../input/TC71_alignment/expDatGood_091319.rda")
sampTabGood = singleCellNet::utils_loadObject("../input/TC71_alignment/sampTabGood_091319.rda")
counts_mat<-expDatGood

# Convert to a matrix of ints, required by DESeq
int_counts_mat <- apply (counts_mat, c (1, 2), function (x) {
  (as.integer(x))
})

#run DESeq2
#deseq_wrapper(int_counts_mat, newSampTab=sampTabGood, group1 = "lung", group2 = "bone_marrow", study_name = "os", mydate, design_col="description1", design=~description1, alpha=0.05)
sampTabGood$description1 = as.factor(sampTabGood$description1)

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = int_counts_mat,
  colData = sampTabGood,
  design = ~description1)  

dds <- DESeq(ddsFullCountTable)

group1 <- "lung"
group2 <- "bone_marrow"
res <- results(dds, contrast=c("description1", group1, group2))

dir.create("../output/TC71_analysis/DE_genes", recursive = TRUE)
saveRDS(res, file = '../output/TC71_analysis/DE_genes/res.rds')

##### run fgsea #####

fgsea_wrapper <- function(res, rows = 20, rank = "stat", path = "./") {
  
  # Detailed tutorial here: https://bioconductor.org/packages/3.7/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
  # Example usage: fgsea_wrapper(int_counts_mat, newSampTab, "PAPD5 CRISPR KO", "WT")
  
  # Load gene sets downloaded from MSigDB
  # C2: Curated gene sets
  c2.all <- gmtPathways("../input/GeneSets/c2.all.v7.0.symbols.gmt")
  # C5: GO gene sets 
  c5.all <- gmtPathways("../input/GeneSets/c5.all.v7.0.symbols.gmt")
  # C6
  c6.all <- gmtPathways("../input/GeneSets/c6.all.v7.0.symbols.gmt")
  # C3 MIR only
  c3.mir <- gmtPathways("../input/GeneSets/c3.mir.v7.0.symbols.gmt")
  
  deseq_stat <- res[,rank]
  names(deseq_stat) <- rownames(res)
  deseq_stat <- deseq_stat[!is.na(deseq_stat)]
  
  #plot 20 up and 20 downreg genes
  #rows <- 20
  
  ## look at c2
  c2.fgseaRes <- as.data.frame(fgsea(pathways = c2.all, 
                                     stats = deseq_stat,
                                     minSize=3, maxSize=500, # min and max size limits on the size of gene sets considered
                                     nperm=1000))
  up <- c2.fgseaRes[c2.fgseaRes["ES"] > 0, ]
  down <- c2.fgseaRes[c2.fgseaRes["ES"] < 0, ]
  topPathwaysUp <- up[head(order(up["padj"]), n=rows), "pathway"]
  topPathwaysDown <- down[head(order(down["padj"]), n=rows), "pathway"]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  topPathways <- topPathways[!is.na(topPathways)]
  pdf(file=paste0(path, "c2_", group1, "_vs_", group2, "_fGSEA_Results1.pdf"), width=16, height=12)
  plotGseaTable(c2.all[topPathways], deseq_stat, c2.fgseaRes, 
                gseaParam = 0.5)
  dev.off()
  
  ## look at c3.mir
  c3.mir.fgseaRes <- as.data.frame(fgsea(pathways = c3.mir, 
                                         stats = deseq_stat,
                                         minSize=3, maxSize=500, # min and max size limits on the size of gene sets considered
                                         nperm=1000))
  up <- c3.mir.fgseaRes[c3.mir.fgseaRes["ES"] > 0, ]
  down <- c3.mir.fgseaRes[c3.mir.fgseaRes["ES"] < 0, ]
  topPathwaysUp <- up[head(order(up["padj"]), n=rows), "pathway"]
  topPathwaysDown <- down[head(order(down["padj"]), n=rows), "pathway"]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  topPathways <- topPathways[!is.na(topPathways)]
  pdf(file=paste0(path, "c3_mir_", group1, "_vs_", group2, "_fGSEA_Results.pdf"), width=16, height=12)
  plotGseaTable(c3.mir[topPathways], deseq_stat, c3.mir.fgseaRes, 
                gseaParam = 0.5)
  dev.off()
  
  ## look at c5
  c5.fgseaRes <- as.data.frame(fgsea(pathways = c5.all, 
                                     stats = deseq_stat,
                                     minSize=3, maxSize=500, # min and max size limits on the size of gene sets considered
                                     nperm=1000))
  up <- c5.fgseaRes[c5.fgseaRes["ES"] > 0, ]
  down <- c5.fgseaRes[c5.fgseaRes["ES"] < 0, ]
  topPathwaysUp <- up[head(order(up["padj"]), n=rows), "pathway"]
  topPathwaysDown <- down[head(order(down["padj"]), n=rows), "pathway"]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  topPathways <- topPathways[!is.na(topPathways)]
  pdf(file=paste0(path, "c5_", group1, "_vs_", group2, "_fGSEA_Results.pdf"), width=16, height=12)
  plotGseaTable(c5.all[topPathways], deseq_stat, c5.fgseaRes, 
                gseaParam = 0.5)
  dev.off()
  
  c6.fgseaRes <- as.data.frame(fgsea(pathways = c6.all, 
                                     stats = deseq_stat,
                                     minSize=3, maxSize=500, # min and max size limits on the size of gene sets considered
                                     nperm=1000)) 
  up <- c6.fgseaRes[c6.fgseaRes["ES"] > 0, ]
  down <- c6.fgseaRes[c6.fgseaRes["ES"] < 0, ]
  topPathwaysUp <- up[head(order(up["padj"]), n=rows), "pathway"]
  topPathwaysDown <- down[head(order(down["padj"]), n=rows), "pathway"]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  topPathways <- topPathways[!is.na(topPathways)]
  pdf(file=paste0(path, "c6_", group1, "_vs_", group2, "_fGSEA_Results.pdf"), width=16, height=12)
  plotGseaTable(c6.all[topPathways], deseq_stat, c6.fgseaRes, 
                gseaParam = 0.5)
  dev.off()
  
  c2.fgseaRes$leadingEdge=vapply(c2.fgseaRes$leadingEdge, paste, collapse = ", ", character(1L))
  write.csv(c2.fgseaRes, file = paste0(path, "c2.fgseaRes.csv"))
  
  c5.fgseaRes$leadingEdge=vapply(c5.fgseaRes$leadingEdge, paste, collapse = ", ", character(1L))
  write.csv(c5.fgseaRes, file = paste0(path, "c5.fgseaRes.csv"))
  
  c6.fgseaRes$leadingEdge=vapply(c6.fgseaRes$leadingEdge, paste, collapse = ", ", character(1L))
  write.csv(c6.fgseaRes, file = paste0(path, "c6.fgseaRes.csv"))
  
  c3.mir.fgseaRes$leadingEdge=vapply(c3.mir.fgseaRes$leadingEdge, paste, collapse = ", ", character(1L))
  write.csv(c3.mir.fgseaRes, file = paste0(path, "c3.mir.fgseaRes.csv"))
}

#GSEA analysis
dir.create("../output/TC71_analysis/gsea", recursive = TRUE)
fgsea_wrapper(res, path = '../output/TC71_analysis/gsea/')

