library(ggplot2)
library(DESeq2)
library(stringr)
library(biomaRt)
library(enrichR)

# get all the samples from DFMO v control quantification 
all_samples = list.dirs("../input/DFMO_control_alignment/quants/", recursive = FALSE, full.names = FALSE)

###### compile all the raw counts into matrix ######
compiled_df = data.frame()
for(samp_name in all_samples) {
  temp_read = read.csv(file.path("../input/DFMO_control_alignment/quants/", samp_name, 'quant.sf'), sep = '\t')
  temp_df = data.frame(row.names = temp_read$Name, 
                       'counts' = round(temp_read$NumReads))
  colnames(temp_df) = stringr::str_remove(samp_name, "_L001_quant")
  if(nrow(compiled_df) == 0) { 
    compiled_df = temp_df
  } else {
    compiled_df = merge(compiled_df, temp_df, by = 'row.names')
    rownames(compiled_df) = compiled_df$Row.names
    compiled_df$Row.names = NULL
  }
}

# construct meta table 
meta_tab = data.frame(row.names = colnames(compiled_df), 
                      'sequence_id' = colnames(compiled_df))
meta_tab$PDX_id = str_split_fixed(rownames(meta_tab), "_", n = 2)[, 1]

# filling in the below based on email patrick sent 
meta_tab$treated = 'control'
meta_tab[grep("E", meta_tab$PDX_id), 'treated'] = 'DFMO'
meta_tab$gender = 'male'
meta_tab[meta_tab$PDX_id %in% c("E1", "E3", 'E4', 'C3'), 'gender'] = 'female'
meta_tab$days = 0
meta_tab[meta_tab$PDX_id %in% c('E1', 'E3'), 'days'] = 47
meta_tab[meta_tab$PDX_id %in% c('E4'), 'days'] = 40
meta_tab[meta_tab$PDX_id %in% c('E6'), 'days'] = 64

meta_tab$type = NA
meta_tab[grep("C", meta_tab$PDX_id), 'type'] = 'control'
meta_tab[grep("E", meta_tab$PDX_id), 'type'] = 'treated'
meta_tab[grep("L", meta_tab$PDX_id), 'type'] = 'late-treated'

# start DE analysis 
# start with the human transcripts 
mouse_counts = compiled_df[grep("ENSMUST", rownames(compiled_df)), ]
human_counts = compiled_df[grep("ENST", rownames(compiled_df)), ]

dir.create("../output/DFMO_control_analyis/human_mouse_counts", recursive = TRUE)
saveRDS(mouse_counts, file = '../output/DFMO_control_analyis/human_mouse_counts/mouse_counts.rds')
saveRDS(human_counts, file = '../output/DFMO_control_analyis/human_mouse_counts/human_counts.rds')
write.csv(meta_tab, file = '../output/DFMO_control_analyis/human_mouse_counts/meta_tab.csv')

stats_tab = data.frame('human_reads' = format(apply(human_counts, MARGIN = 2, FUN = sum), scientific = TRUE), 
                       'mouse_reads' = format(apply(mouse_counts, MARGIN = 2, FUN = sum), scientific = TRUE))
write.csv(stats_tab, file = '../output/DFMO_control_analyis/human_mouse_counts/reads_num_counts.csv')

##### find ways to convert the ensemble transcripts into genes ######
# load in the ensemble transcript to gene converting table 
header_tab = read.csv("../input/DFMO_control_alignment/reference_genomes/hg38_mm10_header.txt", sep = ' ', header = FALSE)
header_tab = header_tab[, c(1, 7)]
colnames(header_tab) = c("ensembl_id", 'gene_id')
header_tab$ensembl_id = str_remove_all(header_tab$ensembl_id, ">")
header_tab$gene_id = str_remove_all(header_tab$gene_id, 'gene_symbol:')
header_tab = header_tab[header_tab$ensembl_id %in% rownames(compiled_df), ]
rownames(header_tab) = header_tab$ensembl_id

##### perform DEseq2 #####
# calculate using just human counts 
dds <- DESeqDataSetFromMatrix(countData=human_counts, 
                              colData=meta_tab, 
                              design=~treated, tidy = FALSE)
keep = rowSums(counts(dds) >= 100) >= 3
dds = dds[keep, ]

dds <- DESeq(dds)

###### perform analysis of DEseq2 results #####
dir.create("../output/DFMO_control_analyis/human_DE_analysis")

withr::with_dir('../output/DFMO_control_analyis/human_DE_analysis', {
  norm_counts = counts(dds, normalize = TRUE)
  rownames(norm_counts) = header_tab[rownames(norm_counts), 'gene_id']
  write.csv(norm_counts, file = 'human_norm_exp.csv')
  saveRDS(dds, file = 'human_DESeq2_object.rds')
})

res <- results(dds)
res_df = as.data.frame(res)
res_df = res_df[complete.cases(res_df), ]

DFMO_up = res_df[res_df$log2FoldChange > 0, ]
DFMO_up = DFMO_up[DFMO_up$padj < 0.05, ]
DFMO_up$gene = header_tab[rownames(DFMO_up), 'gene_id']

enrichR::setEnrichrSite('Enrichr')
databases = c('GO_Biological_Process_2021', 'GO_Molecular_Function_2021','KEGG_2021_Human', 'WikiPathway_2021_Human', 'Reactome_2016', 'BioCarta_2016', 'MSigDB_Hallmark_2020')
enrichment_analysis = enrichR::enrichr(DFMO_up$gene, databases)
dir.create("../output/DFMO_control_analyis/human_DE_analysis/DFMO_genes")

withr::with_dir('../output/DFMO_control_analyis/human_DE_analysis/DFMO_genes', {
  write.csv(DFMO_up, file = 'human_DFMO_up_genes.csv')
  saveRDS(enrichment_analysis, file = 'DFMO_enrichment_analysis.rds')
  
  for(analysis_name in names(enrichment_analysis)) { 
    temp_df = enrichment_analysis[[analysis_name]]
    temp_df = temp_df[temp_df$Adjusted.P.value < 0.05, ]
    write.csv(temp_df, paste0(analysis_name, "_enrichment_analysis.csv"))
    
    plot_df = data.frame("category" = temp_df$Term, 
                         'log_adj_p' = -log(temp_df$Adjusted.P.value))
    plot_df = plot_df[order(plot_df$log_adj_p, decreasing = TRUE), ]
    
    p<-ggplot(data=plot_df, aes(x=reorder(category, log_adj_p), y=log_adj_p)) +
      geom_bar(stat="identity") + 
      ylab("-log adj-p") +
      xlab("Terms")+
      ggtitle(analysis_name) + 
      theme_bw() +
      coord_flip()
    ggsave(filename = paste0("significant_terms_", analysis_name, '.png'), width = 12, height = 10)
  }
})

dir.create("../output/DFMO_control_analyis/human_DE_analysis/control_genes")

DFMO_down = res_df[res_df$log2FoldChange < 0, ]
DFMO_down = DFMO_down[DFMO_down$padj < 0.05, ]
DFMO_down$gene = header_tab[rownames(DFMO_down), 'gene_id']

enrichment_analysis = enrichR::enrichr(DFMO_down$gene, databases)

withr::with_dir('../output/DFMO_control_analyis/human_DE_analysis/control_genes', {
  write.csv(DFMO_down, file = 'human_control_up_genes.csv')
  saveRDS(enrichment_analysis, file = 'control_enrichment_analysis.rds')
  
  for(analysis_name in names(enrichment_analysis)) { 
    temp_df = enrichment_analysis[[analysis_name]]
    temp_df = temp_df[temp_df$Adjusted.P.value < 0.05, ]
    write.csv(temp_df, paste0(analysis_name, "_enrichment_analysis.csv"))
    
    plot_df = data.frame("category" = temp_df$Term, 
                         'log_adj_p' = -log(temp_df$Adjusted.P.value))
    plot_df = plot_df[order(plot_df$log_adj_p, decreasing = TRUE), ]
    
    p<-ggplot(data=plot_df, aes(x=reorder(category, log_adj_p), y=log_adj_p)) +
      geom_bar(stat="identity") + 
      ylab("-log adj-p") +
      xlab("Terms")+
      ggtitle(analysis_name) + 
      theme_bw() +
      coord_flip()
    ggsave(filename = paste0("significant_terms_", analysis_name, '.png'), width = 12, height = 10)
  }
})

###### this is for the mouse counts data (NOT USED IN THE MANUSCRIPT) ######
dds <- DESeqDataSetFromMatrix(countData=mouse_counts, 
                              colData=meta_tab, 
                              design=~treated, tidy = FALSE)
keep = rowSums(counts(dds) >= 100) >= 3
dds = dds[keep, ]

dds <- DESeq(dds)

dir.create("../output/DFMO_control_analyis/mouse_DE_analysis")

withr::with_dir('../output/DFMO_control_analyis/mouse_DE_analysis', {
  norm_counts = counts(dds, normalize = TRUE)
  rownames(norm_counts) = header_tab[rownames(norm_counts), 'gene_id']
  write.csv(norm_counts, file = 'mouse_norm_exp.csv')
  saveRDS(dds, file = 'mouse_DESeq2_object.rds')
})

dds = readRDS("../output/DFMO_control_analyis/mouse_DE_analysis/mouse_DESeq2_object.rds")

res <- results(dds)
res_df = as.data.frame(res)
res_df = res_df[complete.cases(res_df), ]

DFMO_up = res_df[res_df$log2FoldChange > 0, ]
DFMO_up = DFMO_up[DFMO_up$padj < 0.05, ]
DFMO_up$gene = header_tab[rownames(DFMO_up), 'gene_id']

databases = c('KEGG_2019_Mouse', 'WikiPathways_2019_Mouse', 'Mouse_Gene_Atlas', 'GO_Biological_Process_2021')
enrichment_analysis = enrichR::enrichr(DFMO_up$gene, databases)
enrichment_df = enrichment_analysis$Mouse_Gene_Atlas

dir.create("../output/DFMO_control_analyis/mouse_DE_analysis/DFMO_genes")
withr::with_dir('../output/DFMO_control_analyis/mouse_DE_analysis/DFMO_genes', {
  write.csv(DFMO_up, file = 'mouse_DFMO_up_genes.csv')
  saveRDS(enrichment_analysis, file = 'DFMO_enrichment_analysis.rds')
  
  for(analysis_name in names(enrichment_analysis)) { 
    temp_df = enrichment_analysis[[analysis_name]]
    temp_df = temp_df[temp_df$Adjusted.P.value < 0.05, ]
    write.csv(temp_df, paste0(analysis_name, "_enrichment_analysis.csv"))
    
    plot_df = data.frame("category" = temp_df$Term, 
                         'log_adj_p' = -log(temp_df$Adjusted.P.value))
    plot_df = plot_df[order(plot_df$log_adj_p, decreasing = TRUE), ]
    plot_df = plot_df[1:20, ]
    plot_df = plot_df[complete.cases(plot_df), ]
    
    p<-ggplot(data=plot_df, aes(x=reorder(category, log_adj_p), y=log_adj_p)) +
      geom_bar(stat="identity") + 
      ylab("-log adj-p") +
      xlab("Terms")+
      ggtitle(analysis_name) + 
      theme_bw() +
      coord_flip()
    ggsave(filename = paste0("significant_terms_", analysis_name, '.png'), width = 12, height = 14)
  }
})

# this is for genesets that are down regulated by DFMO 
DFMO_down = res_df[res_df$log2FoldChange < 0, ]
DFMO_down = DFMO_down[DFMO_down$padj < 0.05, ]
DFMO_down$gene = header_tab[rownames(DFMO_down), 'gene_id']

enrichment_analysis = enrichR::enrichr(DFMO_down$gene, databases)

dir.create("../output/DFMO_control_analyis/mouse_DE_analysis/control_genes")
withr::with_dir('../output/DFMO_control_analyis/mouse_DE_analysis/control_genes', {
  write.csv(DFMO_down, file = 'mouse_control_up_genes.csv')
  saveRDS(enrichment_analysis, file = 'control_enrichment_analysis.rds')
  
  for(analysis_name in names(enrichment_analysis)) { 
    temp_df = enrichment_analysis[[analysis_name]]
    temp_df = temp_df[temp_df$Adjusted.P.value < 0.05, ]
    write.csv(temp_df, paste0(analysis_name, "_enrichment_analysis.csv"))
    
    plot_df = data.frame("category" = temp_df$Term, 
                         'log_adj_p' = -log(temp_df$Adjusted.P.value))
    plot_df = plot_df[order(plot_df$log_adj_p, decreasing = TRUE), ]
    plot_df = plot_df[1:15, ]
    plot_df = plot_df[complete.cases(plot_df), ]
    p<-ggplot(data=plot_df, aes(x=reorder(category, log_adj_p), y=log_adj_p)) +
      geom_bar(stat="identity") + 
      ylab("-log adj-p") +
      xlab("Terms")+
      ggtitle(analysis_name) + 
      theme_bw() +
      coord_flip()
    ggsave(filename = paste0("significant_terms_", analysis_name, '.png'), width = 12, height = 10)
  }
})

