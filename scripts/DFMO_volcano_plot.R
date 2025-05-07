library(EnhancedVolcano)
library(DESeq2)
library(ggplot2)

###### make the plot for the PDX data #####
dds = readRDS('../output/DFMO_control_analyis/human_DE_analysis/human_DESeq2_object.rds')
res <- results(dds)
res_df = as.data.frame(res)
res_df = res_df[complete.cases(res_df), ]

header_tab = read.csv("../input/DFMO_control_alignment/header_tab.csv", row.names = 1)
header_tab = header_tab[header_tab$ensembl_id %in% rownames(res_df), ]
rownames(header_tab) = header_tab$ensembl_id
res_df$gene = header_tab[rownames(res_df), 'gene_id']
res_df = res_df[order(res_df$pvalue), ]
res_df = res_df[!duplicated(res_df$gene), ]

p = EnhancedVolcano(res_df,
                    lab = res_df$gene,
                    x = 'log2FoldChange',
                    title = 'DFMO vs control',
                    subtitle = "Deseq2 results",
                    y = 'pvalue',
                    caption = bquote(~Log[2]~ "fold change cutoff, 1 ; p-value cutoff, 0.05"),
                    drawConnectors = TRUE,
                    widthConnectors = 0.75,
                    pCutoff = 0.05, FCcutoff = 1)
ggsave(filename = '../output/DFMO_control_analyis/human_DE_analysis/volcano_plot.png', height = 10, width = 14)

###### check expression direction ##### 
norm_exp = read.csv('../output/DFMO_control_analyis/human_DE_analysis/human_norm_exp.csv')
apply(norm_exp[norm_exp$X == 'SHC1', grepl('E', colnames(norm_exp))], FUN = mean, MARGIN = 1)
apply(norm_exp[norm_exp$X == 'SHC1', grepl('C', colnames(norm_exp)) | grepl('L', colnames(norm_exp))], FUN = mean, MARGIN = 1)

grepl('E', colnames(norm_exp))
