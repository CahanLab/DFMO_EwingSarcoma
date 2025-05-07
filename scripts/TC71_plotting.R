library(ggplot2)
library(DESeq2)
library(stringr)
library(enrichR)
library(pheatmap)
library(singleCellNet)

# Primary Author: Yuqi Tan 

#### load in the data #####
exp_dat = utils_loadObject("../input/TC71_alignment/expDatGood_091319.rda")
samp_tab = utils_loadObject("../input/TC71_alignment/sampTabGood_091319.rda")

samp_tab = samp_tab[samp_tab$description1 != 'primary_tumor', ]
samp_tab = samp_tab[order(samp_tab$description1), ]
exp_dat = exp_dat[, rownames(samp_tab)]

exp_dat = round(exp_dat)
dds <- DESeqDataSetFromMatrix(countData=exp_dat, 
                              colData=samp_tab, 
                              design=~description1, tidy = FALSE)
#dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res = results(dds)

raw_res_df = data.frame(res)
res_df = data.frame(res)
res_df = res_df[!is.na(res_df$padj), ]
res_df = res_df[res_df$padj < 0.05, ]

gsea_results = read.csv("../output/TC71_analysis/gsea/c2.fgseaRes.csv", row.names = 1)

make_plots_data <- function(gsea_results, dds, pathway_name, save_path) {
  dir.create(file.path(save_path, pathway_name))
  target_genes = gsea_results[gsea_results$pathway == pathway_name, 'leadingEdge']
  
  target_genes = stringr::str_split(target_genes, ', ')[[1]]
  
  norm_exp = counts(dds, normalized = TRUE)
  my_sample_col <- data.frame(sample = samp_tab$description1)
  row.names(my_sample_col) <- rownames(samp_tab)
  
  norm_exp = norm_exp[, rownames(my_sample_col)]
  sub_norm_exp = norm_exp[target_genes, ]
  sub_norm_exp = t(scale(t(sub_norm_exp)))
  
  withr::with_dir(file.path(save_path, pathway_name), {
    png(filename = 'heatmap.png', width = 8000, height = 5000, res = 600)
    pheatmap(sub_norm_exp,
             annotation_col = my_sample_col, 
             legend = TRUE, cluster_cols = FALSE, cluster_rows = TRUE)
    dev.off()
  })
  
  for(gene in target_genes) {
    print(gene)
    plot_df = data.frame('norm_expression' = log(norm_exp[gene, ] + 1), 
                         'type' = samp_tab$description1)
    p <- ggplot(plot_df, aes(x=type, y=norm_expression, fill = type)) + 
      geom_violin() + geom_boxplot(width=0.05, fill="white") + theme_bw() + 
      scale_fill_brewer(palette="Dark2") +
      xlab("Data Type") + 
      ylab("Log Normalized Expression") + 
      ggtitle(gene) + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    withr::with_dir(file.path(save_path, pathway_name), {
      ggsave(paste0(gene, "_exp.png"), plot = p, width = 6, height = 4)
    })
  }
  
}

make_plots_data(gsea_results, dds, 'KEGG_PENTOSE_PHOSPHATE_PATHWAY', save_path = "../output/TC71_analysis/gsea")  
make_plots_data(gsea_results, dds, 'KEGG_GLUTATHIONE_METABOLISM', save_path = "../output/TC71_analysis/gsea")  
make_plots_data(gsea_results, dds, 'KEGG_FATTY_ACID_METABOLISM', save_path = "../output/TC71_analysis/gsea")  

##### make individual plots for ACSL3 and PTGS2 #####
pathway_name = 'ACSL3_PTGS2'
dir.create(file.path('../output/TC71_analysis/gsea', pathway_name))

target_genes = c("ACSL3", 'PTGS2')

norm_exp = counts(dds, normalized = TRUE)
my_sample_col <- data.frame(sample = samp_tab$description1)
row.names(my_sample_col) <- rownames(samp_tab)

norm_exp = norm_exp[, rownames(my_sample_col)]
sub_norm_exp = norm_exp[target_genes, ]
sub_norm_exp = t(scale(t(sub_norm_exp)))

withr::with_dir(file.path('../output/TC71_analysis/gsea', pathway_name), {
  png(filename = 'heatmap.png', width = 10000, height = 5000, res = 1200)
  pheatmap(sub_norm_exp,
           annotation_col = my_sample_col, 
           legend = TRUE, cluster_cols = FALSE, cluster_rows = FALSE)
  dev.off()
})

for(gene in target_genes) {
  print(gene)
  plot_df = data.frame('norm_expression' = log(norm_exp[gene, ] + 1), 
                       'type' = samp_tab$description1)
  p <- ggplot(plot_df, aes(x=type, y=norm_expression, fill = type)) + 
    geom_violin() + geom_boxplot(width=0.05, fill="white") + theme_bw() + 
    scale_fill_brewer(palette="Dark2") +
    xlab("Data Type") + 
    ylab("Log Normalized Expression") + 
    ggtitle(gene) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  withr::with_dir(file.path('../output/TC71_analysis/gsea', pathway_name), {
    ggsave(paste0(gene, "_exp.png"), plot = p, width = 6, height = 4)
  })
}

# plot DECR1
norm_exp = counts(dds, normalized = TRUE)
plot_df = data.frame('norm_expression' = log(norm_exp['DECR1', ] + 1), 
                     'type' = samp_tab$description1)
p <- ggplot(plot_df, aes(x=type, y=norm_expression, fill = type)) + 
  geom_violin() + geom_boxplot(width=0.05, fill="white") + theme_bw() + 
  scale_fill_brewer(palette="Dark2") +
  xlab("Data Type") + 
  ylab("Log Normalized Expression") + 
  ggtitle('DECR1') + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(file.path('../output/TC71_analysis/gsea', pathway_name, paste0("DECR1_exp.png")), plot = p, width = 6, height = 4)

###### remake the manuscript figures ###### 
pathway_name = 'Supp_figure_S10'
dir.create(file.path('../output/TC71_analysis/gsea', pathway_name))

norm_exp = counts(dds, normalized = TRUE)
target_genes = c("DECR1", "GSTM2", 'GSTM3', 'ACSL3', 'PTGS2', 'G6PD')
for(gene in target_genes) {
  plot_df = data.frame('norm_expression' = norm_exp[gene, ], 
                       'type' = samp_tab$description1)
  p <- ggplot(plot_df, aes(x=type, y=norm_expression, fill = type)) + 
    geom_violin() + geom_boxplot(width=0.05, fill="white") + theme_bw() + 
    scale_fill_brewer(palette="Dark2") +
    xlab("Data Type") + 
    ylab("Deseq2 Normalized Expression") + 
    ggtitle(paste0(gene, " p-value: ", round(raw_res_df[gene, 'pvalue'], 2))) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(file.path('../output/TC71_analysis/gsea', pathway_name, paste0(gene, "_exp.png")), plot = p, width = 6, height = 4)
}

