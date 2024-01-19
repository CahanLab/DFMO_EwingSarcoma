library(DESeq2)
library(stringr)
library(ggplot2)
header_tab = read.csv("../input/DFMO_control_alignment/reference_genomes/hg38_mm10_header.txt", sep = ' ', header = FALSE)
header_tab = header_tab[, c(1, 7)]
colnames(header_tab) = c("ensembl_id", 'gene_id')
header_tab$ensembl_id = str_remove_all(header_tab$ensembl_id, ">")
header_tab$gene_id = str_remove_all(header_tab$gene_id, 'gene_symbol:')

dds = readRDS("../output/DFMO_control_analyis/human_DE_analysis/human_DESeq2_object.rds")
res <- results(dds)

header_tab = header_tab[header_tab$ensembl_id %in% rownames(res), ]
rownames(header_tab) = header_tab$ensembl_id

rownames(res) = header_tab[rownames(res), 'gene_id']
res_df = data.frame(res)
res_df$log10_padj = -log10(res_df$padj) 
res_df$genes = rownames(res_df)
res_df$Type = NA
res_df[res_df$log2FoldChange > 0, 'Type'] = 'Up in DFMO'
res_df[res_df$log2FoldChange < 0, 'Type'] = 'Up in control'

##### look at ferroptosis genes #####
gene_sets = read.csv("../input/ferroptosis_metabolic_genesets//ferroptosis.txt", header = FALSE)
gene_sets = gene_sets$V1
gene_sets = intersect(gene_sets, rownames(res))
sub_res_df = res_df[gene_sets, ]

p<-ggplot(data=sub_res_df, aes(x=reorder(genes, abs(log10_padj)), y=log10_padj, fill = Type)) +
  geom_bar(stat="identity") + coord_flip() + geom_hline(yintercept = -log10(0.5)) + theme_bw() + 
  xlab('ferroptosis genes') +
  ylab('-log10 adj p-value') + 
  scale_fill_brewer(palette = 'Set2') +
  ggtitle('ferroptosis genes') + 
  theme(text = element_text(size = 20))
ggsave(filename = '../output/DFMO_control_analyis/human_DE_analysis/ferroptosis_genes_bar.png', width = 7, height = 6)

p<-ggplot(data=sub_res_df, aes(x=reorder(genes, abs(log2FoldChange)), y=log2FoldChange, fill = Type)) +
  geom_bar(stat="identity") + coord_flip() + theme_bw() +
  xlab('ferroptosis genes') +
  ylab('log2fold Change') + 
  ggtitle('ferroptosis genes') + 
  theme(text = element_text(size = 20)) + 
  scale_fill_brewer(palette = 'Set2') 
ggsave(filename = '../output/DFMO_control_analyis/human_DE_analysis/ferroptosis_genes_bar_log2Fold.png', width = 7, height = 6)

##### look at metabolic genes #####
gene_sets = read.csv("../input/ferroptosis_metabolic_genesets/metabolic.txt", header = FALSE)
gene_sets = gene_sets$V1
gene_sets = intersect(gene_sets, rownames(res))
sub_res_df = res_df[gene_sets, ]

p<-ggplot(data=sub_res_df, aes(x=reorder(genes, abs(log10_padj)), y=log10_padj, fill = Type)) +
  geom_bar(stat="identity") + coord_flip() + geom_hline(yintercept = -log10(0.5)) + theme_bw() + 
  xlab('metabolic genes') +
  ylab('-log10 adj p-value') + 
  ggtitle('metabolic genes')
ggsave(filename = '../output/DFMO_control_analyis/human_DE_analysis/metabolic_genes_bar.png', width = 8, height = 6)

p<-ggplot(data=sub_res_df, aes(x=reorder(genes, abs(log2FoldChange)), y=log2FoldChange, fill = Type)) +
  geom_bar(stat="identity") + coord_flip() + 
  xlab('metabolic genes') +
  ylab('log2fold Change') + 
  ggtitle('metabolic genes')
ggsave(filename = '../output/DFMO_control_analyis/metabolic_genes_bar_log2Fold.png', width = 8, height = 6)

##### look at the genes individually (NOT USED) ######
human_counts = read.csv("../output/DFMO_control_analyis/human_DE_analysis/human_norm_exp.csv")
human_counts = human_counts[!duplicated(human_counts$X), ]
rownames(human_counts) = human_counts$X
human_counts$X = NULL

meta_tab = read.csv("../output/DFMO_control_analyis/human_mouse_counts/meta_tab.csv", row.names = 1)

make_plot <- function(exp_df, meta_tab, gene_interest) { 
  exp_df = exp_df[, rownames(meta_tab)]
  plot_df = data.frame("type" = meta_tab$treated, 
                       "norm_exp" = unlist(exp_df[gene_interest, ]))
  
  p = ggplot(plot_df, aes(x=type, y=norm_exp, fill = type)) +
    geom_boxplot()+
    scale_fill_brewer(palette="Set2") +
    theme_classic() + 
    ggtitle(gene_interest) 
  return(p)
}
gene_interest = 'IDH1'
plot = make_plot(human_counts, meta_tab, gene_interest)

gene_interest = 'HMGCR'
make_plot(human_counts, meta_tab, gene_interest)
