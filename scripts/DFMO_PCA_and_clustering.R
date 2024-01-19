library(ggplot2)
library(DESeq2)
library(amap)
library(factoextra)

dds = readRDS("../output/DFMO_control_analyis/human_DE_analysis/human_DESeq2_object.rds")
human_counts = counts(dds, normalized = TRUE)
human_meta = read.csv("../output/DFMO_control_analyis/human_mouse_counts/meta_tab.csv", row.names = 1)
human_meta$id = NA
human_meta[human_meta$treated == 'DFMO', 'id'] = paste0("DFMO", seq(1, 4))
human_meta[human_meta$treated == 'control', 'id'] = paste0("control", seq(1, 6))

##### make plots for PCA plots ######
variance_gene = apply(human_counts, FUN = var, MARGIN = 1)
variance_gene = sort(variance_gene, decreasing = TRUE)

filtered_counts = human_counts[names(variance_gene[1:2000]), ]
PCA_coord = prcomp(t(filtered_counts))
PCA_coord = PCA_coord[['x']]

PCA_coord = as.data.frame(PCA_coord)
human_meta = human_meta[rownames(PCA_coord), ]
rownames(PCA_coord) = human_meta$id
PCA_coord$sample = human_meta$treated
PCA_coord$id = human_meta$id
p = ggplot(PCA_coord, aes(x=PC1, y=PC2, color = sample)) + 
  geom_point(size = 4) + 
  geom_text(
    label=rownames(PCA_coord), 
    color = 'black',
    nudge_x = 0, vjust = -0.9, 
    check_overlap = T
  ) +
  scale_colour_brewer(palette = 'Set2') +
  theme_bw() +
  theme(axis.text = element_blank(), text = element_text(size = 20))

ggsave("../output/DFMO_control_analyis/human_DE_analysis/PCA_plot.png", plot = p, height = 5, width = 8)  

##### plot out clustering (NOT USED IN MANUSCRIPT) ######
human_counts = human_counts[, rownames(human_meta)]
colnames(human_counts) = human_meta$id
my_correlation = Dist(t(human_counts), method = "pearson", nbproc = 2, diag = FALSE, upper = FALSE)
clusters <- hclust(my_correlation, method = 'average')
p = fviz_dend(clusters, k = 3,                 # Cut in four groups
          cex = 0.5,                 # label size
          k_colors = c(RColorBrewer::brewer.pal(3, 'Set2')[1], RColorBrewer::brewer.pal(3, 'Set2')[2], RColorBrewer::brewer.pal(3, 'Set2')[1]),
          color_labels_by_k = FALSE,  # color labels by groups
          ggtheme = theme_minimal()     # Change theme
)
p = p + ylab("Height") + 
  ggtitle("Hierarchal Clustering ") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("../output/DFMO_control_analyis/human_DE_analysis/clustering.png", plot = p, height = 4, width = 6)  
