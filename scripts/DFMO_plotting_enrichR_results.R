library(ggplot2)

enrichment_results = read.csv("../output/DFMO_control_analyis/human_DE_analysis/control_genes/Reactome_2016_enrichment_analysis.csv", row.names = 1)
enrichment_results = enrichment_results[1:40, ]

plot_df = data.frame("category" = enrichment_results$Term, 
                     'log_adj_p' = -log(enrichment_results$Adjusted.P.value))
plot_df = plot_df[order(plot_df$log_adj_p, decreasing = TRUE), ]
plot_df$category = stringr::str_remove_all(plot_df$category, " Homo sapiens.*")
#plot_df = plot_df[1:20, ]
plot_df = plot_df[complete.cases(plot_df), ]
plot_df$color = 'black'
blue_list = c("Cell Cycle, Mitotic", 
             "Cell Cycle", 
             "Mitotic G1-G1/S phases", 
             "Activation of ATR in response to replication stress", 
             "S Phase", 
             "G1/S Transition", 
             "M/G1 Transition", 
             "Regulation of DNA replication", 
             "G2/M Checkpoints")

plot_df[plot_df$category %in% blue_list, 'color'] = 'blue'

red_list = c("Cholesterol biosynthesis", 
              "Regulation of cholesterol biosynthesis by SREBP (SREBF)", 
             "Activation of gene expression by SREBF (SREBP)", 
             "Metabolism of lipids and lipoproteins")

plot_df[plot_df$category %in% red_list, 'color'] = 'red'

yellow_list = c("Apoptotic cleavage of cellular proteins", 
             "Apoptotic execution  phase", 
             "Apoptosis")

plot_df[plot_df$category %in% yellow_list, 'color'] = 'yellow'

p<-ggplot(data=plot_df, aes(x=reorder(category, log_adj_p), y=log_adj_p)) +
  geom_bar(stat="identity") + 
  ylab("-log adj-p") +
  xlab("Terms")+
  ggtitle('Enriched in Control') + 
  theme_bw() +
  theme(text = element_text(size = 20)) +
  coord_flip()
ggsave("../output/DFMO_control_analyis/human_DE_analysis/enrichr_black.png", plot = p, height = 9.5, width = 12)  

p<-ggplot(data=plot_df, aes(x=reorder(category, log_adj_p), y=log_adj_p, fill = color)) +
  geom_bar(stat="identity") + 
  ylab("-log adj-p") +
  xlab("Terms")+
  ggtitle('Enriched in Control') + 
  scale_fill_manual(values=c("#999999", "#2E9FDF", "#FC4E07", "#E7B800")) +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = 'none') +
  coord_flip()
ggsave("../output/DFMO_control_analyis/human_DE_analysis/enrichr_color.png", plot = p, height = 9.5, width = 12)  

sub_plot_df = plot_df[plot_df$color != 'black', ]
p<-ggplot(data=sub_plot_df, aes(x=reorder(category, log_adj_p), y=log_adj_p)) +
  geom_bar(stat="identity") + 
  ylab("-log adj-p") +
  xlab("Terms")+
  ggtitle('Enriched in Control') + 
  theme_bw() +
  theme(text = element_text(size = 20)) +
  coord_flip()
ggsave("../output/DFMO_control_analyis/human_DE_analysis/selected_enrichr_black.png", plot = p, height = 9.5, width = 12)  

p<-ggplot(data=sub_plot_df, aes(x=reorder(category, log_adj_p), y=log_adj_p, fill = color)) +
  geom_bar(stat="identity") + 
  ylab("-log adj-p") +
  xlab("Terms")+
  ggtitle('Enriched in Control') + 
  scale_fill_manual(values=c("#2E9FDF", "#FC4E07", "#E7B800")) +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = 'none') +
  coord_flip()
ggsave("../output/DFMO_control_analyis/human_DE_analysis/selected_enrichr_color.png", plot = p, height = 9.5, width = 12)  
