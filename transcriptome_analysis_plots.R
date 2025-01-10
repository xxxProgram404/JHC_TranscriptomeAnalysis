setwd("./data_analysis/") # Set working directory
options(scipen = 200) # Avoid scientific notation

# Load required libraries
library(readr)
library(tidyverse)
library(openxlsx)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(scales)
library(RColorBrewer)

# Volcano Plot =============================================================
transcript_data <- read.xlsx("data/transcriptome_data.xlsx", rowNames = FALSE, sheet = "Total_Data")

gene_exp <- transcript_data %>%
  select(contains("FPKM"), contains("gene_name")) %>%
  rename_all(~ sub("_FPKM", "", .))

deg <- read.xlsx("data/transcriptome_data.xlsx", rowNames = FALSE, sheet = "Volcano_Plot")
data <- left_join(deg, gene_exp, by = "gene_name") %>%
  filter(!`Condition1_vs_Condition2_FDR` %in% "--") %>%
  mutate(FDR = as.numeric(`Condition1_vs_Condition2_FDR`),
         log2FC = as.numeric(`Condition1_vs_Condition2_log2FC`))

volcano <- ggplot(data, aes(x = log2FC, y = -log10(FDR), color = `Condition1_vs_Condition2_Regulation`)) +
  geom_point(alpha = 1, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  geom_vline(xintercept = c(-1, 1), linetype = 'dashed') +
  scale_color_manual(values = c('up' = '#F03214', 'down' = '#1E4198', 'normal' = '#778899')) +
  xlim(-10, 10) +
  xlab("log2 Fold Change") + ylab("-log10(FDR)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = 28, hjust = 0.5),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        legend.position = c(0.88, 0.88),
        legend.title = element_blank(),
        legend.text = element_text(size = 21)) +
  guides(colour = guide_legend(override.aes = list(size = 10)))

# Save volcano plot
ggsave('results/volcano_plot.pdf', volcano, width = 8, height = 7)

# Heatmap for Differential Genes ============================================
# Prepare expression data
deg_heatmap <- read.xlsx("data/transcriptome_data.xlsx", rowNames = FALSE, sheet = "Heatmap_Differential_Genes")
data_heatmap <- left_join(deg_heatmap, gene_exp, by = "gene_name")
rownames(data_heatmap) <- data_heatmap$`#ID`
map <- data_heatmap[, 8:17] %>% filter(apply(., 1, function(x) sd(x) != 0))

# Draw heatmap
heatmap_diff <- pheatmap(map,
                         scale = "row",
                         angle_col = 315,
                         color = colorRampPalette(c("#B1AFD7", "white", "#F03214"))(100),
                         fontsize_row = 12,
                         fontsize_col = 12,
                         show_rownames = FALSE)
# Save heatmap
ggsave('results/heatmap_differential_genes.pdf', heatmap_diff, width = 8, height = 7)

# KEGG Enrichment Plot ======================================================
kegg_data <- read.xlsx("data/transcriptome_data.xlsx", rowNames = FALSE, sheet = "KEGG_Enrichment") %>%
  arrange(-Count) %>%
  mutate(Description = factor(Description, levels = rev(Description)),
         Rich_Factor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

kegg_plot <- ggplot(kegg_data, aes(x = Rich_Factor, y = Description)) +
  geom_point(aes(size = Count, colour = pvalue)) +
  scale_color_distiller(palette = "RdYlGn") +
  labs(title = 'Enriched KEGG Pathways',
       x = 'Rich Factor',
       y = 'Description') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 10))

# Save KEGG plot
ggsave('results/kegg_enrichment.pdf', kegg_plot, width = 8, height = 7)

# GO Enrichment Plot ========================================================
go_data <- read.xlsx("data/transcriptome_data.xlsx", rowNames = FALSE, sheet = "GO_Enrichment") %>%
  arrange(ONTOLOGY, -Count)
COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")

GO_term_order <- factor(as.integer(rownames(go_data)), labels = go_data$Description)
go_plot <- ggplot(go_data, aes(x = Count, y = GO_term_order, fill = ONTOLOGY)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = COLS) +
  theme_bw() +
  coord_flip() +
  xlab("Number of Genes") + ylab("GO Term") +
  labs(title = "Enriched GO Terms", fill = "Ontology") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 12),
        legend.position = 'top',
        plot.title = element_text(hjust = 0.5, size = 20))

# Save GO plot
ggsave('results/go_enrichment.pdf', go_plot, width = 15, height = 10)

# General Heatmap ===========================================================
general_id <- read.xlsx("data/general_heatmap_ids.xlsx", rowNames = FALSE, colNames = FALSE, sheet = "General_Heatmap") %>%
  select(X1) %>% unique()

general_map <- left_join(general_id, gene_exp, by = c("X1" = "gene_name"))
rownames(general_map) <- general_map$X1
map <- general_map[, -1] %>% filter(apply(., 1, function(x) sd(x) != 0))

heatmap_general <- pheatmap(map,
                            scale = "row",
                            angle_col = 315,
                            color = colorRampPalette(c("#B1AFD7", "white", "#F03214"))(100),
                            fontsize_row = 12,
                            fontsize_col = 12,
                            show_rownames = FALSE)
# Save general heatmap
ggsave('results/general_heatmap.pdf', heatmap_general, width = 8, height = 7)
