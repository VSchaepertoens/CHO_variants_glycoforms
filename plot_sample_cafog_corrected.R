library(tidyverse)
library("scales")    
library(ComplexHeatmap)
library(circlize)

# load cafog corrected data -----------------------------------------------
# load abundances using a for loop  ---------------------------------------

abundance_data <- NULL
coefs <- c("A16_1", "A16_2", "A2_1", "A2_2", "A3_1", "A3_2", "A4_1", "A4_2", "A8_1", "A8_2", "A19_1", "A19_2")

for (coef in coefs) {
  file_path <- paste0("analysis/cafog/",coef,"/results.csv")
  
  abundance_data <- rbind(abundance_data,
                          read_csv(file_path,
                                   n_max = 12) %>%
                            mutate(CHO_cell_variant_bio_replicate = coef) %>%
                            {.}
  )
}

data_to_plot <- abundance_data %>%
  separate(
    glycoform,
    into = c("glycoform1", "glycoform2",  "glycoform3"),
    sep = "\\s+or\\s+",
    remove = FALSE
  ) %>%
  select(glycoform1, corr_abundance, corr_abundance_error, CHO_cell_variant_bio_replicate) %>%
  mutate(glycoform1 = str_replace_all(glycoform1, c("G0F/G2F" = "G1F/G1F", "G2F/none" = "none/G2F", "G1F/none" = "none/G1F", "G0F/none" = "none/G0F", "G0/none" = "none/G0"))) %>%
  mutate(glycoform1 = factor(glycoform1, levels = c("G1F/S1G1F","G2F/G2F","G1F/G2F","G1F/G1F","G0F/G1F","G0F/G0F","G0F/G0", "none/G2F", "none/G1F", "none/G0F", "none/G0","none/none")))


data_to_plot %>%
  ggplot(aes(glycoform1, corr_abundance)) +
  geom_col(
    aes(y = corr_abundance, fill = CHO_cell_variant_bio_replicate),
    position = position_dodge(.9),
  ) +
  geom_errorbar(
    aes(
      ymin = corr_abundance - corr_abundance_error,
      ymax = corr_abundance + corr_abundance_error,
      group = CHO_cell_variant_bio_replicate
    ),
    position = position_dodge(.9),
    width = .5,
    linewidth = .25
  ) +
  scale_fill_brewer(palette = "Paired") +
  xlab("") +
  ylim(0, 60) +
  ylab("fractional abundance (%)") +
  labs(title = "Hexose bias corrected glycoforms - intact") +
  geom_hline(yintercept = 0, linewidth = .35) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 16, 
                            face = "bold", 
                            family = "sans"),
        axis.text.y = element_text(colour = "black", hjust = 0.5),
        axis.text = element_text(colour = "black"),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
  )


ggsave(filename = "figures/Jan_2024/frac_ab_barplot_cafog_corrected_reordered.png",    
       height = 160,
       width = 160,
       units = "mm",
       dpi = 600)


# plot data as a heatmap --------------------------------------------------

data.matrix <- data_to_plot %>%
  select(glycoform1, corr_abundance, CHO_cell_variant_bio_replicate) %>%
  pivot_wider(names_from = CHO_cell_variant_bio_replicate, values_from = corr_abundance) %>%
  column_to_rownames("glycoform1") %>%
  arrange(c(7,8, 9, 10, 11, 12, 4, 3, 6, 5, 1, 2)) %>%
  as.matrix()
  
## calculate z-score & plot heatmap -------------------------------------

scaled.data.matrix = t(scale(t(data.matrix))) # for scaling by row  

#check for sanity
mean(data.matrix[1,])
sd(data.matrix[1,])
(data.matrix[1] - mean(data.matrix[1,]))/sd(data.matrix[1,])
(data.matrix[1,2] - mean(data.matrix[1,]))/sd(data.matrix[1,])


BASE_TEXT_SIZE_PT <- 5
ht_opt(
  simple_anno_size = unit(1.5, "mm"),
  COLUMN_ANNO_PADDING = unit(1, "pt"),
  DENDROGRAM_PADDING = unit(1, "pt"),
  HEATMAP_LEGEND_PADDING = unit(1, "mm"),
  ROW_ANNO_PADDING = unit(1, "pt"),
  TITLE_PADDING = unit(2, "mm"),
  heatmap_row_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_row_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_column_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  heatmap_column_names_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_labels_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_title_gp = gpar(fontsize = BASE_TEXT_SIZE_PT),
  legend_border = FALSE
)

#set the correct color scheme
min(scaled.data.matrix)
max(scaled.data.matrix)
f1 = colorRamp2(seq(-max(abs(scaled.data.matrix)),
                    max(abs(scaled.data.matrix)),
                    length = 9),
                c("seagreen4",
                  "seagreen3",
                  "seagreen2",
                  "seagreen1",
                  "gold",
                  "darkorchid1",
                  "darkorchid2",
                  "darkorchid3",
                  "darkorchid4"),
                space = "RGB")
#set the correct color scheme
png(filename = "figures/Jan_2024/heatmap_scaled_cafog_corrected_reordered.png",    
    height = 5,
    width = 5,
    units = "cm",
    res = 600)


Heatmap(scaled.data.matrix,
        col = f1,
        cluster_rows = FALSE,
        rect_gp = gpar(col = "white", lwd = 2),
        name = "z-score",
        row_gap = unit(2, "pt"),
        column_gap = unit(2, "pt"),
        width = unit(2, "mm") * ncol(scaled.data.matrix) + 5 * unit(2, "pt"), # to make each cell a square
        height = unit(2, "mm") * nrow(scaled.data.matrix) + 5 * unit(2, "pt"), # to make each cell a square
        show_row_names = TRUE
)

dev.off()
