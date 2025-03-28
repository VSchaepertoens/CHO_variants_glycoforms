library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(viridis)


# load_data ---------------------------------------------------------------

data <- read_csv('data/1_subunit_quantification/RelQuantIntact01.csv') %>%
  column_to_rownames(var = "...1")


# build correct table -----------------------------------------------------

data <- data %>%
  rownames_to_column("cell_variant") %>%
  pivot_longer(cols = colnames(data),
               names_to = c("subunit","replicate"),
               names_sep = "_",
               values_to = "peak_area")

data$replicate <- str_replace(data$replicate, "PeakArea", "")

#calculate group means
data_averaged <- data %>%
  group_by(cell_variant,subunit) %>%
  summarise(mean_peak_area = mean(peak_area)) %>%
  group_by(cell_variant) %>%
  mutate(percent = (mean_peak_area/sum(mean_peak_area))) %>%
  mutate(subunit = str_replace_all(subunit, c('LC/LC' = 'LC-LC'))) %>%
  mutate(cell_variant = factor(cell_variant, levels = c("A25_1", "A25_2", 
                                                        "A2_1", "A2_2", 
                                                        "A3_1", "A3_2", 
                                                        "A4_1", "A4_2", 
                                                        "A8_1", "A8_2", 
                                                        "A16_1", "A16_2", 
                                                        "A19_1", "A19_2"))) %>%
  filter(cell_variant %in% c( "A2_1", "A2_2", 
                              "A3_1", "A3_2", 
                              "A4_1", "A4_2", 
                              "A8_1", "A8_2", 
                              "A16_1", "A16_2", 
                              "A19_1", "A19_2"))


#convert 'subunit' to factor and specify level order ## important for plotting order
data_averaged$subunit <- factor(data_averaged$subunit, levels = c('Intact', 'LC-LC', 'LC'))

# plot stacked bar chart -------------------------------------------------

ggplot(data_averaged, aes(y = cell_variant, x = mean_peak_area, fill = subunit)) + 
  geom_bar(stat = "identity", position = "fill", width = .7) +
  xlab("peak_area (%)") +
  ylab("CHO cell variant") +
  scale_fill_brewer(palette = "Accent") +
  scale_y_discrete(limits = rev) +
  theme(text = element_text(size = 9, 
                            # face = "bold", 
                            family = "sans"),
        axis.text = element_text(colour = "black"),
        panel.background = element_blank(),
        axis.text.y = element_text(margin = margin(r = 0)),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
  )

ggsave("figures/1_subunit_quantification/figure_8.png",        
       width = 8.89,
       height = 8.89,
       units = c("cm"),
       dpi = 600)

# ggsave("figures/1_subunit_quantification/figure_8.svg",        
#        width = 8.89,
#        height = 8.89,
#        units = c("cm"),
#        dpi = 600)

#plot stacked barchart with labels
ggplot(data_averaged, aes(x = cell_variant, y = mean_peak_area, fill = subunit)) + 
  geom_col(position = "fill") +
  geom_text(aes(label = ifelse(mean_peak_area == 0, "", round(percent*100,0))),
            position = position_fill(vjust = 0.5)) +
  geom_hline(yintercept = 0, linewidth = .5) +
  labs(y = "peak area (%)") +
  labs(x = "CHO cell variant") +
  scale_fill_brewer(palette = "Accent") +
  coord_flip() +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank()
    ) 
 
# ggsave("figures/stacked_bar/stacked_bar_percent_all_SC.png", 
#        width = 8.89,
#        height = 6.5,
#        units = c("cm"),
#        dpi = 600)
