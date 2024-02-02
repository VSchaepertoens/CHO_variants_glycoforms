library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(viridis)


# load_data ---------------------------------------------------------------

data <- read_csv('data/RelQuantIntact.csv') %>%
  column_to_rownames(var = "...1")

# plot heatmap of raw data ------------------------------------------------

data.matrix <- as.matrix(data)

Heatmap(data.matrix)
Heatmap(data.matrix, col = plasma(100))
Heatmap(data.matrix, col = rev(rainbow(10)))

# calculate z-score & plot heatmap -------------------------------------------------------

scaled.data.matrix = t(scale(t(data.matrix))) # for scaling by row  

#check for sanity
mean(data.matrix[1,])
sd(data.matrix[1,])
(data.matrix[1] - mean(data.matrix[1,]))/sd(data.matrix[1,])
(data.matrix[1,2] - mean(data.matrix[1,]))/sd(data.matrix[1,])

min(scaled.data.matrix)
max(scaled.data.matrix)
col_fun = colorRamp2(c(-1.3, 0, 1.6), c("green", "white", "red"))
Heatmap(scaled.data.matrix, col = col_fun)
Heatmap(scaled.data.matrix, col = plasma(100))
Heatmap(scaled.data.matrix, col = rev(rainbow(10)))


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
  mutate(percent = (mean_peak_area/sum(mean_peak_area)))

#convert 'subunit' to factor and specify level order
data_averaged$subunit <- factor(data_averaged$subunit, levels = c('Intact', 'LC/LC', 'LC'))

# plot stacked bar chart -------------------------------------------------

ggplot(data_averaged, aes(y = cell_variant, x = mean_peak_area, fill = subunit)) + 
  geom_bar(stat = "identity", position = "fill") +
  xlab("peak_area (%)") +
  scale_fill_brewer(palette = "Accent") +
  scale_y_discrete(limits = rev)

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
 
ggsave("figures/stacked_bar/stacked_bar_percent.png")






















