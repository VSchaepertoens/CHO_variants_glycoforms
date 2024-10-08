library(tidyverse, warn.conflicts = FALSE)
library(RColorBrewer)
library(circlize)
library("scales")    
library(ComplexHeatmap)

# load an overview table of data & analysis paths -------------------------

samples_table <- read_csv("data/Dec_2023/overview_hc_dimer_pngase_merged.csv") 
 

# load abundances using a for loop  ---------------------------------------

abundance_data <- NULL
for (i in 1:nrow(samples_table)) {
  file_path <- paste0(samples_table[i, "analysis_path"], "/frac_ab_tb_cs50.csv")

abundance_data <- rbind(abundance_data,
                        read_delim(file_path) %>%
                          mutate(CHO_cell_variant_bio_replicate = samples_table$CHO_cell_variant_bio_replicate[i],
                                 tech_replicate = samples_table$tech_replicate[i]
                                 )
                        )
}

# plot heatmap ------------------------------------------------------------

## wrangle dataframe into matrix -------------------------------------------
data.matrix <- abundance_data %>%
  mutate(sample_name = paste(CHO_cell_variant_bio_replicate,tech_replicate, sep = "_")) %>%
  select("sample_name", "modcom_name", "frac_ab") %>%
  spread(key = sample_name, value = frac_ab) %>%
  column_to_rownames("modcom_name") %>%
  as.matrix()


## calculate z-score & plot heatmap -------------------------------------

scaled.data.matrix = t(scale(t(data.matrix))) # for scaling by row  

#check for sanity
mean(data.matrix[1,])
sd(data.matrix[1,])
(data.matrix[1] - mean(data.matrix[1,]))/sd(data.matrix[1,])
(data.matrix[1,2] - mean(data.matrix[1,]))/sd(data.matrix[1,])


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
png(filename = "figures/Dec_2023/heatmap_scaled_hc_dimer.png",    
    height = 2400,
    width = 2400,
    units = "px",
    res = 300)

Heatmap(scaled.data.matrix,
        col = f1,
        cluster_rows = TRUE,
        rect_gp = gpar(col = "white", lwd = 2))

dev.off()

# calculate mean and sd and plot ------------------------------------------
abundance_data_averaged <- abundance_data %>% 
  group_by(modcom_name, CHO_cell_variant_bio_replicate) %>%
  summarise(frac_abundance = mean(frac_ab),
            error = sd(frac_ab)) %>%   
  mutate(modcom_name = factor(modcom_name, levels = c("2xHex","1xHex","Lys/Lys","none/Lys","none/none")),
            ) 

save(abundance_data, 
     abundance_data_averaged, 
     file = "analysis/Dec_2023/abundance_data_hc_dimer.RData")


# plot bar chart ----------------------------------------------------------

abundance_data_averaged %>%
  ggplot(aes(modcom_name, frac_abundance)) +
  geom_col(
    aes(y = frac_abundance, fill = CHO_cell_variant_bio_replicate),
    position = position_dodge(.9),
  ) +
  geom_errorbar(
    aes(
      ymin = frac_abundance - error,
      ymax = frac_abundance + error,
      group = CHO_cell_variant_bio_replicate
    ),
    position = position_dodge(.9),
    width = .5,
    linewidth = .25
  ) +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(breaks = c(0,10,20,30,35),
                     labels = number_format(accuracy = 1)) +
  xlab("") +
  ylim(0, 100) +
  ylab("fractional abundance (%)") +
  labs(title = "Lysine and glucose variants - hc_dimer") +
  geom_hline(yintercept = 0, linewidth = .35) +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 16, 
                            face = "bold", 
                            family = "sans"),
        axis.text = element_text(colour = "black"),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
)


ggsave(filename = "figures/Dec_2023/frac_ab_barplot_hc_dimer.png",    
       height = 160,
       width = 160,
       units = "mm",
       dpi = 600)

