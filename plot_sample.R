library(tidyverse, warn.conflicts = FALSE)
library(RColorBrewer)
#install.packages("scales")                  # Install scales package
library("scales")    

# load an overview table of data & analysis paths -------------------------

samples_table <- read_delim("data/Jan_2024/overview_sample_merged.csv", delim = ",") 
 

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


# calculate mean and sd and plot ------------------------------------------
abundance_data_averaged <- abundance_data %>% 
  group_by(modcom_name, CHO_cell_variant_bio_replicate) %>%
  summarise(frac_abundance = mean(frac_ab),
            error = sd(frac_ab)) 


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
  ylim(0, 40) +
  ylab("fractional abundance (%)") +
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


ggsave(filename = "figures/Jan_2024/frac_ab_barplot.png",    
       height = 160,
       width = 160,
       units = "mm",
       dpi = 600)

