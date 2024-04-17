

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
