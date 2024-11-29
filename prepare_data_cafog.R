#prepare data for cafog analysis
library(tidyverse)
library(dplyr)
library(stringr)
library(fs)

# load abundance data -----------------------------------------------------

load("analysis/2_nglycans_quantification/2_1_not_corrected_glycans/intact_output_tables/abundance_data.RData")
glycosylation <-  abundance_data_averaged

# load glycation data -----------------------------------------------------


load("analysis/2_nglycans_quantification/2_2_pngaseF_cpb/pngase_cpb_output_tables/abundance_data_cpb_pngase.RData")

glycation <-  abundance_data_averaged


# Define composition mapping ------------------------------------------------------------------

composition_mapping <- list(
  G1F = "4 HexNAc, 4 Hex, 1 Fuc",
  G2F = "5 HexNAc, 4 Hex",
  S1G1F = "4 HexNAc, 4 Hex, 1 Fuc, 1 Sa",
  none = "0 Hex"
)

# for each coef make files for cafog analysis --------------------------------------
coefs <-  unique(glycosylation$CHO_cell_variant_bio_replicate)

fs::dir_create(paste0("analysis/2_nglycans_quantification/2_3_cafog_corrected_glycans/",coefs))

#for loop to create cafog files for each CHO_cell_variant_bio_repliacte
for (coef in coefs) {
  print(coef)
  glycosylation %>%
    filter(CHO_cell_variant_bio_replicate == coef) %>%
    select(modcom_name, frac_abundance, error) %>%
    rename(`#glycoform` = modcom_name,
           abundance = frac_abundance) %>%
    write_csv(paste0("analysis/2_nglycans_quantification/2_3_cafog_corrected_glycans/",coef,"/glycosylation.csv"),
              col_names = TRUE)

  glycation %>%
    filter(CHO_cell_variant_bio_replicate == coef) %>%
    select(modcom_name, frac_abundance, error) %>%
    rename(`#count` = modcom_name,
           abundance = frac_abundance) %>%
    mutate(`#count` = as.character(`#count`)) %>%
    mutate(`#count` = str_replace_all(`#count`, c("1xHex" = "1", "2xHex" = "2", "none/none" = "0"))) %>%
    write_csv(paste0("analysis/2_nglycans_quantification/2_3_cafog_corrected_glycans/",coef,"/glycation.csv"),
              col_names = TRUE)

  glycosylation %>%
    filter(CHO_cell_variant_bio_replicate == coef) %>%
    select(modcom_name) %>%
    separate(modcom_name, into = c("glycoform_1", "glycoform_2"), sep = "/") %>%
    pivot_longer(cols = c("glycoform_1", "glycoform_2"), names_to = "names", values_to = "glycoforms") %>%
    select(glycoforms) %>%
    unique()  %>%
    mutate(
      composition = case_when(
        # glycoforms == "G1F" ~ composition_mapping[["G1F"]],
        glycoforms == "none" ~ composition_mapping[["none"]],
        glycoforms == "S1G1F" ~ composition_mapping[["S1G1F"]],
        TRUE ~ NA_character_  # Handle unmatched cases
      )) %>%
    write_csv(paste0("analysis/2_nglycans_quantification/2_3_cafog_corrected_glycans/",coef,"/glycan_library.csv"),
                        col_names = TRUE)
} 

## Remove NAs from glycan library manually
## Run CAFOG analysis using subprocess_cafog.ipynb from Anaconda --> vs studio
## Continue with plotting the corrected results --> plot_cafog_corrected.R
