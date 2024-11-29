## ---------------------------
##
## Script name: Quantification of fractional abundances of glycans found in 
## BOKU Herceptin project 2023
##
## Purpose of script: Using fragquaxi package to quantify glycans in samples
##
## Author: Dr. Veronika Schäpertöns
##
## Date Created: 01.02.2024 
##
## Copyright (c) Veronika Schäpertöns, 2023
## Email: veronika.schaepertoens@plus.ac.at
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(fs)
library(tidyverse, warn.conflicts = FALSE)
library(janitor, warn.conflicts = FALSE)
library(fragquaxi)

# define to either run analysis of sample mab or pngase digested
product <- "pngase" #OR product <- "pngase" 

directory <- "data/3_lysine_glycation_quantification/pngase_input_data"

# define mab_sequence  --------------------------------

mab_sequence <- "mab_sequence/trastuzumab_noLysine.txt"

proteins <- define_proteins(
  mab = mab_sequence,
  .disulfides = 16
)

# depending on product, define regexp pattern and modcoms -----------------

if (product == "sample") {
  #specify pattern to match regexp for filename
  #pattern <-  "ambr.*\\.mzML"
  pattern <-  "CpB.*\\.mzML"
  
  modcoms <- tribble(
    ~modcom_name, ~Hex, ~HexNAc, ~Fuc, ~Neu5Ac,
    "none/none",     0,       0,    0,     0,
    "G0F/G0F",       6,       8,    2,     0,
    "none/G0",       3,       4,    0,     0,
    "none/G0F",       3,       4,    1,     0,
    "G0F/G1F",       7,       8,    2,     0, 
    "G1F/G1F",       8,       8,    2,     0,
    "none/G1F",       4,       4,    1,     0,
    "G1F/G2F",       9,       8,    2,     0,
    "G2F/G2F",      10,       8,    2,     0,
    "none/G2F",      5,       4,    1,     0,
    "G0F/G0",        6,       8,    1,     0,
    # "G0/G0",         6,       8,    0,     0,
    #other:
    # "M5/M5",        10,       4,    0,     0,
    "G1F/S1G1F",      9,       8,    2,     1,     #4,
    # "G2F/S1G1F",     10,       8,    2,     1,     #4,
    # "G0F/M5",         8,       6,    1,     0,     #4,
    # "G0/M5",          8,       6,    0,     0,     #4,
    # "G0F",           3,       4,    1,     0,
    # "G0",            3,       4,    0,     0,
    # "G1F",           4,       4,    1,     0,
    # "G2F",           5,       4,    1,     0,
    #"G0F/M5",             8,       6,    1,     0,
    
  ) %>% 
    define_ptm_compositions()
   
 } else if (product == "pngase")  {
   #specify pattern to match regexp for filename 
   pattern <- "pngase.*\\.mzML"
   
   #specify modifications composition
   modcoms <- tribble(
     ~modcom_name, ~Hex, ~HexNAc, ~Fuc, ~Neu5Ac, ~Lys,
     "none/none",     0,       0,    0,     0,  0,
     "none/Lys",     0,       0,    0,     0,  1,
     "Lys/Lys",     0,       0,    0,     0,  2,
     "1xHex",     1,       0,    0,     0,  0,
     "2xHex",     2,       0,    0,     0,  0,
   ) %>% 
     define_ptm_compositions(c(Lys = "C6H12N2O"))
}

# specify paths ------------------------------------------------------------

df <- tibble(mzml_full_path = dir_ls(path = directory,regexp = pattern),) %>%
  separate(mzml_full_path,
           into = c("data",  "subdir1", "subdir2","filename"),
           sep = "/",
           remove = FALSE) %>%
  mutate(analysis_path = fs::path("analysis",
                                  subdir1, 
                                  "pngase_output_tables",
                                  gsub("\\..*$", "", filename))) 
  
fs::dir_create(df$analysis_path)

# wrangle out info ---------------------------------------------------

df <- separate_wider_delim(df,
                           filename,
                           delim = "_",
                           # names = c("ymd", "initials", "CHO_cell_variant", "bio_replicate","enzyme2","tech_replicate","aquisition_number"),
                           names = c("ymd", "initials", "CHO_cell_variant", "bio_replicate","cpb", "enzyme2","aquisition_number"),
                           # too_few = "debug"
                           ) %>%
  unite(CHO_cell_variant_bio_replicate, 
        c("CHO_cell_variant","bio_replicate"), 
        remove = FALSE) %>%
  mutate(tech_replicate = rep(c(1,2,3), 14)) %>% #perhaps is redundant
  # mutate(tech_replicate = rep(c(1,2,3), 12)) %>% #perhaps is redundant
  {.}

# write_csv(x = df,
#           file = "data/Dec_2023/overview_sample_pngase.csv")
write_csv(x = df,
          file = "data/3_lysine_glycation_quantification/overview_sample_pngase_lysine.csv")

# import data with information on charge states and rt limits  --------

# cs_rt_data <- read_csv('data/Dec_2023/rt_seconds_pngase.csv') %>%
cs_rt_data <- read_csv('data/3_lysine_glycation_quantification/rt_seconds_pngase_lysine.csv') %>%
  separate(sample_name,
           into = c("ymd", "initials", "CHO_cell_variant", "bio_replicate2","enzyme","aquisition_number"),
           sep = "_",
           remove = FALSE) %>%
  filter(subunit == "intact")
  # %>%
  # unite(CHO_cell_variant_bio_replicate, 
  #       c("CHO_cell_variant","bio_replicate2"), 
  #       remove = FALSE)


# merge data sample with data cs and rt -----------------------------------
data_merged <- df %>% 
  select(mzml_full_path, CHO_cell_variant, analysis_path, tech_replicate, bio_replicate, CHO_cell_variant_bio_replicate) %>%
  left_join(cs_rt_data, by = "CHO_cell_variant") %>%
  #cs_rt_data %>%  select(CHO_cell_variant_bio_replicate, )
  filter(CHO_cell_variant != "A25") %>%# nearly no signal
  {.}

# write_csv(x = data_merged,
#           file = "data/Dec_2023/overview_sample_pngase_merged.csv")

write_csv(x = data_merged,
          file = "data/3_lysine_glycation_quantification/overview_sample_pngase_lysine_merged.csv")

# fragquaxi analysis ------------------------------------------------------


## custom function ---------------------------------------------------------
calculate_abundance <- function(mzml_full_path,
                                charge_state_50_min, 
                                charge_state_50_max, 
                                rt_start_sec,
                                rt_end_sec,
                                scan_number_start,
                                scan_number_end,
                                analysis_path, 
                                ...){
  
  ms_data <- mzR::openMSfile(mzml_full_path)
  print(ms_data)
  print(c(rt_start_sec,rt_end_sec))
  
  pfm_ions <-
    assemble_proteoforms(proteins, modcoms) %>%
    ionize(charge_states = c(charge_state_50_min:charge_state_50_max), ppm = 300)
  print(dim(pfm_ions)) #check that for every file the correct # of charge states was used
  
  extracted_filename <- str_extract(ms_data@fileName, "(?<=data\\/).*(?=\\.mzML)")
  # single charge states
  plot_ions(
    ms_data,
    ions = pfm_ions,
    scans = scan_number_start:scan_number_end,
    xlim = c(3270, 3330)
  )
  
  ggsave(filename = paste0("figures/",extracted_filename,"one_cs.png"),    
         height = 200,
         width = 300,
         units = "mm",
         dpi = 600)

  abundances <- quantify_ions(ms_data,
                              ions = pfm_ions,
                              rt_limits = c(rt_start_sec,rt_end_sec)
                              ) %>%
    as_tibble() %>% 
    mutate(modcom_name = factor(modcom_name) %>% 
    fct_inorder()
    ) %>% 
    unnest(abundance_data) %>% 
    group_by(modcom_name) %>%
    summarise(abundance = sum(abundance)) %>% 
    mutate(frac_ab = abundance / sum(abundance) * 100,
           file_name = mzml_full_path)

    write_csv(x = abundances,
              file = paste(analysis_path,"frac_ab_tb_cs50.csv",sep = "/")
              )
  
  print('Analysis finished')
}


## apply custom function to dfr --------------------------------------------

pwalk(data_merged, calculate_abundance, .progress = TRUE)





          
