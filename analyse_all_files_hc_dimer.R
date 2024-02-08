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
product <- "hc_dimer" #OR product <- "pngase" OR product <- "sample" 

directory <- fs::dir_create("data/Jan_2024")

# depending on product, define regexp pattern and modcoms -----------------

if (product == "sample") {
  #specify pattern to match regexp for filename
  #pattern <-  "ambr.*\\.mzML"
  pattern <-  "CpB.*\\.mzML"
  
  # define mab_sequence  --------------------------------
  
  mab_sequence <- ("mab_sequence/trastuzumab_noLysine.txt")
  
  proteins <- define_proteins(
    mab = mab_sequence,
    .disulfides = 16
  )
  
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
    "G1F/S1G1F",      9,       8,    2,     1,     
  ) %>% 
    define_ptm_compositions()
   
 } else if (product == "pngase")  {
   #specify pattern to match regexp for filename 
   pattern <- "PNGaseF.*\\.mzML"
   
   # define mab_sequence  --------------------------------
   
   mab_sequence <- ("mab_sequence/trastuzumab_noLysine.txt")
   
   proteins <- define_proteins(
     mab = mab_sequence,
     .disulfides = 16
   )
   #specify modifications composition
   modcoms <- tribble(
     ~modcom_name, ~Hex, ~HexNAc, ~Fuc, ~Neu5Ac,
     "none", 0, 0, 0, 0,
     "1xHex", 1, 0, 0, 0,
     "2xHex", 2, 0, 0, 0,
     "3xHex", 3, 0, 0, 0, 
   ) %>%
     define_ptm_compositions()
   
 } else if (product == "hc_dimer")  {
   #specify pattern to match regexp for filename 
   pattern <- ".*(A19|A16).*CpB.*\\.mzML"
   
   # define mab_sequence  --------------------------------
   
   mab_sequence <- ("mab_sequence/trastuzumab_two_hc_noLysine.txt")
   
   proteins <- define_proteins(
     hc_dimer = mab_sequence,
     .disulfides = 10
   )
   
   #obtained from Kathi's script Code_trastuzumab.R
   #modified to fit with two heavy chain, doubly glycosylated
   modcoms <- tribble(
     ~modcom_name, ~Hex, ~HexNAc, ~Fuc, ~Neu5Ac,
     "G0F/G0F",       6,       8,    2,     0,
     "G0F/G1F",       7,       8,    2,     0, 
     "G1F/G1F",       8,       8,    2,     0,
     "G1F/G2F",       9,       8,    2,     0,
     "G2F/G2F",      10,       8,    2,     0,
     "G0F/G0",        6,       8,    1,     0,
   ) %>%
     define_ptm_compositions()
 }

# specify paths ------------------------------------------------------------

df <- tibble(mzml_full_path = dir_ls(path = directory,regexp = pattern),) %>%
  separate(mzml_full_path,
           into = c("data", "exp_month","filename"),
           sep = "/",
           remove = FALSE) %>%
  mutate(analysis_path = fs::path("analysis",
                                  exp_month,
                                  gsub("\\..*$", "", filename),
                                  product) #added for specifically hc_dimer
         ) 
  
fs::dir_create(df$analysis_path)

# wrangle out info ---------------------------------------------------

df <- separate_wider_delim(df,
                           filename,
                           delim = "_",
                           names = c("ymd", "initials", "CHO_cell_variant", "bio_replicate","enzyme","aquisition_number"),
                           too_few = "debug") %>%
  unite(CHO_cell_variant_bio_replicate, 
        c("CHO_cell_variant","bio_replicate"), 
        remove = FALSE) %>%
  mutate(tech_replicate = rep(c(1,2,3), 4)) #perhaps is redundant

write_csv(x = df,
          file = "data/Jan_2024/overview_hc_dimer.csv")


# import data with information on charge states and rt limits  --------

cs_rt_data <- read_csv('data/Jan_2024/rt_seconds_Jan2024.csv') %>%
  separate(sample_name,
           into = c("ymd", "initials", "CHO_cell_variant", "bio_replicate","enzyme","aquisition_number"),
           sep = "_",
           remove = FALSE) %>%
  filter(subunit == "hc") %>%
  unite(CHO_cell_variant_bio_replicate, 
        c("CHO_cell_variant","bio_replicate"), 
        remove = FALSE)



# merge data sample with data cs and rt -----------------------------------
data_merged <- df %>% 
  left_join(cs_rt_data, by = "CHO_cell_variant_bio_replicate") %>%
  #cs_rt_data %>%  select(CHO_cell_variant_bio_replicate, )
  filter(CHO_cell_variant.x != "A25") # nearly no signal

write_csv(x = data_merged,
          file = "data/Jan_2024/overview_hc_dimer_merged.csv")

# fragquaxi analysis ------------------------------------------------------


## custom function ---------------------------------------------------------
calculate_abundance <- function(mzml_full_path,
                                charge_state_50_min, 
                                charge_state_50_max, 
                                rt_start_sec,
                                rt_end_sec, 
                                analysis_path, 
                                ...){
  
  ms_data <- mzR::openMSfile(mzml_full_path)
  print(ms_data)
  print(c(rt_start_sec,rt_end_sec))
  
  pfm_ions <-
    assemble_proteoforms(proteins, modcoms) %>%
    ionize(charge_states = c(charge_state_50_min:charge_state_50_max), ppm = 300)
  print(dim(pfm_ions)) #check that for every file the correct # of charge states was used
  
  # #need to filter out overlapping charge states
  # pfm_ions <- pfm_ions %>% 
  #   filter((modcom_name == "none/none" & !(z %in% c(53, 52, 50, 48, 47, 45, 43))) | 
  #            (modcom_name == "G0F/G0" & !(z %in% c(53))) |
  #            (modcom_name == "G0F/G0F" & !(z %in% c(51))) |
  #            (modcom_name == "G0F/G1F" & !(z %in% c(49, 48))) |
  #            (modcom_name == "G1F/G1F" & !(z %in% c(46))) |
  #            (modcom_name == "G1F/G2F" & !(z %in% c(44))) |
  #            (modcom_name %in% c("G2F/G2F", "none/G0", "none/G0F", "none/G1F", "none/G2F", "G2F/G2F", "G1F/S1G1F"))
  #   )
  
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
  
  write.table(x = abundances,
              file = paste(analysis_path,"frac_ab_tb_cs50.csv",sep = "/"),
              sep = "\t",
              row.names = TRUE,
              col.names = NA)
  
  print('Analysis finished')
}


## apply custom function to dfr --------------------------------------------

pwalk(data_merged, calculate_abundance, .progress = TRUE)











          
