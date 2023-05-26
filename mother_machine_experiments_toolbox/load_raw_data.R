# Load necessary packages
source("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/mother_machine_experiments_toolbox/load_functions_and_packages.R")

# Path to denoised data, gfp and q already corrected for autofluorescence
raw_data_folder <- "/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/raw_curated_data_complete_cycles"

# Setting conditions, promoters and datatypes
conditions <- c("acetate005","glycerol040","glucose020","glucoseaa020")
conditions_labels <- c("acetate","glycerol","glucose","glucose+a.a")
promoters <- c("hi1","hi3","med2","med3","rrnB","rplN","6300")
datatypes <- c("allnoisesources","nodivisionnoise","nogrowthnoise","noprodnoise","nonoise")
datatypes_labels <- c("all noise sources","no division noise","no growth noise","no production noise","no noise")

condition_label_df <- tibble(condition=conditions,condition_label=conditions_labels)
datatype_label_df <- tibble(datatype=datatypes,datatype_label=datatypes_labels)

# Importing raw experimental data
raw_data <- lapply(raw_data_folder,function(.l){
  list.files(path=.l,pattern="rawdata.csv$",recursive = TRUE,full.names=TRUE)}) %>% 
  unlist() %>% 
  #lapply(function(.l) readr::read_csv(.l) %>% 
  mclapply(function(.l) readr::read_csv(.l) %>%
             mutate(path=.l) %>% 
             mutate(datatype="rawdata")) %>% 
  do.call(rbind, .)