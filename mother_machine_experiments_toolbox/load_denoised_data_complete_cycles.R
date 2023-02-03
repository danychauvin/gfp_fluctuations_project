# We assume here, that raw data were curated, that denoising procedure was performed, and denoised data completed with raw data and saved using import_denoised_data.R
denoised_data_dir <- "/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoised_data_complete_cycles"
# The following script will load all denoised data in the specified folder.

myframes <- lapply(denoised_data_dir,function(.l){
  list.files(path=.l,pattern="denoised_complete.csv$",recursive = TRUE,full.names=TRUE)}) %>% 
  unlist() %>% 
  lapply(function(.l) readr::read_csv(.l) %>% 
           mutate(denoised_data_path=.l)) %>% 
  do.call(rbind, .)