# We assume here, that raw data were curated, that denoising procedure was performed, in the following ouput folder
denoised_data_files_dir <- "/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoising_raw_data/denoising_20221215091650"

# The following script will import the denoised data together with the raw data to output a tidied table containing all the info about the datasets.
# Complete denoised data will be saved in ./denoised_data and available for further analysis.

data_folder <- c("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/raw_curated_data")

source("./mother_machine_experiments_toolbox/load_functions_and_packages.R")

list_of_data_files <- lapply(data_folder,function(.l){
  list.files(path=.l,pattern="rawdata.csv$",recursive = TRUE,full.names=TRUE)}) %>% 
  unlist()

df_of_data_files <- tibble(path=list_of_data_files) %>% 
  tidyr::extract(path,"condition","[0-9a-zA-Z]{1,}_([0-9a-zA-Z]{1,})_rawdata.csv$",remove=FALSE,convert=FALSE) %>% 
  tidyr::extract(path,"promoter","([0-9a-z]{1,})_[0-9a-zA-Z]{1,}_rawdata.csv$",remove=FALSE,convert=FALSE) %>% 
  rename(raw_curated_data_path=path) %>% 
  tidyr::extract(condition,"carbon_source_concentration","^[a-z]{1,}0([0-9]{2})",remove=FALSE,convert=FALSE) %>% 
  mutate(carbon_source_concentration=as.double(carbon_source_concentration)/100) %>% 
  tidyr::extract(condition,"carbon_source","^([a-z]{1,})0[0-9]{2}",remove=FALSE,convert=FALSE)

myframes_complete_raw <- df_of_data_files %>% 
  ungroup() %>% 
  distinct(raw_curated_data_path) %>% 
  .$raw_curated_data_path %>% 
  lapply(function(.l) readr::read_csv(.l) %>% 
           mutate(raw_curated_data_path=.l)) %>% 
  do.call(rbind, .) %>% 
  #rename(carbon_source=condition) %>% 
  left_join(df_of_data_files,by=c("condition","promoter","raw_curated_data_path"))


find_command <- sprintf("find %s -type f -name '*prediction.csv'",denoised_data_files_dir)
list_of_denoised_data_files <- system(find_command, intern=TRUE)

myframes_complete_denoised <- list_of_denoised_data_files %>% 
  lapply(function(.l) readr::read_csv(.l,skip=13) %>% 
           mutate(denoised_data_path=.l)) %>% 
  do.call(rbind,.) %>% 
  extract(denoised_data_path,"promotercondition","/([a-z0-9A-Z]{1,}_[0-9a-zA-Z]{1,})_rawdata_[a-z0-9]{1,}_b_prediction.csv$",remove=FALSE,convert=FALSE) %>% 
  separate(promotercondition,into=c("promoter","condition"),sep="_") %>% 
  separate(condition,into=c("carbon_source","carbon_source_concentration"),sep="(?<=[A-Za-z])(?=[0-9])",remove=FALSE) %>% 
  mutate(carbon_source_concentration=as.double(carbon_source_concentration)/100)

find_command <- sprintf("find %s -type f -name '*final.csv'",denoised_data_files_dir)

myparameters <- system(find_command, intern=TRUE) %>% 
  lapply(function(.l) readr::read_csv(.l,n_max=10) %>% 
           mutate(denoised_data_path=.l)) %>% 
  do.call(rbind,.) %>% 
  extract(denoised_data_path,"promotercondition","/([a-z0-9A-Z]{1,}_[0-9a-zA-Z]{1,})_rawdata_[a-z0-9]{1,}_b_final.csv$",remove=FALSE,convert=FALSE) %>% 
  separate(promotercondition,into=c("promoter","condition"),sep="_") %>% 
  separate(condition,into=c("carbon_source","carbon_source_concentration"),sep="(?<=[A-Za-z])(?=[0-9])",remove=FALSE) %>% 
  mutate(carbon_source_concentration=as.double(carbon_source_concentration)/100)

myerrors <- system(find_command, intern=TRUE) %>% 
  lapply(function(.l) readr::read_csv(.l,n_max=3,skip=14) %>% 
           mutate(denoised_data_path=.l)) %>% 
  do.call(rbind,.) %>% 
  extract(denoised_data_path,"promotercondition","/([a-z0-9A-Z]{1,}_[0-9a-zA-Z]{1,})_rawdata_[a-z0-9]{1,}_b_final.csv$",remove=FALSE,convert=FALSE) %>% 
  separate(promotercondition,into=c("promoter","condition"),sep="_") %>% 
  separate(condition,into=c("carbon_source","carbon_source_concentration"),sep="(?<=[A-Za-z])(?=[0-9])",remove=FALSE) %>% 
  mutate(carbon_source_concentration=as.double(carbon_source_concentration)/100)

myframes_complete_raw <- myframes_complete_raw %>% 
  select(-c(length_um,gfp_nb,carbon_source_concentration,carbon_source,parent)) %>% 
  ungroup()

myframes_complete_denoised <- myframes_complete_denoised %>% 
  ungroup() %>% 
  rename(cell=cell_id,
         parent=parent_id,
         time_min=time,
         gfp_nb=mean_g,
         lambda=mean_l,
         q=mean_q) %>% 
  mutate(time_sec=time_min*60,
         length_um_raw=exp(log_length),
         gfp_nb_raw=fp,
         length_um=exp(mean_x)) %>% 
  select(-c(log_length,fp,mean_x))

myframes_final <- left_join(myframes_complete_raw,myframes_complete_denoised,by=c("condition","promoter","cell","time_sec"))

myframes_final %>% 
  group_by(promoter,condition) %>% 
  do((function(.df) {
    readr::write_csv(.df,sprintf("./denoised_data/%s_%s_denoised.csv",unique(.df$promoter),unique(.df$condition)))
    return(tibble())
  })(.)) %>% 
  ungroup()