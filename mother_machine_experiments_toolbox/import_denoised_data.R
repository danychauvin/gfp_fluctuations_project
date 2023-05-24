#Define where raw data are 
data_folder <- c("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/raw_curated_data_complete_cycles")

#Define denoising output folder
denoised_data_files_dir <- "/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoising_raw_data/denoising_20230512104746"

#Import raw data
source("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/mother_machine_experiments_toolbox/load_functions_and_packages.R")
list_of_data_files <- lapply(data_folder,function(.l){
  list.files(path=.l,pattern="rawdata.csv$",recursive = TRUE,full.names=TRUE)}) %>% 
  unlist()

df_of_data_files <- tibble(path=list_of_data_files) %>% 
    tidyr::extract(path,"condition","([0-9a-zA-Z]{1,})_[0-9a-zA-Z]{1,}_[0-9]{8}_rawdata.csv$",remove=FALSE,convert=FALSE) %>% 
    tidyr::extract(path,"promoter","[0-9a-z]{1,}_([0-9a-zA-Z]{1,})_[0-9]{8}_rawdata.csv$",remove=FALSE,convert=FALSE) %>% 
    tidyr::extract(path,"date","[0-9a-z]{1,}_[0-9a-zA-Z]{1,}_([0-9]{8})_rawdata.csv$",remove=FALSE,convert=FALSE) %>% 
  rename(raw_curated_data_path=path) %>% 
  tidyr::extract(condition,"carbon_source_concentration","^[a-z]{1,}0([0-9]{2})",remove=FALSE,convert=FALSE) %>% 
  mutate(carbon_source_concentration=as.double(carbon_source_concentration)/100) %>% 
  tidyr::extract(condition,"carbon_source","^([a-z]{1,})0[0-9]{2}",remove=FALSE,convert=FALSE)

myframes_raw <- df_of_data_files %>% 
  ungroup() %>% 
  distinct(raw_curated_data_path) %>% 
  .$raw_curated_data_path %>% 
  lapply(function(.l) readr::read_csv(.l) %>% 
  mutate(raw_curated_data_path=.l)) %>% 
  do.call(rbind, .) %>% 
  mutate(date=as.character(date)) %>% 
  #rename(carbon_source=condition) %>% 
  left_join(df_of_data_files,by=c("condition","promoter","date","raw_curated_data_path"))

# Import denoised data (as well as parameters and errors)
find_command <- sprintf("find %s -type f -name '*prediction.csv'",denoised_data_files_dir)
list_of_denoised_data_files <- system(find_command, intern=TRUE)

myframes_denoised <- list_of_denoised_data_files %>% 
  lapply(function(.l) readr::read_csv(.l,skip=13) %>% 
           mutate(denoised_data_path=.l)) %>% 
  do.call(rbind,.) %>% 
  extract(denoised_data_path,"conditionpromoter","/([a-z0-9A-Z]{1,}_[0-9a-zA-Z]{1,})_[0-9]{8}_rawdata_[a-z0-9]{1,}_b_prediction.csv$",remove=FALSE,convert=FALSE) %>% 
  extract(denoised_data_path,"date","/[a-z0-9A-Z]{1,}_[0-9a-zA-Z]{1,}_([0-9]{8})_rawdata_[a-z0-9]{1,}_b_prediction.csv$",remove=FALSE,convert=FALSE) %>% 
  separate(conditionpromoter,into=c("condition","promoter"),sep="_") %>% 
  separate(condition,into=c("carbon_source","carbon_source_concentration"),sep="(?<=[A-Za-z])(?=[0-9])",remove=FALSE) %>% 
  mutate(carbon_source_concentration=as.double(carbon_source_concentration)/100)
  
find_command <- sprintf("find %s -type f -name '*final.csv'",denoised_data_files_dir)
myparameters <- system(find_command, intern=TRUE) %>% 
  lapply(function(.l) readr::read_csv(.l,n_max=11) %>% 
           mutate(denoised_data_path=.l)) %>% 
  do.call(rbind,.) %>% 
  extract(denoised_data_path,"conditionpromoter","/([a-z0-9A-Z]{1,}_[0-9a-zA-Z]{1,})_[0-9]{8}_rawdata_[a-z0-9]{1,}_b_final.csv$",remove=FALSE,convert=FALSE) %>% 
  extract(denoised_data_path,"date","/[a-z0-9A-Z]{1,}_[0-9a-zA-Z]{1,}_([0-9]{8})_rawdata_[a-z0-9]{1,}_b_final.csv$",remove=FALSE,convert=FALSE) %>% 
  separate(conditionpromoter,into=c("condition","promoter"),sep="_") %>% 
  separate(condition,into=c("carbon_source","carbon_source_concentration"),sep="(?<=[A-Za-z])(?=[0-9])",remove=FALSE) %>% 
  mutate(carbon_source_concentration=as.double(carbon_source_concentration)/100)

myerrors <- system(find_command, intern=TRUE) %>% 
  lapply(function(.l) readr::read_csv(.l,n_max=3,skip=14) %>% 
           mutate(denoised_data_path=.l)) %>% 
  do.call(rbind,.) %>% 
 extract(denoised_data_path,"conditionpromoter","/([a-z0-9A-Z]{1,}_[0-9a-zA-Z]{1,})_[0-9]{8}_rawdata_[a-z0-9]{1,}_b_final.csv$",remove=FALSE,convert=FALSE) %>% 
  extract(denoised_data_path,"date","/[a-z0-9A-Z]{1,}_[0-9a-zA-Z]{1,}_([0-9]{8})_rawdata_[a-z0-9]{1,}_b_final.csv$",remove=FALSE,convert=FALSE) %>% 
  separate(conditionpromoter,into=c("condition","promoter"),sep="_") %>% 
  separate(condition,into=c("carbon_source","carbon_source_concentration"),sep="(?<=[A-Za-z])(?=[0-9])",remove=FALSE) %>% 
  mutate(carbon_source_concentration=as.double(carbon_source_concentration)/100)

# Concatenate denoised data and raw data
myframes_raw <- myframes_raw %>% 
  select(-c(length_um,gfp_nb,carbon_source_concentration,carbon_source,parent)) %>% 
  ungroup()

myframes_denoised <- myframes_denoised %>% 
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

myframes_final <- left_join(myframes_raw,myframes_denoised,by=c("condition","promoter","date","cell","time_sec"))

# Save the data on the disk for further analysis

## Filtering out cells which have bad segmentation issue.
bad_cells <- myframes_final %>%
  mutate(lambda_high_boundary=lambda+sqrt(cov_ll)) %>% 
  #filter(lambda_high_boundary<0) %>% 
  distinct(cell)

myframes_final <- myframes_final %>% 
  anti_join(bad_cells,by=c("cell")) %>% 
  mutate(date=as.character(date))

# importing auto-fluorescence correction table and correcting gfp_nb and q for auto-fluorescence
autofluorescence_table <- readr::read_csv("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/autofluorescence_correction/autofluorescence_table.csv") %>% 
  mutate(date=as.character(date),
         control_experiment_date=as.character(control_experiment_date))

# importing auto-fluorescence correction table and correcting gfp_nb and q for auto-fluorescence
autofluo_df <- myframes_final %>% 
  filter(promoter=="6300") %>% 
  mutate(concentration=gfp_nb/length_um) %>% 
  group_by(date,condition,promoter,cell) %>%
  sample_n(1) %>% #one random observation by cell
  ungroup() %>% 
  group_by(date,condition,promoter) %>%
  mutate(sd_g_bar=sd(concentration),
         g_bar=mean(concentration,na.rm=TRUE),
         q_bar=mean(q,na.rm=TRUE),
         sd_q_bar=sd(q,na.rm=TRUE)) %>% 
  ungroup() %>% 
  rename("control_experiment_date"=date) %>% 
  distinct(condition,control_experiment_date,g_bar,sd_g_bar,q_bar,sd_q_bar)

myframes_final <- myframes_final %>% 
  left_join(autofluorescence_table,by=c("date")) %>% 
  left_join(autofluo_df, by=c("condition","control_experiment_date")) %>% 
  mutate(sd_g_autofluo=sd_g_bar*length_um+exp(sqrt(cov_xx))*g_bar,
         gfp_nb_corrected=gfp_nb-g_bar*length_um,
         q_corrected=q-q_bar,
         sd_q_autofluo=sd_q_bar)

# Save the data on the disk for further analysis
myframes_final %>% 
  group_by(date,promoter,condition) %>% 
  do((function(.df) {
    #readr::write_csv(.df,sprintf("../denoised_data/%s_%s_denoised.csv",unique(.df$promoter),unique(.df$condition)))
    readr::write_csv(.df,sprintf("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoised_data_complete_cycles/%s_%s_%s_date_denoised_complete.csv",unique(.df$condition),unique(.df$promoter),unique(.df$date)))
    return(tibble())
  })(.)) %>% 
  ungroup()

# Save the data on the disk for forward integration, with mean_q corrected for autofluorescence molecules correction

myparameters_corrected <- myparameters %>%
  left_join(autofluo_df %>% select(condition,q_bar),by=c("condition")) %>%
  mutate(final=ifelse(name=="mean_q",final-q_bar,final)) %>%
  select(c(no,name,type,init,step,lower_bound,upper_bound,final,condition,promoter,date)) %>%
  mutate(lower_bound=" ",
         upper_bound=" ") %>%
  mutate(step=ifelse(is.na(step)," ",as.character(step)))

myframes_final %>%
  filter(promoter!="6300") %>%
  select(-c(parent_id,mean_q)) %>%
  rename(mean_g=gfp_nb_corrected,
         mean_l=lambda,
         mean_q=q_corrected,
         parent_id=parent,
         cell_id=cell) %>%
  mutate(fp=mean_g) %>%
  mutate(mean_x=log(length_um),
         time=time_sec/60,
         log_length=mean_x) %>%
  select(c(date,condition,promoter,cell_id,parent_id,time,log_length,fp,mean_x,mean_g,mean_l,mean_q,cov_xx,cov_xg,cov_xl,cov_xq,cov_gg,cov_gl,cov_gq,cov_ll,cov_lq,cov_qq)) %>%
  group_by(date,promoter,condition) %>%
  do((function(.df) {
    .cond <- unique(.df$condition)
    .prom <- unique(.df$promoter)
    .date <- unique(.df$date)
    #readr::write_csv(.df,sprintf("../denoised_data/%s_%s_denoised.csv",unique(.df$promoter),unique(.df$condition)))
    readr::write_csv(myparameters_corrected %>% filter(condition==.cond,promoter==.prom,date==.date) %>% select(no,name,type,init,step,lower_bound,upper_bound,final),sprintf("../denoised_data_complete_cycles_corrected/%s_%s_%s_date_denoised_complete.csv",unique(.df$condition),unique(.df$promoter),unique(.df$date)))
    readr::write_file("\n",sprintf("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoised_data_complete_cycles_corrected/%s_%s_%s_date_denoised_complete.csv",unique(.df$condition),unique(.df$promoter),unique(.df$date)),append = TRUE)
    readr::write_csv(.df %>% select(-c(date,condition,promoter)),sprintf("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoised_data_complete_cycles_corrected/%s_%s_%s_date_denoised_complete.csv",unique(.df$condition),unique(.df$promoter),unique(.df$date)),append=TRUE,col_names = TRUE)
    return(tibble())
  })(.)) %>%
  ungroup()