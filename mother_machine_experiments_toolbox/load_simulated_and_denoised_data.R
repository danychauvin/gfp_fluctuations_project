# Load necessary functions
source("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/mother_machine_experiments_toolbox/load_functions_and_packages.R")

# Path to simulated data, where inferred mean q was already corrected for autofluorescence
simulation_data_folder <- "/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoised_data_complete_cycles_corrected/integration"

# Path to denoised data, gfp and q already corrected for autofluorescence
experimental_data_folder <- "/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoised_data_complete_cycles"

# Path to correlation data
correlation_data_folder <- "/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoising_raw_data/denoising_20230512104746"

# Setting conditions, promoters and datatypes
conditions <- c("acetate005","glycerol040","glucose020","glucoseaa020")
conditions_labels <- c("acetate","glycerol","glucose","glucose+a.a")
promoters <- c("hi1","hi3","med2","med3","rrnB","rplN")
datatypes <- c("allnoisesources","nodivisionnoise","nogrowthnoise","noprodnoise","nonoise")
datatypes_labels <- c("all noise sources","no division noise","no growth noise","no production noise","no noise")

condition_label_df <- tibble(condition=conditions,condition_label=conditions_labels)
datatype_label_df <- tibble(datatype=datatypes,datatype_label=datatypes_labels)

# Importing experimental data (denoised)
experimental_data <- lapply(experimental_data_folder,function(.l){
  list.files(path=.l,pattern="denoised[_complete]*.csv$",recursive = TRUE,full.names=TRUE)}) %>% 
  unlist() %>% 
  #lapply(function(.l) readr::read_csv(.l) %>% 
  mclapply(function(.l) readr::read_csv(.l) %>%
             mutate(path=.l) %>% 
             mutate(datatype="experimentaldata") %>% 
             mutate(time_sec=time_min*60) %>% 
             mutate(log_length=log(length_um)) %>% 
             rename(gfp_nb_uncorrected=gfp_nb) %>% 
             rename(gfp_nb=gfp_nb_corrected) %>% 
             rename(q_uncorrected=q) %>% 
             rename(q=q_corrected)) %>% 
  do.call(rbind, .) %>% 
  mutate(concentration=gfp_nb/length_um,
         l_q=log(q),
         l_concentration=log(concentration)) %>% 
  mutate(sd_concentration=sqrt((cov_gg+gfp_nb**2*cov_xx-2*gfp_nb*cov_xg)/length_um**2),
         sd_q=sqrt(abs(cov_qq)),
         sd_lambda=sqrt(abs(cov_ll)),
         sd_length_um=length_um*sqrt(abs(cov_xx)),
         sd_gfp_nb=sqrt(abs(cov_gg)),
         sd_l_q=sqrt(abs(1/q**2*cov_qq)),
         sd_l_c=sqrt(abs(cov_gg/(gfp_nb**2)+cov_xx-(2/gfp_nb)*cov_xg)))

# Importing noise contribution data (forward integration), and binding with experimental data
noise_contribution_data <- lapply(simulation_data_folder,function(.l){
  list.files(path=.l,pattern=".csv$",recursive = TRUE,full.names=TRUE)}) %>% 
  unlist() %>% 
  #lapply(prepare_simulation_csv) %>% 
  mclapply(prepare_simulation_csv,mc.cores=(detectCores()-1)) %>% 
  do.call(rbind, .) %>% 
  select(-c("...1"))

mean_cell_rank_df <- experimental_data %>% 
  distinct(cell,mean_cell_rank,gl_id,replica)

noise_contribution_data <- rbind(noise_contribution_data,
                                 experimental_data %>% 
                                   select(c("date","cell","time_min","parent","log_length","gfp_nb","lambda","q","promoter","condition","carbon_source","carbon_source_concentration","time_sec","length_um","datatype","path"))) %>% 
  left_join(mean_cell_rank_df,by=c("cell")) %>% 
  mutate(concentration=gfp_nb/length_um)


# Compute mean doubling times (in hour) for each dataset  
mycluster <- load_cluster()
copy_cluster(mycluster,c("compute_weighted_mean"))

doubling_times_df <- experimental_data %>% 
  filter(condition %in% conditions) %>% 
  filter(promoter %in% promoters) %>% 
  mutate(dataset=paste(condition,promoter,replica,date,sep=".")) %>% 
  group_by(dataset) %>% 
  partition(mycluster) %>% 
  mutate(mean_l=compute_weighted_mean(lambda,sqrt(abs(cov_ll))),
         mean_doubling_time=1/(mean_l*60/log(2))) %>% 
  collect() %>% 
  ungroup() %>% 
  distinct(dataset,mean_doubling_time) %>% 
  mutate(time_threshold=3*mean_doubling_time) %>% 
  separate(dataset,into=c("condition","promoter","replica","date"),remove=FALSE)

# Filter noise contribution data to keep only cells that are characterized by at least 3 generations. They should have forgotten about the initial concentration in their lineage.
old_cells <- experimental_data %>% 
  distinct(cell,parent) %>% 
  rename(mother=parent)

mothers <- experimental_data %>% 
  distinct(cell,parent) %>% 
  rename(g_mother=parent,
         mother=cell)

g_mothers <- experimental_data %>% 
  distinct(cell,parent) %>% 
  rename(gg_mother=parent,
         g_mother=cell)

gg_mothers <- experimental_data %>% 
  distinct(cell,parent) %>% 
  rename(ggg_mother=parent,
         gg_mother=cell)

# t <- old_cells_filtered <- old_cells %>% 
#   left_join(mothers,by=c("mother")) %>% 
#   left_join(g_mothers,by=c("g_mother")) %>% 
#   left_join(gg_mothers,by=c("gg_mother")) #%>% 
  #na.omit() %>% 
  #distinct(cell)

old_cells_filtered <- old_cells %>% 
  left_join(mothers,by=c("mother")) %>%
  left_join(g_mothers,by=c("g_mother")) %>%
  left_join(gg_mothers,by=c("gg_mother")) %>%
  na.omit() %>%
  distinct(cell)

noise_contribution_data_filtered <- noise_contribution_data %>%
  filter(promoter %in% promoters) %>% 
  filter(condition %in% conditions) %>% 
  mutate(dataset=paste(condition,promoter,replica,date,sep=".")) %>% 
  mutate(replica=as.character(replica)) %>% 
  group_by(dataset) %>% 
  partition(mycluster) %>% 
  mutate(time_sec_min=min(time_sec,na.rm=TRUE)) %>%
  mutate(time_h_rescaled=time_sec/3600-time_sec_min/3600) %>% 
  collect() %>% 
  ungroup() %>% 
  left_join(doubling_times_df,by=c("dataset","condition","promoter","date","replica")) %>% 
  filter(time_h_rescaled>time_threshold) %>% 
  semi_join(old_cells_filtered,by=c("cell"))  

# Load correlation data, parameters and errors from denoising
list_of_corr_files <- lapply(correlation_data_folder,function(.l){
  list.files(path=.l,pattern="correlations.csv$",recursive = TRUE,full.names=TRUE)}) %>% 
  unlist()

df_of_data_files <- tibble(path=list_of_corr_files) %>% 
  tidyr::extract(path,"condition","([0-9a-z]{1,})_[0-9a-zA-Z]{1,}_[0-9]{8}_rawdata_[a-z0-9]{1,}_b_correlations.csv$",remove=FALSE,convert=FALSE) %>% 
  tidyr::extract(path,"promoter","[0-9a-z]{1,}_([0-9a-zA-Z]{1,})_[0-9]{8}_rawdata_[a-z0-9]{1,}_b_correlations.csv$",remove=FALSE,convert=FALSE) %>% 
  tidyr::extract(path,"date","[0-9a-z]{1,}_[0-9a-zA-Z]{1,}_([0-9]{8})_rawdata_[a-z0-9]{1,}_b_correlations.csv$",remove=FALSE,convert=FALSE) %>% 
  rename(raw_curated_data_path=path) %>% 
  tidyr::extract(condition,"carbon_source_concentration","^[a-z]{1,}0([0-9]{2})",remove=FALSE,convert=FALSE) %>% 
  mutate(carbon_source_concentration=as.double(carbon_source_concentration)/100) %>% 
  tidyr::extract(condition,"carbon_source","^([a-z]{1,})0[0-9]{2}",remove=FALSE,convert=FALSE)

correlations_df <- df_of_data_files %>% 
  ungroup() %>% 
  distinct(raw_curated_data_path) %>% 
  .$raw_curated_data_path %>% 
  lapply(function(.l) readr::read_csv(.l) %>% 
           mutate(raw_curated_data_path=.l)) %>% 
  do.call(rbind, .) %>% 
  left_join(df_of_data_files,by=c("raw_curated_data_path"))

correlations_df <- correlations_df %>% 
  rename(cov_l_l="cov_l(t+dt)l(t)",
         cov_l_l_err="cov_l(t+dt)l(t)_err",
         col_q_q="cov_q(t+dt)q(t)",
         cov_q_q_err="cov_q(t+dt)q(t)_err",
         cov_q_l="cov_q(t+dt)l(t)",
         cov_q_l_err="cov_q(t+dt)l(t)_err",
         cov_l_q="cov_l(t+dt)q(t)",
         cov_l_q_err="cov_l(t+dt)q(t)_err",
         cov_c_c="cov_c(t+dt)c(t)",
         cov_c_c_err="cov_c(t+dt)c(t)_err",
         corr_l_l="corr_l(t+dt)l(t)",
         corr_l_l_err="corr_l(t+dt)l(t)_err",
         corr_q_q="corr_q(t+dt)q(t)",
         corr_q_q_err="corr_q(t+dt)q(t)_err",
         corr_q_l="corr_q(t+dt)l(t)",
         corr_q_l_err="corr_q(t+dt)l(t)_err",
         corr_l_q="corr_l(t+dt)q(t)",
         corr_l_q_err="corr_l(t+dt)q(t)_err",
         corr_c_c="corr_c(t+dt)c(t)",
         corr_c_c_err="corr_c(t+dt)c(t)_err",
         corr_naive_l_l="corr_naive_l(t+dt)l(t)",
         corr_naive_l_q="corr_naive_l(t+dt)q(t)",
         corr_naive_q_l="corr_naive_q(t+dt)l(t)",
         corr_naive_q_q="corr_naive_q(t+dt)q(t)",
         corr_naive_c_c="corr_naive_c(t+dt)c(t)")

#Import parameters and errors, inferred from denoising
find_command <- sprintf("find %s -type f -name '*final.csv'",correlation_data_folder)

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
