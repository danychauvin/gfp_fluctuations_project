noise_contribution_data <- lapply(simulation_data_folder,function(.l){
  list.files(path=.l,pattern=".csv$",recursive = TRUE,full.names=TRUE)}) %>% 
  unlist() %>% 
  lapply(function(.l) readr::read_csv(.l) %>% 
           mutate(path=.l)) %>% 
  do.call(rbind, .) %>% 
  tidyr::extract(path,"promoter","/([0-9a-zA-Z]{1,})_[0-9a-zA-Z]{1,}_[a-z]{1,}.csv$",remove=FALSE,convert=FALSE) %>% 
  tidyr::extract(path,"condition","/[0-9a-zA-Z]{1,}_([0-9a-zA-Z]{1,})_[a-z]{1,}.csv$",remove=FALSE,convert=FALSE) %>%
  tidyr::extract(condition,"carbon_source_concentration","^[a-z]{1,}0([0-9]{2})",remove=FALSE,convert=FALSE) %>% 
  mutate(carbon_source_concentration=as.double(carbon_source_concentration)/100) %>% 
  tidyr::extract(condition,"carbon_source","^([a-z]{1,})0[0-9]{2}",remove=FALSE,convert=FALSE) %>% 
  tidyr::extract(path,"datatype","/[0-9a-z]{1,}_[0-9a-zA-Z]{1,}_([a-z]{1,}).csv$",remove=FALSE,convert=FALSE) %>% 
  rename(cell=cell_id,parent=parent_id,gfp_nb=gfp,lambda=lt,q=qt) %>% 
  mutate(time_sec=time_min*60,
  length_um=exp(log_length)) %>%
  select(-c("...1"))

experimental_data <- lapply(experimental_data_folder,function(.l){
  list.files(path=.l,pattern="denoised[_complete]*.csv$",recursive = TRUE,full.names=TRUE)}) %>% 
  unlist() %>% 
  lapply(function(.l) readr::read_csv(.l) %>% 
           mutate(path=.l)) %>% 
  do.call(rbind, .) %>% 
  mutate(datatype="experimentaldata") %>% 
  mutate(time_sec=time_min*60) %>% 
  mutate(log_length=log(length_um)) %>% 
  select(c("cell","time_min","parent","log_length","gfp_nb","lambda","q","promoter","condition","carbon_source","carbon_source_concentration","time_sec","length_um","datatype","path","mean_cell_rank","gl_id"))

mean_cell_rank_df <- experimental_data %>% 
  distinct(cell,mean_cell_rank,gl_id)

noise_contribution_data <- rbind(noise_contribution_data,
                                 experimental_data %>% 
                                   select(-c(mean_cell_rank,gl_id))) %>% 
  left_join(mean_cell_rank_df,by=c("cell")) %>% 
  mutate(concentration=gfp_nb/length_um)
