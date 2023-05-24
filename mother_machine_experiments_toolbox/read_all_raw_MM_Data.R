source("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/mother_machine_experiments_toolbox/load_functions_and_packages.R")
path_to_MM_data_summary <- "/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/mother_machine_experiments_toolbox/all_mmexperiments_datalist.csv"
source("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/mother_machine_experiments_toolbox/import_tidy_filter_raw_MMData.R")

myframes_complete_goodcells <- myframes_complete %>% 
  group_by(date,condition,promoter,cell) %>% 
  arrange(frame) %>% 
  filter(n()>4) %>% 
  mutate(l_birth=first(length_um_oriented_bbox)) %>% 
  mutate(l_div=last(length_um_oriented_bbox)) %>% 
  mutate(added_length=l_div-l_birth) %>% 
  mutate(g_birth=first(gfp_nb)) %>% 
  mutate(g_div=last(gfp_nb)) %>% 
  mutate(added_gfp=g_div-g_birth) %>% 
  mutate(mean_q=added_gfp/added_length*alpha_oriented_bbox) %>% 
  mutate(mean_c=mean(gfp_nb/length_um_oriented_bbox,na.rm=TRUE)) %>% 
  ungroup() %>% 
  arrange(date,condition,promoter,cell,time_sec)

myframes_complete_goodcells$mean_cell_rank <- factor(myframes_complete_goodcells$mean_cell_rank,levels=as.character(sort(unique(as.integer(myframes_complete_goodcells$mean_cell_rank)))))

myframes_complete_goodcells %>% 
  #anti_join(wrong_gl_df,by=c("gl_id")) %>% 
  group_by(date,promoter,condition) %>% 
  do((function(.df){
    .promoter <- unique(.df$promoter)
    .condition <- unique(.df$condition) 
    .date <- unique(.df$date)
    readr::write_csv(.df,sprintf("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/raw_curated_data_complete_cycles/%s_%s_%s_rawdata.csv",.condition,.promoter,.date))})(.)) %>% 
  ungroup()
