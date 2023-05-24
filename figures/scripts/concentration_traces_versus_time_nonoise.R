#Compare number of cells in noise_contribution data before and after filtering

noise_contribution_data_filtered %>% 
  filter(datatype=="nonoise") %>% 
  ggplot(aes(time_sec,concentration))+
  geom_line(aes(group=cell))+
  facet_wrap(~dataset,ncol=1,scale="free")

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/concentration_traces_versus_time_nonoise.pdf",width=5,height=50,limitsize = FALSE)

cells <- unique(noise_contribution_data_filtered %>% 
  filter(dataset=="glycerol040.rplN.2.20230413") %>% 
  filter(datatype=="nonoise") %>% 
  filter(concentration>45000) %>% 
  distinct(cell) %>% 
  .$cell)

old_cells_filtered %>%
  filter(cell %in% cells)

old_cells_filtered_filtered <- old_cells_filtered %>% 
  filter(cell %in% cells) %>% 
  mutate(lineage=row_number()) %>% 
  pivot_longer(cols=c("cell","mother","g_mother","gg_mother","ggg_mother"),values_to = c("cells")) %>% 
  rename(cell=cells,type=name)
  
noise_contribution_data %>% 
  mutate(dataset=paste(condition,promoter,replica,date,sep=".")) %>% 
  filter(dataset=="glycerol040.rplN.2.20230413") %>% 
  filter(datatype=="nonoise") %>% 
  semi_join(old_cells_filtered_filtered,by=c("cell")) %>% 
  left_join(old_cells_filtered_filtered,by=c("cell")) %>% 
  ggplot()+
  geom_line(aes(time_sec,concentration,col=cell),alpha=1)
  
