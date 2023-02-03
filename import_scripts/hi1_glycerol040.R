## ----message=FALSE,warning=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source("../mother_machine_experiments_toolbox/load_functions_and_packages.R")
path_to_MM_data_summary <- "./datalist_hi1_glycerol040.csv"
source("../mother_machine_experiments_toolbox/read_MM_data.R")
source("../mother_machine_experiments_toolbox/transformMMData.R") 


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
.promoter <- "hi1"
.condition <- "glycerol040"
.date <- 20190515


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
myframes %>% 
  distinct(promoter,condition,date,pos,gl_number) %>% 
  group_by(condition,promoter,date,pos) %>% 
  arrange(pos) %>% 
  summarise(nb_growth_lanes=n())


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_dataset_statistics(.condition,.promoter,FALSE)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_dataset_statistics(.condition,.promoter,TRUE)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
wrong_cells <- myframes_complete %>% 
  ungroup() %>% 
  mutate(wrong=FALSE) %>% 
  group_by(cell) %>% 
  arrange(frame) %>% 
  #mutate(dlength=length_um_oriented_bbox-lag(length_um_oriented_bbox)) %>% 
  #mutate(dlength=ifelse(is.na(dlength),0,dlength)) %>% 
  #mutate(wrong=ifelse(dplyr::last(dlength)<0,TRUE,wrong)) %>% 
  ungroup() %>% 
  mutate(wrong=ifelse(if_any(.cols=everything(),  ~is.na(.)),TRUE,wrong)) %>% 
  filter(wrong==TRUE) %>% 
  distinct(cell)

mycells_complete <- myframes_complete %>% 
  anti_join(wrong_cells,by=c("cell")) %>% 
  #filter(gl_number<=30) %>% 
  group_by(cell) %>% 
  arrange(frame) %>% 
  filter(n()>8) %>% 
  #filter(logl_time_r2>0.9) %>% 
  #filter(alpha_oriented_bbox*3600/log(2)<0.5) %>% 
  mutate(added_length=last(length_um_oriented_bbox)-first(length_um_oriented_bbox)) %>% 
  mutate(l_birth=first(length_um_oriented_bbox)) %>% 
  mutate(l_div=first(length_um_oriented_bbox)) %>% 
  mutate(mean_c=mean(gfp_nb/length_um_oriented_bbox,na.rm=TRUE)) %>% 
  ungroup() %>% 
  distinct(cell,.keep_all = TRUE)
mycells_complete$mean_cell_rank <- factor(mycells_complete$mean_cell_rank,levels=as.character(sort(unique(as.integer(mycells_complete$mean_cell_rank)))))

#Filtering out cells whose size is decreasing at the end of the cell cycle.
cells_discarded_manually <- c("1.4.223","4.263","4.6.275","0.9.260","4.7.178","5.2.291","0.12.265")
cells_discarded_manually <- lapply(cells_discarded_manually,function(.l) paste(as.character(.date),as.character(.l),sep=".")) %>% unlist()
cells_discarded_manually <- tibble(cell=cells_discarded_manually)
mycells_complete <- anti_join(mycells_complete,cells_discarded_manually,by="cell",copy=TRUE)

good_cells <- mycells_complete %>% 
  select(cell)

# Keep only the good cells
myframes_complete_goodcells <- myframes_complete %>% 
  semi_join(good_cells,by=c("cell")) %>% 
  arrange(cell,time_sec)

# Full cell cycle observables
mycells_complete_goodcells <- mycells_complete %>% 
  semi_join(good_cells,by=c("cell"))

readr::write_csv(myframes_complete_goodcells,sprintf("../raw_curated_data_complete_cycles/%s_%s_rawdata.csv",.promoter,.condition))



## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cells_discarded_manually <- c("1.4.223","4.263","4.6.275","0.9.260","4.7.178","5.2.291","0.12.265")

cells_discarded_manually <- lapply(cells_discarded_manually,function(.l) paste(as.character(.date),as.character(.l),sep=".")) %>% unlist()
cells_discarded_manually <- tibble(cell=cells_discarded_manually)

# Discarding cells with too few observations or very bad r2
myframes_all <- myframes %>% 
  anti_join(cells_discarded_manually,by="cell",copy=TRUE) %>% 
  group_by(cell) %>% 
  arrange(frame) %>% 
  filter(n()>8) %>% 
  na.omit() %>% 
  ungroup() %>% 
  arrange(cell,frame)

readr::write_csv(myframes_all,sprintf("../raw_curated_data/%s_%s_rawdata.csv",.promoter,.condition))


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mean_exponential_gr_min <- unique(myframes_complete %>% 
  semi_join(good_cells,by=c("cell")) %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  ungroup() %>% 
  mutate(mean_exp_gr=mean(alpha_oriented_bbox*60)) %>% 
  .$mean_exp_gr)

mean_concentration <- unique(myframes_complete %>% 
  semi_join(good_cells,by=c("cell")) %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  ungroup() %>% 
  group_by(cell) %>% 
  mutate(mean_concentration=mean(gfp_nb/length_um_oriented_bbox)) %>% 
  ungroup() %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  mutate(mean_mean_concentration=mean(mean_concentration)) %>% 
  .$mean_mean_concentration)

mean_production <- mean_concentration*mean_exponential_gr_min

print(sprintf("Mean exponential growth-rate: %s min-1",signif(mean_exponential_gr_min,3)))
print(sprintf("Mean production: %s GFP/min/um",signif(mean_production,3)))
print(sprintf("Mean concentration: %s GFP/um",signif(mean_concentration,3)))


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(plot==TRUE){
for (i in c(0:4)){
  good_cells_fraction <- good_cells %>% 
    filter(dplyr::between(row_number(),nrow(good_cells)/5*i,nrow(good_cells)/5*(i+1))==TRUE)
           
  good_data <- myframes_complete %>% 
    semi_join(good_cells_fraction,by=c("cell"))
  
  print(
  good_data %>%
    ggplot()+
    geom_point(aes(frame,log(length_um_oriented_bbox))) +
    facet_wrap(~cell,scale="free"))
  
ggsave(sprintf("../raw_curated_data_complete_cycles/%s_%s_%i.pdf",.promoter,.condition,i),width=50,height=50,limitsize=FALSE)}}


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
plot_dataset_statistics_curated(.condition,.promoter,TRUE)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mycells_complete_goodcells %>% 
  group_by(mean_cell_rank) %>% 
  count()


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mycells_complete_goodcells %>% 
  ggplot()+
  stat_ecdf(aes(alpha_oriented_bbox*3600/log(2),col=pos))+
  facet_wrap(~paste("Cell rank",mean_cell_rank,sep =" "),ncol=1)+
  guides(col = guide_legend(title = "Position"))+
  xlim(0,1)+
  xlab("Growth rate (h-1)")+
  ylab("CDF")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mycells_complete_goodcells %>% 
  ggplot()+
  stat_ecdf(aes(alpha_oriented_bbox*3600/log(2),col=mean_cell_rank))+
  facet_wrap(~paste("Pos",pos,sep =" "),ncol=1)+
  guides(col = guide_legend(title = "Cell rank"))+
  xlim(0,1)+
  xlab("Growth rate (h-1)")+
  ylab("CDF")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mycells_complete_goodcells %>% 
  ggplot()+
  stat_ecdf(aes(alpha_oriented_bbox*3600/log(2),col=paste(mean_cell_rank,"x",pos,sep=" ")))+
  #facet_wrap(~paste("Pos",pos,sep =" "),ncol=1)+
  guides(col = guide_legend(title = "Cell rank x Position"))+
  scale_color_manual(values=c25)+
  xlim(0,0.7)+
  xlab("Growth rate (h-1)")+
  ylab("CDF")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mycells_complete_goodcells %>% 
  group_by(mean_cell_rank) %>% 
  mutate(mean_gr=mean(alpha_oriented_bbox*3600/log(2))) %>% 
  mutate(sd_gr=sd(alpha_oriented_bbox*3600/log(2))) %>% 
  mutate(cv_gr=sd_gr/mean_gr) %>% 
  ungroup() %>% 
  distinct(mean_cell_rank,mean_gr,sd_gr,cv_gr) %>% 
  mutate(rank_0_gr=dplyr::first(mean_gr),
         gr_ratio=mean_gr/rank_0_gr) %>% 
  select(-rank_0_gr)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mycells_complete_goodcells %>% 
  ggplot()+
  stat_ecdf(aes(logl_time_r2,col=paste(mean_cell_rank,"x",pos,sep=" ")))+
  #facet_wrap(~paste("Pos",pos,sep =" "),ncol=1)+
  guides(col = guide_legend(title = "Cell rank x Position"))+
  scale_color_manual(values=c25)+
  xlim(0.9,1)+
  xlab("log(length) versus time correlation (R2)")+
  ylab("CDF")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mycells_complete_goodcells %>% 
  ggplot()+
  geom_point(aes(alpha_oriented_bbox*3600/log(2),logl_time_r2,col=paste(mean_cell_rank,"x",pos,sep=" ")))+
  #facet_wrap(~paste("Pos",pos,sep =" "),ncol=1)+
  guides(col = guide_legend(title = "Cell rank x Position"))+
  scale_color_manual(values=c25)+
  #xlim(0,1)+
  xlab("Growth rate (h-1)")+
  ylab("Correlation")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mycells_complete_goodcells %>% 
  ggplot()+
  geom_point(aes(l_birth,added_length,col=paste(mean_cell_rank,"x",pos,sep=" ")))+
  #facet_wrap(~paste("Pos",pos,sep =" "),ncol=1)+
  guides(col = guide_legend(title = "Cell rank x Position"))+
  scale_color_manual(values=c25)+
  xlim(0,4)+
  ylim(0,4)+
  xlab("Length at birth")+
  ylab("Added length")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mycells_complete_goodcells %>% 
  ggplot()+
  geom_point(aes(l_birth,added_length,col=paste(mean_cell_rank,"x",pos,sep=" ")))+
  #facet_wrap(~paste("Pos",pos,sep =" "),ncol=1)+
  guides(col = guide_legend(title = "Cell rank x Position"))+
  scale_color_manual(values=c25)+
  facet_wrap(~paste(mean_cell_rank,"x",pos,sep=" "))+
  xlim(0,4)+
  ylim(0,4)+
  xlab("Length at birth")+
  ylab("Added length")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mycells_complete_goodcells %>% 
  ggplot()+
  geom_point(aes(alpha_oriented_bbox*3600/log(2),added_length,col=paste(mean_cell_rank,"x",pos,sep=" ")))+
  #facet_wrap(~paste("Pos",pos,sep =" "),ncol=1)+
  guides(col = guide_legend(title = "Cell rank x Position"))+
  scale_color_manual(values=c25)+
  xlim(0,1)+
  ylim(0,4)+
  xlab("Growth-rate (h-1)")+
  ylab("Added length")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
unique(mycells_complete_goodcells %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  ungroup() %>% 
  mutate(mean_doubling_time=mean(1/(alpha_oriented_bbox*3600/log(2)))) %>% 
  .$mean_doubling_time)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mycells_complete_goodcells %>% 
  group_by(cell) %>% 
  arrange(frame) %>% 
  mutate(median_time=first(time_sec)+(last(time_sec)-first(time_sec))/2) %>% 
  ungroup() %>% 
  select(median_time,alpha_oriented_bbox,pos,mean_cell_rank,gl_number) %>% 
  ggplot()+
  geom_point(aes(median_time/3600,alpha_oriented_bbox*3600/log(2),col=gl_number)) +
  facet_wrap(~pos,ncol=1)+
  ylim(0,1)+
  xlab("Time (h)")+
  ylab("Growth rate (h-1)")


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
myframes_complete_goodcells %>% 
  group_by(cell) %>% 
  arrange(frame) %>% 
  mutate(median_time=first(time_sec)+(last(time_sec)-first(time_sec))/2) %>% 
  ungroup() %>% 
  select(cell,median_time,alpha_oriented_bbox,pos,mean_cell_rank,gl_number) %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  ggplot()+
  geom_point(aes(median_time/3600,alpha_oriented_bbox*3600/log(2),col=gl_number)) +
  facet_wrap(~paste(pos,gl_number),ncol=1)+
  ylim(0,1)+
  xlab("Time (h)")+
  ylab("Growth rate (h-1)")

#ggsave("./hi1_acetate2_growth_rate_vs_time_all_gl.pdf",width=20,height=50,limitsize=FALSE)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
.gl_min <- 0 
.gl_max <- 30
plot_concentration_distribution(.date,.gl_min,.gl_max)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
selected_cells <- myframes_all %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  group_by(condition,promoter) %>% 
  sample_n(50) %>% 
  ungroup() %>% 
  distinct(cell)

myframes_all %>% 
  semi_join(selected_cells,by=c("cell")) %>% 
  ggplot()+
  geom_point(aes(-vertical_top,length_um))+
  #geom_line(aes(time_sec,length_um))+
  facet_wrap(~paste(cell,condition),scale="free")

ggsave("./cell_detection_offset_glycerol040.pdf",height=20,width=20)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
myframes_all %>% 
  mutate(wrong=FALSE) %>%
  group_by(cell) %>% 
  arrange(frame) %>% 
  mutate(dlength=length_um_oriented_bbox-lag(length_um_oriented_bbox)) %>% 
  mutate(dlength=ifelse(is.na(dlength),0,dlength)) %>% 
  mutate(wrong=ifelse(dplyr::last(dlength)<0,TRUE,wrong)) %>% 
  ungroup() %>% 
  mutate(wrong=ifelse(if_any(.cols=everything(),  ~is.na(.)),TRUE,wrong)) %>% 
  filter(wrong==TRUE) %>% 
  #distinct(cell)
  ggplot()+
  geom_point(aes(-vertical_top,length_um))+
  facet_wrap(~cell)
#facet_wrap(~paste(cell,condition),scale="free")

ggsave("./length_time_glycerol040.pdf",height=20,width=20)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
myframes_all %>% 
  filter(cell=="20190515.0.9.260") %>% 
  mutate(wrong=FALSE) %>%
  group_by(cell) %>% 
  arrange(frame) %>% 
  mutate(dlength=length_um_oriented_bbox-lag(length_um_oriented_bbox)) %>% 
  mutate(dlength=ifelse(is.na(dlength),0,dlength)) %>% 
  mutate(wrong=ifelse(dplyr::last(dlength)<0,TRUE,wrong)) %>% 
  ungroup() %>% 
  mutate(wrong=ifelse(if_any(.cols=everything(),  ~is.na(.)),TRUE,wrong)) %>% 
  filter(wrong==TRUE) %>% 
  #distinct(cell)
  ggplot()+
  geom_point(aes(-vertical_top,length_um))+
  facet_wrap(~cell)
#facet_wrap(~paste(cell,condition),scale="free")

