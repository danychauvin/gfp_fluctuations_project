# IMPORT DEEP MOMA DATA INTO R DATAFRAME
# Author: Dany Chauvin

# Set defaults variables
dl <- 0.065 # µm/pixel

# Set up a parallel environments for multidyplr

if(!exists("mycluster")){
mycluster <- min(30, parallel::detectCores()-1) %>%  # do not use more than 30 cores
  default_cluster() %>% cluster_library( # load currently loaded packages on each core
    names(sessionInfo()$otherPkgs))
  
cluster_copy(mycluster,c("fit_exp_elongation_predict","fit_exp_elongation_slope","fit_exp_elongation_sd_slope","fit_exp_elongation_intercept"))
}

# Importing conditions from dataList

#path_to_MM_data_summary <- '/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/mother_machine_experiments_toolbox/all_mmexperiments_datalist.csv'

myconditions <- readr::read_csv(path_to_MM_data_summary,
                               col_types = cols(
                                 date=col_character(),
                                 description=col_character(),
                                 f_start=col_character(),
                                 f_end=col_character(),
                                 condition=col_character(),
                                 t_interval=col_double(),
                                 data_path=col_character(),
                                 #folder=col_character(),
                                 promoter=col_character(),
                                 vector=col_character(),
                                 cell_detection_offset=col_double(),
                                 experimental_setup=col_character()))

# Generate the list of curated data files to fetch
# Computing t_start and t_end in seconds
myconditions <- myconditions %>% 
  mutate(f_start=as.double(f_start),
         f_end=as.double(f_end)) %>% 
  group_by(date,description,data_path) %>% 
  arrange(f_start) %>%
  mutate(t_end=ifelse(row_number()==1,(f_end)*t_interval*60,NA)) %>% 
  mutate(t_start=ifelse(row_number()==1,(f_start)*t_interval*60,NA)) %>% 
  do((function(.df){
    new_df <- .df
    if(nrow(new_df)>=2){
      for(i in c(2:nrow(new_df))) 
      {
        new_df$t_start[i]=new_df$t_end[i-1]
        new_df$t_end[i]=new_df$t_start[i]+(new_df$f_end[i]-new_df$f_start[i])*new_df$t_interval[i]*60}
    }
    return(new_df)
  })(.)) %>%
  ungroup()

#Find the files
myconditions <- myconditions %>% 
  group_by(date,description,data_path) %>% 
  do((function(.df){
    data_folder <- unique(.df$data_path)
    new_df <- find.files(data_folder, "*CellStats_*.csv",.mindepth=1,.maxdepth = 2) %>% 
      data.frame(file_path=., stringsAsFactors=FALSE)
    new_df <- crossing(.df,new_df)
    return(new_df)
  } )(.)) %>% 
  ungroup()

# Fetching information about positions to propagate them in the final data table
#mypositions <- myconditions %>% 
#  ungroup() %>% 
#  distinct(data_path) %>% 
#  .$data_path %>% 
#  lapply(function(.l) as_tibble(data.table::fread(paste(.l,"/positions.csv",sep=""),skip=0,sep=","))) %>% 
#  do.call(rbind, .) %>% 
#  mutate_at(vars(date), factor) %>% 
#  mutate_at(vars(pos), factor)
  

#Define proper column names
nFiles <- myconditions %>% distinct(file_path) %>% nrow()
properColNames <- c("lane_ID","cell_ID","frame","cell_rank","genealogy","type_of_end","parent_ID","cells_in_lane","bbox_top_px","bbox_bottom_px","touches_detection_roi_top","center_x_px","center_y_px","width_px","length_px","tilt_rad","area_px","bgmask_area_px","phc_total_intensity_au","phc_intensity_coefficient_of_variation","label:dead","label:dying","label:fading","label:shrinking","label:dividing","fluo_cellmask_ch_1","fluo_bgmask_ch_1","fluo_ampl_ch_1","fluo_bg_ch_1","fluo_cellmask_ch_2","fluo_bgmask_ch_2","fluo_ampl_ch_2","fluo_bg_ch_2","oriented_bbox_center_x_px","oriented_bbox_center_y_px","oriented_bbox_width_px","oriented_bbox_length_px","oriented_bbox_orientation_angle_rad","contour_area__px","contour_centroid_x__px","contour_centroid_y__px","contour_variance_x__px2","contour_variance_y__px2","contour_covariance__px2")
properColNames <- c(properColNames,"file_path")

nr <- function(.l){
  .df <- as_tibble(data.table::fread(.l,skip=4,sep=","))
  return(c(.l,nrow(.df)))}

#Reading curated data table content
myframes <- myconditions %>% 
  ungroup() %>% 
  distinct(file_path) %>% 
  .$file_path %>% 
  lapply(function(.l) .df <- as_tibble(data.table::fread(.l,skip=4,sep=",")) %>%
  mutate(file_path=.l) %>% 
  mutate(type_of_end=ifelse(type_of_end=="",NA,type_of_end))) %>% 
  lapply(function(.l) .l %>% select(all_of(properColNames))) %>% 
           do.call(rbind, .)

#Parsing positions and growth lanes number
myframes <- myframes %>%
  left_join(myconditions,by=c('file_path')) %>% 
  group_by(file_path,condition) %>% 
  #filter(data.table::between(frame,f_start,f_end-1)) %>% 
  filter(frame>=f_start) %>% 
  tidyr::extract(file_path, c("pos", "gl"), "__\\d{8}[glugly]*_.*[Pp]os(\\d+)_GL(\\d+).*.csv$", remove=FALSE, convert=TRUE) %>%
  mutate(gl=as.integer(gl)) %>% 
  mutate(pos=as.integer(pos)) %>% 
  mutate(date=as.integer(date)) %>% 
  ungroup()

#Renaming columns and computing useful variables
myframes <- myframes %>% 
  rename(id=cell_ID,
         end_type=type_of_end,
         parent_id=parent_ID,
         vertical_bottom='bbox_bottom_px',
         vertical_top='bbox_top_px',
         center_x_pixel='center_x_px',
         center_y_pixel='center_y_px',
         width_pixel='width_px',
         length_pixel='length_px',
         tilt_radian='tilt_rad',
         area_pixel2='area_px',
         background_mask_aread_pixel2='bgmask_area_px',
         fluo_cell_mask_ch1='fluo_cellmask_ch_1',
         fluo_background_mask_ch1='fluo_bgmask_ch_1',
         fluo_amplitude='fluo_ampl_ch_1',
         fluo_background_ch1='fluo_bg_ch_1') %>% 
  #Convert date to factors
  mutate_at(vars(date), factor) %>% 
  mutate_at(vars(pos), factor) %>% 
  mutate_at(vars(gl), factor) %>% 
  #Set cell ref
  #mutate(cell=paste(date,pos,gl,id,sep='.'),
  mutate(cell=paste(date,pos,gl,id,sep='.'),
         parent=paste(date,pos,gl,parent_id,sep='.'),
         lane_ID=paste(date,pos,gl,sep='.')) %>% 
  #Compute biomass estimations
  mutate(length_um_oriented_bbox=oriented_bbox_length_px*dl,
         width_um_oriented_bbox=oriented_bbox_width_px*dl,
         length_um=length_um_oriented_bbox) %>%
  #Compute vertical center
  mutate(vertical_center=vertical_top+(vertical_bottom - vertical_top)/2) %>%
  # Propagating end_type information
  mutate(gl_id=paste(date, pos,gl, sep='.')) %>% 
  tidyr::extract(gl_id,"gl_number","^20[0-9]{6}.[0-9]{1,}.([0-9]{1,})$",remove="FALSE",convert="FALSE") %>% 
  mutate(gl_number=as.double(gl_number)) %>% 
  group_by(cell,condition) %>% 
  mutate(pruned=ifelse(end_type=="pruned",TRUE,FALSE)) %>% 
  mutate(pruned=ifelse(is.na(end_type)==TRUE,FALSE,pruned)) %>%
  fill(end_type,.direction="up") %>%
  mutate(discard_top=(vertical_top <= cell_detection_offset)) %>%
  #mutate(discard_top=ifelse(vertical_center < 150,TRUE,discard_top)) %>%
  mutate(end_type=ifelse(any(discard_top), 'touchtop', end_type)) %>%
  ungroup()

#Computing single-cell growth-rate based on oriented bounding box length, also define cell_rank_beg, the rank of the cell at birth
myframes <- myframes %>% 
  #filter(!(cell %in% wrong_cells)) %>%
  mutate(time_sec=t_start+(frame-f_start)*t_interval*60) %>% 
  group_by(cell,condition) %>% 
  partition(cluster=mycluster) %>%
  mutate(cell_rank_beg=first(cell_rank),
         mean_cell_rank=round(mean(cell_rank),0)) %>% 
  mutate(start_time=first(time_sec), end_time=last(time_sec),
         length_predict_oriented_bbox=fit_exp_elongation_predict(time_sec, length_um_oriented_bbox),
         alpha_oriented_bbox=fit_exp_elongation_slope(time_sec,length_um_oriented_bbox),
         logl_time_r2=cor(time_sec, log(length_um_oriented_bbox))^2,
         sd_alpha_oriented_bbox=fit_exp_elongation_sd_slope(time_sec,length_um_oriented_bbox)
         ) %>%
  collect() %>% 
  arrange(cell, frame) %>% # sort data after `partition()`
  ungroup()

#Add the mean_position of the cell between two divisions in the growth lane.
myframes <- myframes %>% 
  group_by(cell,condition) %>%
  mutate(mean_position_um=mean(vertical_center)*dl) %>% 
  ungroup()

# Add column switch to cells that are experiencing a switch
## LARGE NUMBER OF CELLS LOST HERE
myframes <- myframes %>%
  group_by(cell,condition) %>% 
  arrange(frame) %>% 
  mutate(switch_frame=ifelse((time_sec-t_start)<time_threshold_h*3600,TRUE,FALSE)) %>% 
  mutate(switch=ifelse(any((time_sec-t_start)<time_threshold_h*3600),TRUE,FALSE)) %>% 
  mutate(switch_frame=ifelse(time_sec>=t_end,TRUE,switch_frame)) %>% 
  mutate(switch=ifelse(any(time_sec>=t_end),TRUE,switch)) %>% 
  mutate(mean_vertical_center=mean(vertical_center,na.rm=TRUE)) %>% 
  ungroup()

# An observation "i" (i.e a row), belonging to a cell cycle "c" is in myframes_complete if:
# - myframes_complete contain all other observations about "c", from birth to division.
# - "c" 

#Searching for cells that show NAs
cells_with_nas <- myframes[!complete.cases(myframes),] %>% 
  distinct(cell)

#used to compute myframes all: frames that should be discarded: touching the top, pruned, or within a switch.
print(sprintf("%i cells initially.",nrow(myframes %>% distinct(cell))))
#print(sprintf("%i cells with NA's.",nrow(cells_with_nas)))

myframes <- myframes %>% 
  filter(switch_frame==FALSE) %>% 
  filter(discard_top==FALSE) %>% 
  filter(pruned==FALSE) %>% 
  #also discarding cell cycles with NA's in
  anti_join(cells_with_nas,by=c("cell"))

print(sprintf("%i cells after filtering for switch, touchtop, pruned and cells with NAs.",nrow(myframes %>% distinct(cell))))

#discard segmentation issue that were not seen during curation
bad_cells <- myframes %>% 
  arrange(date,condition,promoter,cell,frame) %>% 
  group_by(date,condition,promoter,cell) %>% 
  mutate(delta_size=abs(log(length_um_oriented_bbox)-log(lag(length_um_oriented_bbox)))) %>% 
  ungroup() %>% 
  drop_na() %>% 
  filter(delta_size>0.2) %>% 
  distinct(cell,condition)

myframes <- myframes %>% 
  anti_join(bad_cells,by=c("cell","condition"))

print(sprintf("%i cells after filtering for probable segmentation issues.",nrow(myframes %>% distinct(cell))))

#cells that should be kept: not belonging to a switch, not touching the top, not exiting, not eod, not pruned, only dividing, we should see its birth.
#with one frame belonging to a switch (switch==TRUE).
myframes_complete <- myframes %>% 
  filter(switch==FALSE) %>% 
  group_by(cell,condition) %>% 
  filter(parent_id!=-1) %>% #we should see the birth of the cell
  filter(end_type=='div') %>%  #%>% #we should see cell division
  #filter(n()>10) %>% 
  #filter(logl_time_r2>0.9) %>% 
  #filter(switch==FALSE) %>% #if an observation occurs during a switch (2h< after the beginning of the experiment, or after a switch), the full cell cycle is discarded.
  ungroup()

print(sprintf("%i cells after filtering for complete cell cycles.",nrow(myframes_complete %>% distinct(cell))))

dividing_times <- myframes_complete %>%
  group_by(cell,condition) %>% 
  arrange(time_sec) %>% 
  mutate(div_time=(last(time_sec)-first(time_sec))/3600) %>% 
  ungroup() %>% 
  distinct(cell,.keep_all = TRUE) %>% 
  group_by(date,promoter,condition) %>% 
  mutate(mean_division_time=mean(div_time,na.rm=TRUE)) %>% 
  ungroup() %>% 
  distinct(date,condition,promoter,mean_division_time)

# Filtering out filamenting and slow growing cells
# myframes_complete <- myframes_complete %>% 
#   left_join(dividing_times,by=c("date","condition","promoter")) %>% 
#   group_by(cell) %>% 
#   arrange(time_sec) %>% 
#   mutate(div_time=(last(time_sec)-first(time_sec))/3600) %>% 
#   ungroup()
#   #filter(div_time<4*mean_division_time)

#print(sprintf("%i cells after filtering for growing cells.",nrow(myframes_complete %>% distinct(cell))))

# myframes_complete <- myframes_complete %>% 
#   group_by(date,promoter,condition) %>% 
#   mutate(mean_cell_length=mean(length_um,na.rm=TRUE),
#          sd_cell_length=sd(length_um,na.rm=TRUE),
#          min_threshold=mean_cell_length-4*sd_cell_length,
#          max_threshold=mean_cell_length+4*sd_cell_length) %>% 
#   ungroup() %>% 
#   mutate(filamenting=ifelse(length_um>max_threshold,TRUE,FALSE)) %>% 
#   group_by(cell) %>% 
#   filter(!any(filamenting)) %>% 
#   ungroup()

#print(sprintf("%i cells after filtering out filamenting cells.",nrow(myframes_complete %>% distinct(cell))))

# Filter cell rank characterized by a too small number of cells
#myframes_complete <- myframes_complete %>% 
#  mutate(gl_type=ifelse(((condition=="acetate") & (gl_number>30)),"long","short"))

# Unused filter
#cell_ranks_to_keep <- myframes_complete %>% 
#  distinct(date,promoter,condition,date,pos,gl_number,cell,mean_cell_rank) %>% 
#  group_by(date,condition,promoter,mean_cell_rank) %>% 
#  arrange(pos) %>% 
#  count() %>% 
#  filter(n>50) %>% 
#  select(-n) %>% 
#  distinct(date,condition,promoter,mean_cell_rank)

#myframes_complete <- myframes_complete %>% 
#  semi_join(cell_ranks_to_keep,by=c("date","promoter","condition","mean_cell_rank"))

mycells <- myframes_complete %>% 
  group_by(date,pos,gl,cell,description,condition,vector) %>% 
  partition(cluster=mycluster) %>%
  do((function(.df) {
    .mod_ll_t <- fastLmPure( cbind(1, .df$time_sec), log(.df$length_um_oriented_bbox) )
    .mod_l_t <- fastLmPure( cbind(1, .df$time_sec), .df$length_um_oriented_bbox)
    .time_birth <- first(.df$time_sec)
    .time_div <- last(.df$time_sec)
    .cell_rank_beg <- first(.df$cell_rank)
    .alpha_oriented_bbox=first(.df$alpha_oriented_bbox)
    .sd_alpha_oriented_bbox=first(.df$sd_alpha_oriented_bbox)
    data.frame(npoints=.mod_ll_t$df.residual+1,
               time_birth=.time_birth, time_div=.time_div, div_time=.time_div-.time_birth,
               l_birth=first(.df$length_predict_oriented_bbox), l_div=last(.df$length_predict_oriented_bbox),
               w_birth=first(.df$width_um_oriented_bbox),w_div=last(.df$width_um_oriented_bbox),
               logl_time_slope=.mod_ll_t$coefficients[2], logl_time_slopesd=.mod_ll_t$stderr[2],
               logl_time_r2=cor(.df$time_sec, log(.df$length_um_oriented_bbox))^2,
               l_time_slope=.mod_l_t$coefficients[2], l_time_slopesd=.mod_l_t$stderr[2], 
               l_time_r2=cor(.df$time_sec, .df$length_um_oriented_bbox)^2,
               cell_rank_beg=.cell_rank_beg,
               alpha_oriented_bbox=.alpha_oriented_bbox,sd_alpha_oriented_bbox=.sd_alpha_oriented_bbox)
  })(.) ) %>% 
  collect() %>% 
  #arrange(condition, date, pos, gl, id) %>% 
  ungroup()

mycells <- mycells %>% 
  # filter by r2 of exponential fit
  #filter(logl_time_r2>0.95) %>% 
  # create new variables of interest
  mutate(dl=l_div - l_birth,
         alpha_len=log(l_div/l_birth) / div_time)
# ok since l_xxx are fitted
#c_birth=g_birth/l_birth,
#c_div=g_div/l_div,
#dg=g_div - g_birth,
#dcdt=(g_div/l_div - g_birth/l_birth) / div_time,
#g=log(g_div/g_birth) / div_time, # would be cleaner from a fit but difficult to compare
#gamma=dg/dl,
#q=dg/dl*alpha)

myframes %>% distinct(cell) %>% nrow()
myframes_complete %>% distinct(cell) %>% nrow()
mycells %>%  distinct(cell) %>% nrow()
myframes_complete <- myframes_complete %>% 
  semi_join(mycells,by=c("cell"))


#myframes <- myframes %>%
#  left_join(mypositions %>% 
#              mutate(pos=as.integer(pos)) %>% 
#              select(date,ch,pos,width) %>%
#              mutate_at(vars(date), factor) %>% 
#              mutate_at(vars(pos), factor) %>% 
#              select(date,ch,pos,width),by=c("date","pos"))
  
#myframes_complete <- myframes_complete %>% 
#  left_join(mypositions %>% 
#              mutate(pos=as.integer(pos)) %>% 
#              select(date,ch,pos,width) %>%
#              mutate_at(vars(date), factor) %>% 
#              mutate_at(vars(pos), factor),
#            by=c("date","pos"))

conditions=c("acetate005","glycerol040","glucose020","glucoseaa020")
myframes$condition <- factor(myframes$condition,levels=conditions)
myframes_complete$condition <- factor(myframes_complete$condition,levels=conditions)

#Propagate fluorescence data
# Setting strains
myframes_complete <- myframes_complete %>% 
  mutate(strain='MG1655')
# Converting GFP units
myframes_complete <- myframes_complete %>%
  # append relevant conversion parameters
  mutate(autofluo_predict = 133.6,
         fp_per_dn = 0.198) %>%  #zaslaver library, mg1655
  # convert to gfp units (after subtracting autofluorescence)
  #mutate(fluogfp_amplitude = fluo_amplitude - autofluo_predict * length_um_oriented_bbox,
  #       gfp_nb = fluogfp_amplitude * fp_per_dn ) %>% 
  mutate(gfp_nb = fluo_amplitude * fp_per_dn) %>% 
  group_by(cell,condition) %>% 
  arrange(time_sec) %>% 
  mutate(g_birth=first(gfp_nb)) %>% 
  mutate(g_div=last(gfp_nb)) %>% 
  ungroup()

myframes <- myframes %>% 
  mutate(strain='MG1655')
# Converting GFP units
myframes_complete <- myframes_complete %>%
  # append relevant conversion parameters
  mutate(autofluo_predict = 133.6,
         fp_per_dn = 0.198) %>%  #zaslaver library, mg1655
  # convert to gfp units (after subtracting autofluorescence)
  #mutate(fluogfp_amplitude = fluo_amplitude - autofluo_predict * length_um_oriented_bbox,
  #       gfp_nb = fluogfp_amplitude * fp_per_dn ) %>% 
  mutate(gfp_nb = fluo_amplitude * fp_per_dn)
  
myframes_complete <- myframes_complete %>% 
  ungroup() %>% 
  mutate(dataset=paste(condition,promoter,date,sep="."))

print(sprintf("%i cells finally.",nrow(myframes_complete %>% distinct(cell))))

