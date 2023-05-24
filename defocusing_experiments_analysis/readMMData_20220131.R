# IMPORT DEEP MOMA DATA INTO R DATAFRAME
# Author: Dany Chauvin
# Date: 20220130
# Version which enables to imporing of the spines.

# Input
# Path to the dataframe in the .csv file which recapitulates experimental details.
# Ouput
# mycells and myframes_to_mycells dataframes

# SET DEFAULT VARIABLES
dl <- 0.065 # µm/pixel
vertical_cutoff <- 150 #Necessary to get rid of cells that are touching the top of the growth lane where segmentation is not that good anymore.

# Parallel environment for multidplyr
if(!exists("mycluster")){
mycluster <- min(30, parallel::detectCores()-1) %>%  # do not use more than 30 cores
  default_cluster() %>% cluster_library( # load currently loaded packages on each core
    names(sessionInfo()$otherPkgs))
  
cluster_copy(mycluster,c("fit_exp_elongation_predict","fit_exp_elongation_slope","fit_exp_elongation_sd_slope","fit_exp_elongation_intercept"))
}

myconditions <- readr::read_csv(path_to_MM_data_summary,
                               col_types = cols(
                                 date=col_character(),
                                 description=col_character(),
                                 f_start=col_character(),
                                 f_end=col_character(),
                                 condition=col_character(),
                                 t_interval=col_double(),
                                 orientation=col_character(),
                                 data_path=col_character(),
                                 promoter=col_character(),
                                 vector=col_character()))

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

myconditions <- myconditions %>% 
  group_by(date,description,data_path) %>% 
  do((function(.df){
    new_df <- find.files(unique(.df$data_path), "ExtendedExportedCellStats_*.csv") %>% 
      data.frame(file_path=., stringsAsFactors=FALSE)
    new_df <- crossing(.df,new_df)
    return(new_df)
  } )(.)) %>% 
  ungroup()

nFiles <- myconditions %>% distinct(file_path) %>% nrow()

#nc <- min(nFiles, length(mycluster)) # this is a dirty hack because multidplyr crashes with less shards than cores
#Proper column names

properColNames <- c("lane_ID","cell_ID","frame","cell_rank","genealogy","type_of_end","parent_ID","cells_in_lane","bbox_top_px","bbox_bottom_px","touches_detection_roi_top","center_x_px","center_y_px","width_px","length_px","tilt_rad","area_px","bgmask_area_px","phc_total_intensity_au","phc_intensity_coefficient_of_variation","label:dead","label:dying","label:fading","label:shrinking","label:dividing","fluo_cellmask_ch_1","fluo_bgmask_ch_1","fluo_ampl_ch_1","fluo_bg_ch_1","fluo_cellmask_ch_2","fluo_bgmask_ch_2","fluo_ampl_ch_2","fluo_bg_ch_2","oriented_bbox_center_x_px","oriented_bbox_center_y_px","oriented_bbox_width_px","oriented_bbox_length_px","oriented_bbox_orientation_angle_rad","contour_area__px","contour_centroid_x__px","contour_centroid_y__px","contour_variance_x__px2","contour_variance_y__px2","contour_covariance__px2","fl_intensity_profile")
properColNames <- c(properColNames,"file_path")


myframes <- myconditions %>% 
  ungroup() %>% 
  distinct(file_path) %>% 
  .$file_path %>% 
  #lapply(function(.l) readr::read_delim(.l,skip=2,delim=";",show_col_types =FALSE) %>% 
           #mutate(file_path=.l)) %>% 
  lapply(function(.l) as_tibble(data.table::fread(.l,skip=0,sep=",")) %>% 
  mutate(file_path=.l) %>% 
  mutate(type_of_end=ifelse(type_of_end=="",NA,type_of_end))) %>% 
  lapply(function(.l) .l %>% select(all_of(properColNames))) %>% 
           do.call(rbind, .)
         
myframes <- myframes %>%
  left_join(myconditions,by=c('file_path')) %>% 
  group_by(file_path,condition) %>% 
  filter(data.table::between(frame,f_start,f_end-1)) %>% 
  tidyr::extract(file_path, c("date","pos", "gl"), ".*(\\d{8})_.*[Pp]os(\\d+).*_GL(\\d+).*", remove=FALSE, convert=TRUE) %>%
  ungroup()

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
  #Set cell ref
  mutate(cell=paste(date,pos,orientation,gl,id,sep='.'),
         parent=paste(date,pos,orientation,gl,parent_id,sep='.'),
         lane_ID=paste(date,pos,orientation,gl,sep='.')) %>% 
  
  #Compute biomass estimations
  mutate(#length_um=length_pixel*dl, #ellipse fitting
         #width_um=width_pixel*dl,
         #surface=pi*width_um*length_um/4,
         
         #length_um_prob=area_probability_sum__px*(dl**2),
         
         #length_um_bbox=(vertical_bottom-vertical_top)*dl, #bbox
         
         length_um_oriented_bbox=oriented_bbox_length_px*dl, #oriented_bbox
         width_um_oriented_bbox=oriented_bbox_width_px*dl,
         
         #length_um_fake=length_um_oriented_bbox*1.1,
         
         #length_um_spine=spine_length__px*dl, #spine
         #length_um_medial_line=medial_line_length__px*dl, #medial_line
         
         #volume_um=compute_vol(length_um,width_um) # compute volume assuming a perfect rod shape
         
         ) %>%
  #Compute vertical center
  mutate(vertical_center=(vertical_bottom + vertical_top)/2) %>%
  # Propagating end_type information
  mutate(gl_id=paste(date, pos, orientation , gl, sep='.')) %>% 
  group_by(cell) %>% 
  fill(end_type,.direction="up") %>%
  mutate(discard_top=(vertical_top < vertical_cutoff)) %>%
  mutate(end_type=ifelse(any(discard_top), 'touchtop', end_type)) %>%
  ungroup()

myframes <- myframes %>% 
  #filter(!(cell %in% wrong_cells)) %>%
  mutate(time_sec=t_start+(frame-f_start)*t_interval*60) %>% 
  group_by(cell) %>% 
  partition(cluster=mycluster) %>%
  mutate(cell_rank_beg=first(cell_rank)) %>% 
  mutate(start_time=first(time_sec), end_time=last(time_sec),
         #cell_cycle=ifelse(end_type=='div', ((1:length(frame))-.5) / length(frame), NA),
         #length_predict=fit_exp_elongation(time_sec,length_um), #,
         
         #length_predict=fit_exp_elongation_predict(time_sec, length_um),
         #alpha=fit_exp_elongation_slope(time_sec,length_um),
         #sd_alpha=fit_exp_elongation_sd_slope(time_sec,length_um),
         
         #alpha_surface=fit_exp_elongation_slope(time_sec,surface), 
         #sd_alpha_surface=fit_exp_elongation_sd_slope(time_sec,surface),
         
         #length_predict_bbox=fit_exp_elongation_predict(time_sec, length_um_bbox),
         #alpha_bbox=fit_exp_elongation_slope(time_sec,length_um_bbox),
         #sd_alpha_bbox=fit_exp_elongation_sd_slope(time_sec,length_um_bbox),
         
         length_predict_oriented_bbox=fit_exp_elongation_predict(time_sec, length_um_oriented_bbox),
         alpha_oriented_bbox=fit_exp_elongation_slope(time_sec,length_um_oriented_bbox),
         sd_alpha_oriented_bbox=fit_exp_elongation_sd_slope(time_sec,length_um_oriented_bbox),
         
         #length_predict_spine=fit_exp_elongation_predict(time_sec, length_um_spine),
         #alpha_spine=fit_exp_elongation_slope(time_sec,length_um_spine),
         #sd_alpha_spine=fit_exp_elongation_sd_slope(time_sec,length_um_spine),
         
         #length_predict_medial_line=fit_exp_elongation_predict(time_sec, length_um_medial_line),
         #alpha_medial_line=fit_exp_elongation_slope(time_sec,length_um_medial_line),
         #sd_alpha_medial_line=fit_exp_elongation_sd_slope(time_sec,length_um_medial_line),
         
         #length_predict_fake=fit_exp_elongation_predict(time_sec, length_um_fake),
         #alpha_fake=fit_exp_elongation_slope(time_sec,length_um_fake),
         #sd_alpha_fake=fit_exp_elongation_sd_slope(time_sec,length_um_fake),
         
         #length_predict_prob=fit_exp_elongation_predict(time_sec, length_um_prob),
         #alpha_prob=fit_exp_elongation_slope(time_sec,length_um_prob),
         #sd_alpha_prob=fit_exp_elongation_sd_slope(time_sec,length_um_prob)
         

         ) %>%
  collect() %>% 
  arrange(cell, frame) %>% # sort data after `partition()`
  ungroup()

# Add the mean_position of the cell between two divisions in the growth lane.
myframes <- myframes %>% 
  group_by(cell) %>%
  mutate(mean_position_um=mean(vertical_center)*dl) %>% 
  ungroup()

# COMPUTE INSTANTANEOUS RATES ####
# discrete fluo derivative
#myframes <- myframes %>% 
#  group_by(date, pos, orientation, gl, id) %>% 
#  mutate(gfp_deriv=(gfp_nb-lag(gfp_nb))/(time_sec-lag(time_sec)),
#         l_deriv=(length_um-lag(length_um))/(time_sec-lag(time_sec))) %>% 
#  ungroup
