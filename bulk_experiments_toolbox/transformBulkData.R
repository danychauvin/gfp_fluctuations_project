# This script imports BULK results and performs the necessary data transformation
# importMMData.R has to be run before running importingBulkData.R
# IMPORTING DATA FROM BULK EXPERIMENTS

mydata <- mydata %>% 
  mutate(valid=NA) %>% 
  mutate(well=paste(plate,row,column,strain,sep="_"))

mean_blank_od <- 0.03840125 #Determined in a previous experiment where we looked ad OD of empty wells for entire plates
mydata <- mydata %>% 
  mutate(corrected_od=od-mean_blank_od)

mean_blank_fluo <- 35.3075 #Idem
mydata <- mydata %>% 
  mutate(corrected_fluo=fluo-mean_blank_fluo) %>% 
  filter(!grepl("H_7",well)) #Discarding this well which is problematic

#path_to_manual_curation<-"./bulk_experiments_toolbox/manualCuration.csv"
mydata_bad <- readr::read_csv(path_to_manual_curation,
                              show_col_types = FALSE)

condition_df <- tibble(old_condition=c("acetate","glycerol","glucose","glucoseaa"),condition=c("acetate020","glycerol040","glucose020","glucoseaa020"))
conditions <- c("acetate020","glycerol040","glucose020","glucoseaa020")

mydata <- mydata %>% 
  mutate(rowcol=paste(row,column,sep="")) %>% 
  rename(old_condition=condition) %>% 
  left_join(condition_df,by=c("old_condition"))
  
replica_df <- tibble(description=c("20210312constitexpr3.0","20210316constitexpr3.1","20210204Testrun2",
                                   "20210430constitexpr3.2",
                                   "20210504constitexpr3.3",
                                   "20210507constitexpr3.4",
                                   "20210512constitexpr3.5"),
                     replica=c(2,3,1,4,5,6,7))
promoters_of_interest <- list("hi1","hi3","med2","med3","rpsB","rpmB","rplN","rrnB")

mydata_curated <-mydata %>% 
  anti_join(mydata_bad,by=c("description","plate","rowcol")) %>% 
  left_join(replica_df,by=c("description")) %>% 
  mutate(well_rep=paste(well,replica,sep="-"))

mydata_curated <- rbind(

  mydata_curated %>% 
    filter((strain %in% c("MG1655","empty"))) %>% 
    mutate(promoter="nopromoter") %>% 
    mutate(vector="novector"),
  
  mydata_curated %>% 
    filter(!(strain %in% c("MG1655","empty"))) %>% 
    group_by(strain) %>% 
    mutate(prom_part1=stringr::str_split(strain,"-")[[1]][[1]],
           prom_part2=stringr::str_split(strain,"-")[[1]][[2]]) %>% 
    mutate(promoter=ifelse(prom_part2 %in% promoters_of_interest,prom_part2,paste(prom_part1,prom_part2,sep=""))) %>% 
    mutate(vector=ifelse(prom_part2 %in% promoters_of_interest,prom_part1,"pl")) %>%
    select(-c(prom_part1,prom_part2,strain)) %>% 
    mutate(strain=paste(vector,promoter,sep="-")) %>% 
    ungroup())

# INFERRING SLOPES AND INTERCEPT OF FLUO VS OD

od_range_df <- tibble(condition=c("acetate020","glycerol040","glucose020","glucoseaa020"),
                      min_corrected_od=c(0.05,0.05,0.05,0.05),
                      max_corrected_od=c(0.15,0.2,0.4,0.6))

hi_or_med_p2 <- tibble(plate=rep(c("pl1","pl3","pl5","pl7"),each=6),strength=rep(rep(c("med","hi"),each=3),times=4),column=rep(c(7,8,9,10,11,12),times=4))
hi_or_med_p3 <- tibble(plate=rep(c("pl2","pl4","pl6","pl8"),each=12),strength=rep(rep(c("med","hi"),each=3),times=8),column=rep(c(1,2,3,4,5,6,7,8,9,10,11,12),times=4))
hi_or_med <- rbind(hi_or_med_p2,hi_or_med_p3)

mydata_summary_od <- mydata_curated %>% 
    ungroup() %>%
    left_join(hi_or_med,by=c("plate","column")) %>% 
    filter(!grepl("empty",strain)) %>% #filter out empty wells
    left_join(od_range_df,by=c("condition")) %>%
    filter(max(corrected_od)>min_corrected_od) %>% #to get rid of wells in which there are no growth at all
    group_by(description,plate,well) %>%
    arrange(time_min) %>% 
    mutate(above_max_time=ifelse(corrected_od>=max_corrected_od,time_min,NA)) %>% 
    mutate(time_limit=min(c(above_max_time,Inf),na.rm=TRUE)) %>% 
    filter(time_min<time_limit) %>% 
    mutate(below_min_time=ifelse(corrected_od<=min_corrected_od,time_min,NA)) %>% 
    mutate(time_limit=max(c(below_min_time,-Inf),na.rm=TRUE)) %>% 
    filter(time_min>time_limit) %>% 
    ungroup() %>% 
    select(-c(above_max_time,time_limit,below_min_time)) %>%
    group_by(description,well) %>%
    mutate(p=n(),
           mean_od=mean(corrected_od),
           var_od=1/p*sum((corrected_od-mean_od)**2),
           mean_fluo=mean(fluo),
           var_fluo=1/p*sum((fluo-mean_fluo)**2),
           covar_od_fluo=1/p*sum(corrected_od*fluo-mean_od*mean_fluo),
           alpha_tot=(var_fluo-var_od)/(2*covar_od_fluo)+sqrt(1+((var_fluo-var_od)/(2*covar_od_fluo))**2),#assuming that this is the best slope
           sd_alpha=compute_error(alpha_tot,var_od,var_fluo,covar_od_fluo,p),
           beta=1/p*sum(fluo-alpha_tot*corrected_od))%>%
    ungroup() %>% 
    distinct(strain, replica, well,.keep_all=TRUE) %>% 
  select(condition,strain,promoter,vector,replica,well,alpha_tot,sd_alpha,beta,strength,plate)

# COMPUTING GROWTH-RATES IN SINGLE WELLS

mydata_growthrates_summary <- mydata_curated %>% 
    ungroup() %>%
    left_join(hi_or_med,by=c("plate","column")) %>% 
    filter(!grepl("empty",strain)) %>% #filter out empty wells
    left_join(od_range_df,by=c("condition")) %>%
    filter(max(corrected_od)>min_corrected_od) %>% #to get rid of wells in which there are no growth at all
    group_by(description,plate,well) %>%
    arrange(time_min) %>% 
    mutate(above_max_time=ifelse(corrected_od>=max_corrected_od,time_min,NA)) %>% 
    mutate(time_limit=min(c(above_max_time,Inf),na.rm=TRUE)) %>% 
    filter(time_min<time_limit) %>% 
    mutate(below_min_time=ifelse(corrected_od<=min_corrected_od,time_min,NA)) %>% 
    mutate(time_limit=max(c(below_min_time,-Inf),na.rm=TRUE)) %>% 
    filter(time_min>time_limit) %>% 
    ungroup() %>% 
    select(-c(above_max_time,time_limit,below_min_time)) %>%
    filter(corrected_od>0) %>% 
    group_by(description,well) %>%
    mutate(
      well_gr=fit_exp_elongation_slope(time_min,corrected_od)*60/log(2), 
      sd_well_gr=fit_exp_elongation_sd_slope(time_min,corrected_od)*60/log(2)) %>%
  ungroup() %>% 
  distinct(strain, replica, well,.keep_all=TRUE) %>% 
  select(condition,strain,promoter,vector,replica,well,well_gr,sd_well_gr,strength,plate)

mydata_growthrates_summary %>% 
  distinct(condition,replica,well,.keep_all = TRUE) %>% 
  filter(!((condition=="acetate020") & (replica %in% c(5,6,7)))) %>% 
  filter(!((condition=="glucoseaa020") & (replica %in% c(4)))) %>% 
  group_by(condition) %>% 
  mutate(p=n()) %>% 
  mutate(mean_doubling=mean(well_gr,na.rm=TRUE)) %>% 
  mutate(sd_doubling=sd(well_gr,na.rm=TRUE)) %>% 
  mutate(mean_error_fit=mean(sd_well_gr,na.rm=TRUE)) %>% 
  ungroup() %>% 
  distinct(condition,.keep_all = TRUE) %>% 
  select(condition,mean_doubling,sd_doubling,mean_error_fit) %>% 
  arrange(mean_doubling,condition)

#Necessary for comparison with MM data
#myconcentrations_condition <- mycells %>%  
#  distinct(condition,strain,mean_concentration_condition,sd_concentration_condition)

#Set strains to select for the comparison, for which I have data both in mother machine and spectro.
#selected_strains_df <- mycells %>% 
#  distinct(strain,condition)

mydata_summary_od$replica <- factor(mydata_summary_od$replica,levels = c(1,2,3,4,5,6,7))
mydata_growthrates_summary$replica <- factor(mydata_growthrates_summary$replica,levels = c(1,2,3,4,5,6,7))

mydata_curated_cond_for_plot <- mydata_summary_od %>% 
  distinct(replica,condition,well,alpha_tot,sd_alpha,strain,promoter,vector,plate,strength)

mydata_curated_gr_for_plot <- mydata_growthrates_summary %>% 
  distinct(replica,condition,well,well_gr,sd_well_gr,strain,promoter,vector,plate,strength)

exp_vs_gr_results_df <- left_join(mydata_curated_cond_for_plot,mydata_curated_gr_for_plot ,by=c("replica","condition","well","strain","promoter","vector","plate","strength")) %>% filter(strain=="MG1655") %>% 
  select(-strength) %>% 
  filter(alpha_tot>0) %>% 
  filter(well_gr>0.2) %>% 
  filter(!((condition=="acetate020") & (alpha_tot<1000))) %>% 
  mutate(l_gr=log(well_gr)) %>% #here growth rate is in doublings per hour, which is more intuitive
  mutate(l_exp=log(alpha_tot)) %>% 
  drop_na() %>% 
  mutate(p=n(),
         m_l_gr=mean(l_gr,na.rm=TRUE),
         m_l_exp=mean(l_exp,na.rm=TRUE),
         var_l_gr=var(l_gr,na.rm=TRUE),
         var_l_exp=var(l_exp,na.rm=TRUE),
         cov_l_gr_l_exp=1/p*sum(l_exp*l_gr-m_l_gr*m_l_exp,na.rm=TRUE)) %>% 
  filter(row_number()==1) %>% 
  select(m_l_gr,m_l_exp,var_l_gr,var_l_exp,cov_l_gr_l_exp,p)

exp_vs_gr_results_df

mydata_gr_exp_summary <- left_join(mydata_curated_cond_for_plot,mydata_curated_gr_for_plot ,by=c("replica","condition","well","strain","promoter","vector","plate","strength")) %>%
  filter(!(strain %in% c("MG1655","empty"))) %>% 
  mutate(l_x=log(well_gr),
         var_l_x=((sd_well_gr)**2)/((well_gr)**2)) %>%
  mutate(m_l_gr=unique(exp_vs_gr_results_df$m_l_gr),
         var_l_gr=unique(exp_vs_gr_results_df$var_l_gr),
         m_l_exp=unique(exp_vs_gr_results_df$m_l_exp),
         var_l_exp=unique(exp_vs_gr_results_df$var_l_exp),
         cov_l_gr_l_exp=unique(exp_vs_gr_results_df$cov_l_gr_l_exp),
         p=unique(exp_vs_gr_results_df$p)) %>% 
  mutate(l_y=m_l_exp+cov_l_gr_l_exp/(var_l_gr)*(l_x-m_l_gr),
         r=cov_l_gr_l_exp/(sqrt(var_l_gr*var_l_exp)),
         var_l_y=var_l_x*(cov_l_gr_l_exp**2/(var_l_gr**2))+var_l_exp*(1-r**2)*(1+(var_l_x+(l_x-m_l_gr)**2)/p))

mydata_gr_exp_summary<- mydata_gr_exp_summary %>%
  mutate(alpha_gfp=alpha_tot-exp(l_y),
         sd_alpha_autofluo=sqrt(l_y**2*var_l_y),
         sd_alpha_gfp=sqrt(sd_alpha_autofluo**2+sd_alpha**2)) %>% 
  filter(alpha_gfp>3*sd_alpha_autofluo)

# Formatting data with more comprehensible column names
strength_to_promoter_type <- tibble(strength=c("hi","med"),promoter_type_synth=c("high","medium"))
promoter_type_df <- tibble(promoter=c("hi1","hi3","med2","med3","rpsB","rpmB","rplN","rrnB"),promoter_type=c("high","high","medium","medium","ribosomal","ribosomal","ribsomomal","ribosomal"))


bulk_data <- mydata_gr_exp_summary %>% 
  left_join(strength_to_promoter_type,by=c("strength")) %>% 
  left_join(promoter_type_df,by=c("promoter")) %>% #fixing a mistake regarding the strength of two promoters
  mutate(promoter_type=ifelse(is.na(promoter_type),promoter_type_synth,promoter_type)) %>% 
  mutate(promoter_type=ifelse(promoter %in% c("p2C11","p3D12"),"medium",promoter_type)) %>% 
  select(condition,promoter,vector,replica,well,well_gr,sd_well_gr,alpha_gfp,sd_alpha_gfp,promoter_type)

#Necessary for comparison with MMData
#mycalibrations <- mydata_gr_exp_summary %>%
#  semi_join(selected_strains_df,by=c("strain","condition")) %>% 
#  filter(log(well_gr)>-2) %>% 
#  left_join(myconcentrations_condition,by=c("condition","strain")) %>% 
#  group_by(condition) %>% 
#  mutate(intercept=mean(log(mean_concentration_condition))-mean(log(alpha_gfp))) %>% 
#  ungroup() %>% 
#  distinct(condition,intercept)

#Necessary for comparison with MMData
#mycalibrations$condition <- factor(mycalibrations$condition, levels=c("acetate","glycerol","glucose","glucoseaa"))
#mycalibrations %>% 
#  arrange(condition) %>% 
#  mutate(calibration_constant=exp(intercept))


#mydata_results <- mydata_gr_exp_summary %>% 
#  filter(log(well_gr)>-2) %>% 
#  left_join(mycalibrations,by=c("condition")) %>% 
#  mutate(l_pred_concentration=log(alpha_gfp)+intercept,
#         max_l_pred_concentration=log(alpha_gfp+sd_alpha)+intercept,
#         min_l_pred_concentration=log(alpha_gfp-sd_alpha)+intercept)

#hi_or_med_p2 <- tibble(plate=rep(c("pl1","pl3","pl5","pl7"),each=6),strength=rep(rep(c("med","hi"),each=3),times=4),column=rep(c(7,8,9,10,11,12),times=4))
#hi_or_med_p3 <- tibble(plate=rep(c("pl2","pl4","pl6","pl8"),each=12),strength=rep(rep(c("med","hi"),each=3),times=8),column=rep(c(1,2,3,4,5,6,7,8,9,10,11,12),times=4))
#hi_or_med <- rbind(hi_or_med_p2,hi_or_med_p3)



