datatypes <- c("allnoisesources","nodivisionnoise","nogrowthnoise","noprodnoise","nonoise")
replicas <- c("1","2")
#Correct for mistake with rrnB third replica
experimental_data <- experimental_data %>% 
  mutate(replica=ifelse(replica==3,1,replica))

data_toplot <- experimental_data %>% 
  #filter(promoter=="med2") %>% 
  filter(!((condition=="glucoseaa020")&(promoter=="med3"))) %>%
  filter(promoter!="6300") %>% 
  #mutate(x=log(gfp_nb/length_um)) %>%
  mutate(x=gfp_nb) %>% 
  #computing error bars for log concentration
  #mutate(sd_x=sqrt(cov_gg/((gfp_nb)**2)+cov_xx-2*cov_xg/gfp_nb)) %>% 
  mutate(sd_x=sqrt(abs(cov_gg))) %>% 
  group_by(condition,promoter,replica,date) %>%
  #sample_n(500) %>%
  mutate(mean_gr=compute_weighted_mean(lambda,sqrt(abs(cov_ll)))) %>% 
  mutate(sd_gr=compute_weighted_sd(lambda,sqrt(abs(cov_ll)))) %>% 
  mutate(q25=quantile(x,0.25)) %>% 
  mutate(q75=quantile(x,0.75)) %>%
  mutate(q50=quantile(x,0.50)) %>%
  mutate(mean_x=compute_weighted_mean(x,sd_x)) %>% 
  mutate(sd_x=compute_weighted_sd(x,sd_x)) %>% 
  mutate(ir=q75-q25) %>% 
  ungroup() %>% 
  distinct(condition,promoter,replica,date,mean_x,sd_x,mean_gr,sd_gr) %>% 
  rename(mean_gfp=mean_x,sd_gfp=sd_x)

readr::write_csv(data_toplot,"/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/mean_tot_GFP_vs_mean_growth_rate_table_all.csv")