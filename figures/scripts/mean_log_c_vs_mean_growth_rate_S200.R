datatypes <- c("allnoisesources","nodivisionnoise","nogrowthnoise","noprodnoise","nonoise")
replicas <- c("1","2")
#Correct for mistake with rrnB third replica
experimental_data <- experimental_data %>% 
  mutate(replica=ifelse(replica==3,1,replica))

data_toplot <- experimental_data %>% 
  filter(!((condition=="glucoseaa020")&(promoter=="med3"))) %>%
  filter(promoter!="6300") %>% 
  mutate(x=log(gfp_nb/length_um)) %>% 
  #computing error bars for log concentration
  mutate(sd_x=sqrt(cov_gg/((gfp_nb)**2)+cov_xx-2*cov_xg/gfp_nb)) %>% 
  group_by(condition,promoter,replica,date) %>%
  sample_n(200) %>%
  mutate(mean_gr=compute_weighted_mean(lambda,sqrt(abs(cov_ll)))) %>% 
  mutate(sd_gr=compute_weighted_sd(lambda,sqrt(abs(cov_ll)))) %>% 
  mutate(q25=quantile(x,0.25)) %>% 
  mutate(q75=quantile(x,0.75)) %>%
  mutate(q50=quantile(x,0.50)) %>%
  mutate(mean_x=compute_weighted_mean(x,sd_x)) %>% 
  mutate(sd_x=compute_weighted_sd(x,sd_x)) %>% 
  mutate(ir=log(q75)-log(q25)) %>% 
  ungroup() %>% 
  distinct(condition,promoter,replica,date,mean_x,sd_x,mean_gr,sd_gr)

data_toplot %>% 
  left_join(label_df,by=c("condition")) %>% 
  ggplot()+
  geom_point(aes(mean_gr,mean_x,col=promoter),alpha=0.5)+
  geom_errorbar(aes(x=mean_gr,ymax=mean_x+sd_x,ymin=mean_x-sd_x,col=promoter),alpha=0.5)+
  geom_errorbar(aes(y=mean_x,xmax=mean_gr+sd_gr,xmin=mean_gr-sd_gr,col=promoter),alpha=0.5)+
  #facet_wrap(~factor(promoter,levels=promoters),scale="free")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Mean growth rate (min-1)")+
  ylab("Mean log-concentration")+
  #scale_y_continuous(limits=c(0,0.8))+
  scale_x_continuous(limits=c(0,0.030))+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/fig2/mean_log_c_vs_mean_growth_rate_all_same_plot_sample200.pdf",width=6,height=4)