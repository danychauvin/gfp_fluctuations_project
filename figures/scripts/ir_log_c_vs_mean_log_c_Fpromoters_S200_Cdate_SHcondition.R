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
  mutate(date=as.character(date)) %>% 
  group_by(condition,promoter,replica,date) %>%
  sample_n(200) %>%
  mutate(mean_gr=compute_weighted_mean(lambda,sqrt(abs(cov_ll)))) %>% 
  mutate(sd_gr=compute_weighted_sd(lambda,sqrt(abs(cov_ll)))) %>% 
  mutate(q25=quantile(x,0.25)) %>% 
  mutate(q75=quantile(x,0.75)) %>%
  mutate(q50=quantile(x,0.50)) %>%
  mutate(mean_x=compute_weighted_mean(x,sd_x)) %>% 
  mutate(sd_x=compute_weighted_sd(x,sd_x)) %>% 
  mutate(ir=q75-q25) %>% 
  ungroup() %>% 
  distinct(condition,promoter,replica,date,ir,mean_x,sd_x)

data_toplot %>% 
  left_join(label_df,by=c("condition")) %>% 
  ggplot()+
  geom_point(aes(mean_x,ir,col=date,shape=factor(condition_label,levels=conditions_labels)),alpha=0.5)+
  facet_wrap(~factor(promoter,levels=promoters),scale="free")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Mean log-concentration")+
  ylab("Interquartile range")+
  scale_y_continuous(limits=c(0,0.8))+
  scale_color_manual(values=c25)+
  #scale_x_continuous(limits=c(0,0.025))+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))+
  labs(shape="condition",col="date")

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/fig2/interquartile_log_c_vs_mean_log_c_Fpromoters_S200_Cdate_SHcondition.pdf",width=6,height=4)