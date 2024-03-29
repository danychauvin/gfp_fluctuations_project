datatypes <- c("allnoisesources","nodivisionnoise","nogrowthnoise","noprodnoise","nonoise")
replicas <- c("1","2")
#Correct for mistake with rrnB third replica
experimental_data <- experimental_data %>% 
  mutate(replica=ifelse(replica==3,1,replica))

data_toplot <- experimental_data %>% 
  filter(!((condition=="glucoseaa020")&(promoter=="med3"))) %>%
  filter(promoter!="6300") %>% 
  filter(q>0) %>% 
  mutate(x=log(q),
         sd_x=sqrt(cov_qq/(q**2))) %>% 
  group_by(condition,promoter,replica,date) %>%
  #sample_n(200) %>%
  mutate(mean_gr=compute_weighted_mean(lambda,sqrt(abs(cov_ll)))) %>% 
  mutate(sd_gr=compute_weighted_sd(lambda,sqrt(abs(cov_ll)))) %>% 
  mutate(q25=quantile(x,0.25),na.rm=TRUE) %>% 
  mutate(q75=quantile(x,0.75),na.rm=TRUE) %>%
  mutate(q50=quantile(x,0.50),na.rm=TRUE) %>%
  mutate(mean_x=compute_weighted_mean(x,sd_x)) %>% 
  mutate(sd_x=compute_weighted_sd(x,sd_x)) %>% 
  mutate(ir=q75-q25) %>% 
  ungroup() %>% 
  distinct(condition,promoter,replica,date,mean_gr,ir,sd_gr,mean_x,sd_x) %>% 
  mutate(date=as.character(date))

data_toplot %>% 
  left_join(label_df,by=c("condition")) %>% 
  ggplot()+
  geom_point(aes(mean_x,ir,col=promoter,shape=factor(condition_label,levels=conditions_labels)),alpha=0.5)+
  #geom_errorbar(aes(x=mean_x,ymax=mean_x+sd_x,ymin=mean_x-sd_x,col=date),alpha=0.5)+
  #geom_errorbar(aes(y=mean_x,xmax=mean_gr+sd_gr,xmin=mean_gr-sd_gr,col=date),alpha=0.5)+
  #facet_wrap(~factor(promoter,levels=promoters))+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Log-production interquartile range")+
  xlab("Mean log-production")+
  expand_limits(y=0)+
  #scale_y_continuous(limits=c(0,0.8))+
  #scale_x_continuous(limits=c(0,0.03))+
  #scale_color_manual(values=c25)+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))+
  labs(shape="condition",col="promoter")

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/fig2/ir_log_production_vs_mean_log_production_nofacet.pdf",width=6,height=4)