#We want to plot, a heat-map of actual growth-rate versus production (mean).For each condition for hi1.
#We assume here that data were already imported. Using load_simulation_data.R

datatypes <- c("allnoisesources","nodivisionnoise","nogrowthnoise","noprodnoise","nonoise")
replicas <- c("1","2")
#Correct for mistake with rrnB third replica
experimental_data <- experimental_data %>% 
  mutate(replica=ifelse(replica==3,1,replica))

toplot <- experimental_data %>% 
  filter(!((condition=="glucoseaa020")&(promoter=="med3"))) %>%
  filter(promoter!="6300") %>% 
  mutate(replica=as.character(replica)) %>% 
  filter(q-sqrt(abs(cov_qq))>0) %>% 
  mutate(l_q=log(q))
  #filter(q-sqrt(abs(cov_qq))>0) %>% 
  #mutate(l_q=log(q))

toplot <- toplot %>% 
  group_by(condition,promoter,date,replica) %>% 
  sample_n(200) %>%
  mutate(corr=cor(lambda,l_q)) %>% 
  mutate(mean_gr=compute_weighted_mean(lambda,sqrt(abs(cov_ll)))) %>% 
  mutate(sd_gr=compute_weighted_sd(lambda,sqrt(abs(cov_ll)))) %>%
  ungroup() %>% 
  mutate(date=as.character(date)) %>%
  mutate(replica=as.character(replica)) %>% 
  distinct(condition,promoter,replica,date,mean_gr,corr)

toplot %>% 
  ggplot()+
  geom_point(aes(mean_gr,corr,col=date))+
  facet_wrap(~factor(promoter,levels=promoters))+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Mean growth rate (min-1)")+
  ylab("Pearson correlation coefficient")+
  #scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(limits=c(0,0.025))+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))+
  scale_colour_manual(values=c25)

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/fig2/corr_log_q_lambda_vs_mean_lambda_S200.pdf",width=6,height=4)
  