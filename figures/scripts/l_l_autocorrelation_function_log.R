gamma_lambda_df <- myparameters %>% 
  filter(name=="gamma_lambda") %>% 
  select(promoter,condition,date,final) %>% 
  #mutate(final=final) %>% 
  rename(gamma_lambda=final) %>% 
  left_join(doubling_times_df,by=c("promoter","condition","date")) %>% 
  mutate(gamma_lambda_normalized=gamma_lambda*(mean_doubling_time*60)) %>% #minutes
  select(condition,promoter,date,gamma_lambda,gamma_lambda_normalized)

gamma_lambda_err_df <- myerrors %>% 
  filter(epsilon==0.05) %>% 
  mutate(gamma_lambda_err=sqrt(gamma_lambda)) %>% 
  select(condition,promoter,date,gamma_lambda_err)

correlations_df %>%
  left_join(doubling_times_df,by=c("promoter","condition","date")) %>% 
  left_join(gamma_lambda_df,c("promoter","condition","date")) %>% 
  left_join(gamma_lambda_err_df,by=c("promoter","condition","date")) %>% 
  mutate(condition=factor(condition,levels=c("acetate005","glycerol040","glucose020","glucoseaa020"))) %>% 
  mutate(dt_norm=dt/(mean_doubling_time*60)) %>% 
  mutate(corr_OU=exp(-gamma_lambda*dt),
         corr_OU_min=exp(-(gamma_lambda+gamma_lambda_err)*dt),
         corr_OU_max=exp(-(gamma_lambda-gamma_lambda_err)*dt)) %>%
  filter(dt_norm<3) %>% 
  ggplot()+
  geom_line(aes(dt_norm,corr_l_l,col=condition,group=interaction(promoter,condition,date)))+
  geom_ribbon(aes(x=dt_norm,y=corr_l_l,ymin=corr_l_l-corr_l_l_err,ymax=corr_l_l+corr_l_l_err,fill=condition,group=interaction(promoter,condition,date)),alpha=0.2)+
  geom_line(aes(dt_norm,corr_OU,col=condition,group=interaction(promoter,condition,date)),linetype="dashed")+
  geom_ribbon(aes(x=dt_norm,y=corr_OU,ymin=corr_OU_min,ymax=corr_OU_max,fill=condition,group=interaction(promoter,condition,date)),alpha=0.2)+
  theme_cowplot()+
  facet_wrap(~promoter)+
  #scale_y_continuous(trans="log10",limits=c(0.1,1))+
  #annotation_logticks(sides = "l") + 
  xlim(c(0,1))+
  #ylim(c(0,1))+
  #coord_cartesian(ylim=c(0.1115,1))+
  xlab("dt / mean doubing time")+
  ylab("Correlation l(t+dt),l(t)")

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/fig2/l_l_autocorrelation_function_linear.pdf",width=14,height=4)

correlations_df %>%
  left_join(doubling_times_df,by=c("promoter","condition","date")) %>% 
  left_join(gamma_lambda_df,c("promoter","condition","date")) %>% 
  left_join(gamma_lambda_err_df,by=c("promoter","condition","date")) %>% 
  mutate(condition=factor(condition,levels=c("acetate005","glycerol040","glucose020","glucoseaa020"))) %>% 
  mutate(dt_norm=dt/(mean_doubling_time*60)) %>% 
  mutate(corr_OU=exp(-gamma_lambda*dt),
         corr_OU_min=exp(-(gamma_lambda+gamma_lambda_err)*dt),
         corr_OU_max=exp(-(gamma_lambda-gamma_lambda_err)*dt)) %>%
  filter(dt_norm<3) %>% 
  ggplot()+
  geom_line(aes(dt_norm,corr_l_l,col=condition,group=interaction(promoter,condition,date)))+
  geom_ribbon(aes(x=dt_norm,y=corr_l_l,ymin=corr_l_l-corr_l_l_err,ymax=corr_l_l+corr_l_l_err,fill=condition,group=interaction(promoter,condition,date)),alpha=0.2)+
  geom_line(aes(dt_norm,corr_OU,col=condition,group=interaction(promoter,condition,date)),linetype="dashed")+
  geom_ribbon(aes(x=dt_norm,y=corr_OU,ymin=corr_OU_min,ymax=corr_OU_max,fill=condition,group=interaction(promoter,condition,date)),alpha=0.2)+
  theme_cowplot()+
  facet_wrap(~promoter)+
  scale_y_continuous(trans="log10",limits=c(0.1,1))+
  annotation_logticks(sides = "l") + 
  xlim(c(0,1))+
  #ylim(c(0,1))+
  coord_cartesian(ylim=c(0.1115,1))+
  xlab("dt / mean doubing time")+
  ylab("Correlation l(t+dt),l(t)")

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/fig2/l_l_autocorrelation_function_log.pdf",width=14,height=4)