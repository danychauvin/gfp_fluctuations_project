#Import parameters first
gamma_q_df <- myparameters %>%
  filter(name=="gamma_q") %>%
  select(promoter,condition,date,final) %>%
  #mutate(final=final) %>%
  rename(gamma_q=final) %>%
  left_join(doubling_times_df,by=c("promoter","date","condition")) %>%
  mutate(gamma_q_normalized=gamma_q*(mean_doubling_time*60)) %>%
  select(promoter,condition,date,gamma_q,gamma_q_normalized)

gamma_q_err_df <- myerrors %>%
  filter(epsilon==0.05) %>%
  mutate(gamma_q_err=sqrt(gamma_q)) %>%
  select(gamma_q_err,condition,promoter,date)

correlations_df %>%
  left_join(doubling_times_df,by=c("condition","date","promoter")) %>%
  left_join(gamma_q_df,by=c("condition","date","promoter")) %>%
  left_join(gamma_q_err_df,by=c("condition","date","promoter")) %>%
  mutate(condition=factor(condition,levels=c("acetate005","glycerol040","glucose020","glucoseaa020"))) %>%
  mutate(dt_norm=dt/(mean_doubling_time*60)) %>%
  mutate(corr_OU=exp(-gamma_q*dt),
         corr_OU_min=exp(-(gamma_q+gamma_q_err)*dt),
         corr_OU_max=exp(-(gamma_q-gamma_q_err)*dt)) %>%
  filter(dt_norm<1) %>%
  ggplot()+
  geom_line(aes(dt_norm,corr_q_q,col=condition,group=interaction(promoter,condition,date)))+
  geom_ribbon(aes(x=dt_norm,y=corr_q_q,ymin=corr_q_q-corr_q_q_err,ymax=corr_q_q+corr_q_q_err,fill=condition,group=interaction(promoter,condition,date)),alpha=0.2)+
  geom_line(aes(dt_norm,corr_OU,col=condition,group=interaction(promoter,condition,date)),linetype="dashed")+
  geom_ribbon(aes(x=dt_norm,y=corr_OU,ymin=corr_OU_min,ymax=corr_OU_max,fill=condition,group=interaction(promoter,condition,date)),alpha=0.2)+
  theme_cowplot()+
  scale_y_continuous(trans="log10",limits=c(0.1,1))+
  annotation_logticks(sides = "l") +
  facet_wrap(~promoter)+
  #xlim(c(0,1))+
  #ylim(c(0,1))+
  #coord_cartesian(ylim=c(0.1115,1))+
  xlab("dt / mean doubing time")+
  ylab("Correlation q(t+dt),q(t)")

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/fig2/q_q_autocorrelation_function_log.pdf",width=14,height=4)


correlations_df %>%
  left_join(doubling_times_df,by=c("condition","date","promoter")) %>%
  left_join(gamma_q_df,by=c("condition","date","promoter")) %>%
  left_join(gamma_q_err_df,by=c("condition","date","promoter")) %>%
  mutate(condition=factor(condition,levels=c("acetate005","glycerol040","glucose020","glucoseaa020"))) %>%
  mutate(dt_norm=dt/(mean_doubling_time*60)) %>%
  mutate(corr_OU=exp(-gamma_q*dt),
         corr_OU_min=exp(-(gamma_q+gamma_q_err)*dt),
         corr_OU_max=exp(-(gamma_q-gamma_q_err)*dt)) %>%
  filter(dt_norm<1) %>%
  ggplot()+
  geom_line(aes(dt_norm,corr_q_q,col=condition,group=interaction(promoter,condition,date)))+
  geom_ribbon(aes(x=dt_norm,y=corr_q_q,ymin=corr_q_q-corr_q_q_err,ymax=corr_q_q+corr_q_q_err,fill=condition,group=interaction(promoter,condition,date)),alpha=0.2)+
  geom_line(aes(dt_norm,corr_OU,col=condition,group=interaction(promoter,condition,date)),linetype="dashed")+
  geom_ribbon(aes(x=dt_norm,y=corr_OU,ymin=corr_OU_min,ymax=corr_OU_max,fill=condition,group=interaction(promoter,condition,date)),alpha=0.2)+
  theme_cowplot()+
  #scale_y_continuous(trans="log10",limits=c(0.1,1))+
  annotation_logticks(sides = "l") +
  facet_wrap(~promoter)+
  #xlim(c(0,1))+
  #ylim(c(0,1))+
  #coord_cartesian(ylim=c(0.1115,1))+
  xlab("dt / mean doubing time")+
  ylab("Correlation q(t+dt),q(t)")

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/fig2/q_q_autocorrelation_function_linear.pdf",width=14,height=4)
