# Correlation between l(t) and q(t+dt)

correlations_df %>%
  left_join(doubling_times_df,by=c("condition","promoter","date")) %>%
  mutate(condition=factor(condition,levels=c("acetate005","acetate020","glycerol040","glucose020","glucoseaa020"))) %>%
  mutate(dt_norm=dt/(mean_doubling_time*60)) %>%
  filter(dt_norm<3) %>%
  ggplot()+
  geom_line(aes(dt_norm,corr_q_l,col=condition,group=interaction(promoter,condition,date)))+
  geom_ribbon(aes(x=dt_norm,y=corr_q_l,ymin=corr_q_l-corr_q_l_err,ymax=corr_q_l+corr_q_l_err,fill=condition,group=interaction(promoter,condition,date)))+
  #geom_line(aes(dt_norm,corr_OU,col=condition),linetype="dashed")+
  #geom_ribbon(aes(x=dt_norm,y=corr_OU,ymin=corr_OU_min,ymax=corr_OU_max,fill=condition),alpha=0.2)+
  theme_cowplot()+
  facet_wrap(~promoter)+
  #scale_y_continuous(trans="log10",limits=c(0.1,1))+
  #annotation_logticks(sides = "l") +
  #xlim(c(0,1))+
  #ylim(c(0,1))+
  #coord_cartesian(ylim=c(0.1115,1))+
  xlab("dt / mean doubing time")+
  ylab("Correlation q(t+dt),l(t)")

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/fig2/l_q_autocorrelation_function_linear.pdf",width=14,height=4)

# Correlation between q(t) and l(t+dt)

#Import parameters first
correlations_df %>%
  left_join(doubling_times_df,by=c("condition","promoter","date")) %>%
  mutate(condition=factor(condition,levels=c("acetate005","acetate020","glycerol040","glucose020","glucoseaa020"))) %>%
  mutate(dt_norm=dt/(mean_doubling_time*60)) %>%
  filter(dt_norm<3) %>%
  ggplot()+
  geom_line(aes(dt_norm,corr_l_q,col=condition,group=interaction(promoter,condition,date)))+
  geom_ribbon(aes(x=dt_norm,y=corr_l_q,ymin=corr_l_q-corr_l_q_err,ymax=corr_l_q+corr_l_q_err,fill=condition,group=interaction(promoter,condition,date)),alpha=0.2)+
  #geom_line(aes(dt_norm,corr_OU,col=condition),linetype="dashed")+
  #geom_ribbon(aes(x=dt_norm,y=corr_OU,ymin=corr_OU_min,ymax=corr_OU_max,fill=condition),alpha=0.2)+
  theme_cowplot()+
  facet_wrap(~promoter)+
  #scale_y_continuous(trans="log10",limits=c(0.1,1))+
  #annotation_logticks(sides = "l") +
  #xlim(c(0,1))+
  #ylim(c(0,1))+
  #coord_cartesian(ylim=c(0.1115,1))+
  xlab("dt / mean doubing time")+
  ylab("Correlation l(t+dt),q(t)")

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/fig2/q_l_autocorrelation_function_linear.pdf",width=14,height=4)