#Import parameters first
correlations_df %>%
  left_join(doubling_times_df,by=c("condition","promoter","date")) %>%
  mutate(condition=factor(condition,levels=c("acetate005","acetate020","glycerol040","glucose020","glucoseaa020"))) %>%
  mutate(dt_norm=dt/(mean_doubling_time*60)) %>%
  filter(dt_norm<4) %>%
  ggplot()+
  geom_line(aes(dt_norm,corr_c_c,col=condition,group=interaction(condition,promoter,date)))+
  geom_ribbon(aes(x=dt_norm,y=corr_c_c,ymin=corr_c_c-corr_c_c_err,ymax=corr_c_c+corr_c_c_err,fill=condition,group=interaction(condition,promoter,date)),alpha=0.2)+
  theme_cowplot()+
  facet_wrap(~promoter)+
  ylim(c(0,1))+
  xlab("dt / mean doubing time")+
  ylab("Correlation c(t+dt),c(t)")

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/fig2/c_c_autocorrelation_function.pdf",width=14,height=4)