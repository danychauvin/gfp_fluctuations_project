# We assume raw data are already loaded

raw_data %>% 
  left_join(condition_label_df,by=c("condition")) %>% 
  mutate(condition=factor(condition,levels=conditions)) %>% 
  mutate(promoter=factor(promoter,levels=promoters)) %>% 
  mutate(conditions_label=factor(condition_label,levels=conditions_labels)) %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  arrange(condition_label) %>% 
  ggplot(aes(x=alpha_oriented_bbox*60))+
  geom_freqpoly(aes(y=after_stat(density),col=as.character(replica)),bins=40)+
  theme_cowplot()+
  scale_y_continuous(trans="log10")+
  facet_wrap(~fct_relevel(condition_label,levels=conditions_labels)+fct_relevel(promoter,levels=promoters),scale="free",ncol=7)+
  scale_color_manual(values=c25)+
  ylab("Density")+
  xlab("Growth-rate (min-1)")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))+
  labs(col="Replica") 
  
ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/growth_rate_distributions_full_cell_cycle_raw_data_Fconditionpromoter_Creplica.pdf",width=24,height=16,units="cm")