# We assume raw data are already loaded

raw_data %>% 
  mutate(concentration=gfp_nb/length_um_oriented_bbox) %>% 
  #filter(concentration>0) %>% 
  left_join(condition_label_df,by=c("condition")) %>% 
  mutate(condition=factor(condition,levels=conditions)) %>% 
  mutate(promoter=factor(promoter,levels=promoters)) %>% 
  mutate(conditions_label=factor(condition_label,levels=conditions_labels)) %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  arrange(condition_label) %>% 
  ggplot(aes(x=concentration))+
  geom_freqpoly(aes(y=after_stat(density),col=as.character(replica)),bins=80)+
  theme_cowplot()+
  #scale_y_continuous(trans="log10")+
  facet_wrap(~fct_relevel(condition_label,levels=conditions_labels)+fct_relevel(promoter,levels=promoters),scale="free",ncol=7)+
  #facet_wrap(~promoter+fct_relevel(condition_label,levels=conditions_labels),scale="free",ncol=4)+
  #facet_grid(fct_relevel(condition_label,levels=conditions_labels)~promoter,scale="free")+
  scale_color_manual(values=c25)+
  ylab("Density")+
  xlab("Concentration (#GFP/um)")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))+
  labs(col='condition') +
  guides(fill = guide_legend(byrow = TRUE))


ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/concentration_distribution_raw_data_Fconditionpromoter.pdf",width=24,height=16,units="cm")

raw_data %>% 
  mutate(concentration=gfp_nb/length_um_oriented_bbox) %>% 
  #filter(concentration>0) %>% 
  left_join(condition_label_df,by=c("condition")) %>% 
  mutate(condition=factor(condition,levels=conditions)) %>% 
  mutate(promoter=factor(promoter,levels=promoters)) %>% 
  mutate(conditions_label=factor(condition_label,levels=conditions_labels)) %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  arrange(condition_label) %>% 
  ggplot(aes(x=concentration))+
  geom_freqpoly(aes(y=after_stat(density),col=as.character(replica)),bins=80)+
  theme_cowplot()+
  scale_y_continuous(trans="log10")+
  #scale_x_continuous(trans="log10")+
  facet_wrap(~fct_relevel(condition_label,levels=conditions_labels)+fct_relevel(promoter,levels=promoters),scale="free",ncol=7)+
  #facet_wrap(~promoter+fct_relevel(condition_label,levels=conditions_labels),scale="free",ncol=4)+
  #facet_grid(fct_relevel(condition_label,levels=conditions_labels)~promoter,scale="free")+
  scale_color_manual(values=c25)+
  ylab("Density")+
  xlab("Concentration (#GFP/um)")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))+
  labs(col='condition') 
guides(fill = guide_legend(byrow = TRUE))


ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/concentration_distribution_raw_data_Fconditionpromoter_logy.pdf",width=24,height=16,units="cm")

raw_data %>% 
  mutate(concentration=gfp_nb/length_um_oriented_bbox) %>% 
  #filter(concentration>0) %>% 
  left_join(condition_label_df,by=c("condition")) %>% 
  mutate(condition=factor(condition,levels=conditions)) %>% 
  mutate(promoter=factor(promoter,levels=promoters)) %>% 
  mutate(conditions_label=factor(condition_label,levels=conditions_labels)) %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  arrange(condition_label) %>% 
  ggplot(aes(x=concentration))+
  geom_freqpoly(aes(y=after_stat(density),col=as.character(replica)),bins=80)+
  theme_cowplot()+
  #scale_y_continuous(trans="log10")+
  scale_x_continuous(trans="log10")+
  facet_wrap(~fct_relevel(condition_label,levels=conditions_labels)+fct_relevel(promoter,levels=promoters),scale="free",ncol=7)+
  #facet_wrap(~promoter+fct_relevel(condition_label,levels=conditions_labels),scale="free",ncol=4)+
  #facet_grid(fct_relevel(condition_label,levels=conditions_labels)~promoter,scale="free")+
  scale_color_manual(values=c25)+
  ylab("Density")+
  xlab("Concentration (#GFP/um)")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))+
  labs(col='condition') 
guides(fill = guide_legend(byrow = TRUE))


ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/concentration_distribution_raw_data_Fconditionpromoter_logx.pdf",width=24,height=16,units="cm")



raw_data %>% 
  mutate(concentration=gfp_nb/length_um_oriented_bbox) %>% 
  #filter(concentration>0) %>% 
  left_join(condition_label_df,by=c("condition")) %>% 
  mutate(condition=factor(condition,levels=conditions)) %>% 
  mutate(promoter=factor(promoter,levels=promoters)) %>% 
  mutate(conditions_label=factor(condition_label,levels=conditions_labels)) %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  arrange(condition_label) %>% 
  ggplot(aes(x=concentration))+
  geom_freqpoly(aes(y=after_stat(density),col=as.character(replica)),bins=80)+
  theme_cowplot()+
  scale_y_continuous(trans="log10")+
  scale_x_continuous(trans="log10")+
  facet_wrap(~fct_relevel(condition_label,levels=conditions_labels)+fct_relevel(promoter,levels=promoters),scale="free",ncol=7)+
  #facet_wrap(~promoter+fct_relevel(condition_label,levels=conditions_labels),scale="free",ncol=4)+
  #facet_grid(fct_relevel(condition_label,levels=conditions_labels)~promoter,scale="free")+
  scale_color_manual(values=c25)+
  ylab("Density")+
  xlab("Concentration (#GFP/um)")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))+
  labs(col='condition') 
guides(fill = guide_legend(byrow = TRUE))


ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/concentration_distribution_raw_data_Fconditionpromoter_logx_logy.pdf",width=24,height=16,units="cm")
