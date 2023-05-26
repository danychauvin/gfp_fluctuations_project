datatypes <- c("allnoisesources","nodivisionnoise","nogrowthnoise","noprodnoise","nonoise")
replicas <- c("1","2")
noise_contribution_data <- noise_contribution_data %>% 
  mutate(replica=ifelse(replica==3,1,replica))

data_toplot <- noise_contribution_data %>% 
  filter(!((condition=="glucoseaa020")&(promoter=="med3"))) %>%
  filter(datatype %in% datatypes) %>% 
  mutate(x=concentration) %>% 
  group_by(condition,promoter,replica,date,datatype) %>%
  mutate(q25=quantile(x,0.25)) %>% 
  mutate(q75=quantile(x,0.75)) %>%
  mutate(q50=quantile(x,0.50)) %>%
  mutate(mean_x=mean(x)) %>% 
  mutate(ir=log(q75)-log(q25)) %>% 
  ungroup() %>% 
  distinct(condition,promoter,replica,date,datatype,ir)

data_toplot %>% 
  left_join(condition_label_df,by=c("condition")) %>% 
  left_join(datatype_label_df,by=c("datatype")) %>% 
  ggplot()+
  geom_point(aes(factor(datatype_label,levels=datatypes_labels),ir,col=promoter),alpha=0.5)+
  facet_wrap(~factor(condition_label,levels=conditions_labels),scale="free")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  expand_limits(y=0)+
  xlab("")+
  ylab("Interquartile range")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/IR_vs_parameters_facet_condition.pdf",width=12,height=6)

data_toplot %>% 
  left_join(condition_label_df,by=c("condition")) %>% 
  ggplot()+
  geom_point(aes(factor(datatype,levels=datatypes),ir,col=promoter),alpha=0.5)+
  geom_line(aes(factor(datatype,levels=datatypes),ir,col=promoter,group=interaction(condition,promoter,replica)),alpha=0.5)+
  facet_wrap(~factor(condition_label,levels=conditions_labels)+promoter,scale="free",ncol=6)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  expand_limits(y=0)+
  xlab("")+
  ylab("Interquartile range")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/IR_vs_parameters_facet_promoter_condition.pdf",width=12,height=8)
