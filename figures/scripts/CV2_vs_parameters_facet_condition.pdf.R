datatypes <- c("allnoisesources","nodivisionnoise","nogrowthnoise","noprodnoise","nonoise")
replicas <- c("1","2")
noise_contribution_data <- noise_contribution_data %>% 
  mutate(replica=ifelse(replica==3,1,replica))

data_toplot <- noise_contribution_data %>% 
  filter(!((condition=="glucoseaa020")&(promoter=="med3"))) %>%
  filter(datatype %in% datatypes) %>% 
  mutate(x=concentration) %>% 
  mutate(l_x=log(concentration)) %>% 
  #mutate(sd_x=sqrt((cov_gg+gfp_nb**2*cov_xx-2*gfp_nb*cov_xg)/length_um**2)) %>% 
  group_by(condition,promoter,replica,date,datatype) %>%
  mutate(q25=quantile(x,0.25)) %>% 
  mutate(q75=quantile(x,0.75)) %>%
  mutate(q50=quantile(x,0.50)) %>%
  mutate(mean_x=mean(x,na.rm=TRUE),
         sd_x=sd(x,na.rm=TRUE),
         cv2_x=(sd_x/mean_x)**2) %>%
  mutate(mean_l_x=mean(l_x,na.rm=TRUE),
         sd_l_x=sd(l_x,na.rm=TRUE),
         cv2_l_x=(sd_l_x/mean_l_x)**2) %>%
  mutate(ir=log(q75)-log(q25)) %>%
  ungroup() %>% 
  distinct(condition,promoter,replica,date,datatype,ir,cv2_x,cv2_l_x)

data_toplot_all <- noise_contribution_data %>% 
  filter(!((condition=="glucoseaa020")&(promoter=="med3"))) %>%
  filter(datatype=="allnoisesources") %>% 
  mutate(x=concentration) %>% 
  mutate(l_x=log(concentration)) %>% 
  #mutate(sd_x=sqrt((cov_gg+gfp_nb**2*cov_xx-2*gfp_nb*cov_xg)/length_um**2)) %>% 
  group_by(condition,promoter,replica,date,datatype) %>%
  mutate(q25=quantile(x,0.25)) %>% 
  mutate(q75=quantile(x,0.75)) %>%
  mutate(q50=quantile(x,0.50)) %>%
  mutate(mean_x=mean(x,na.rm=TRUE),
         sd_x=sd(x,na.rm=TRUE),
         cv2_x_all=(sd_x/mean_x)**2) %>% 
  mutate(mean_l_x=mean(l_x,na.rm=TRUE),
         sd_l_x=sd(l_x,na.rm=TRUE),
         cv2_l_x_all=(sd_l_x/mean_l_x)**2) %>%
  mutate(ir_all=log(q75)-log(q25)) %>% 
  ungroup() %>% 
  distinct(condition,promoter,replica,date,ir_all,cv2_x_all,cv2_l_x_all)

data_toplot <- data_toplot %>% 
  left_join(data_toplot_all,by=c("condition","promoter","replica","date")) %>% 
  mutate(ir_norm=ir/ir_all) %>% 
  mutate(cv2_x_norm=cv2_x/cv2_x_all) %>% 
  mutate(cv2_l_x_norm=cv2_l_x/cv2_l_x_all)

data_toplot %>% 
  left_join(condition_label_df,by=c("condition")) %>% 
  left_join(datatype_label_df,by=c("datatype")) %>% 
  ggplot()+
  geom_point(aes(factor(datatype_label,levels=datatypes_labels),cv2_x,col=promoter),alpha=0.5)+
  facet_wrap(~factor(condition_label,levels=conditions_labels),scale="free")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  expand_limits(y=0)+
  xlab("")+
  ylab("CV2 concentration")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/cv2_c_vs_parameters_facet_condition.pdf",width=12,height=6)

data_toplot %>% 
  left_join(condition_label_df,by=c("condition")) %>% 
  ggplot()+
  geom_point(aes(factor(datatype,levels=datatypes),cv2_x,col=promoter),alpha=0.5)+
  geom_line(aes(factor(datatype,levels=datatypes),cv2_x,col=promoter,group=interaction(condition,promoter,replica)),alpha=0.5)+
  facet_wrap(~factor(condition_label,levels=conditions_labels)+promoter,scale="free",ncol=6)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  expand_limits(y=0)+
  xlab("")+
  ylab("CV2 concentration")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/cv2_c_vs_parameters_facet_promoter_condition.pdf",width=12,height=8)

data_toplot %>% 
  left_join(condition_label_df,by=c("condition")) %>% 
  ggplot()+
  geom_point(aes(factor(datatype,levels=datatypes),sqrt(cv2_x),col=promoter),alpha=0.5)+
  geom_line(aes(factor(datatype,levels=datatypes),sqrt(cv2_x),col=promoter,group=interaction(condition,promoter,replica)),alpha=0.5)+
  facet_wrap(~factor(condition_label,levels=conditions_labels)+promoter,scale="free",ncol=6)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  expand_limits(y=0)+
  xlab("")+
  ylab("CV concentration")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/cv_c_vs_parameters_facet_promoter_condition.pdf",width=12,height=8)


data_toplot %>% 
  left_join(condition_label_df,by=c("condition")) %>% 
  ggplot()+
  geom_point(aes(factor(datatype,levels=datatypes),cv2_l_x,col=promoter),alpha=0.5)+
  geom_line(aes(factor(datatype,levels=datatypes),cv2_l_x,col=promoter,group=interaction(condition,promoter,replica)),alpha=0.5)+
  facet_wrap(~factor(condition_label,levels=conditions_labels)+promoter,scale="free",ncol=6)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  expand_limits(y=0)+
  xlab("")+
  ylab("CV2 log-concentration")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/cv2_l_c_vs_parameters_facet_promoter_condition.pdf",width=12,height=8)

data_toplot %>% 
  left_join(condition_label_df,by=c("condition")) %>% 
  ggplot()+
  geom_point(aes(factor(datatype,levels=datatypes),sqrt(cv2_l_x),col=promoter),alpha=0.5)+
  geom_line(aes(factor(datatype,levels=datatypes),sqrt(cv2_l_x),col=promoter,group=interaction(condition,promoter,replica)),alpha=0.5)+
  facet_wrap(~factor(condition_label,levels=conditions_labels)+promoter,scale="free",ncol=6)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  expand_limits(y=0)+
  xlab("")+
  ylab("CV log-concentration")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/cv_l_c_vs_parameters_facet_promoter_condition.pdf",width=12,height=8)

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
  ylab("IR log-concentration")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/IR_l_c_vs_parameters_facet_promoter_condition.pdf",width=12,height=8)


# Normalized version

data_toplot %>% 
  left_join(condition_label_df,by=c("condition")) %>% 
  ggplot()+
  geom_point(aes(factor(datatype,levels=datatypes),cv2_l_x_norm,col=promoter),alpha=0.5)+
  geom_line(aes(factor(datatype,levels=datatypes),cv2_l_x_norm,col=promoter,group=interaction(condition,promoter,replica)),alpha=0.5)+
  facet_wrap(~factor(condition_label,levels=conditions_labels)+promoter,scale="free",ncol=6)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  expand_limits(y=0)+
  xlab("")+
  ylab("Normalized CV2 log-concentration")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/cv2_norm_l_c_vs_parameters_facet_promoter_condition.pdf",width=12,height=8)

data_toplot %>% 
  left_join(condition_label_df,by=c("condition")) %>% 
  ggplot()+
  geom_point(aes(factor(datatype,levels=datatypes),sqrt(cv2_l_x)/sqrt(cv2_l_x_all),col=promoter),alpha=0.5)+
  geom_line(aes(factor(datatype,levels=datatypes),sqrt(cv2_l_x)/sqrt(cv2_l_x_all),col=promoter,group=interaction(condition,promoter,replica)),alpha=0.5)+
  facet_wrap(~factor(condition_label,levels=conditions_labels)+promoter,scale="free",ncol=6)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  expand_limits(y=0)+
  xlab("")+
  ylab("Normalized CV log-concentration")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/cv_norm_l_c_vs_parameters_facet_promoter_condition.pdf",width=12,height=8)

data_toplot %>% 
  left_join(condition_label_df,by=c("condition")) %>% 
  ggplot()+
  geom_point(aes(factor(datatype,levels=datatypes),ir_norm,col=promoter),alpha=0.5)+
  geom_line(aes(factor(datatype,levels=datatypes),ir_norm,col=promoter,group=interaction(condition,promoter,replica)),alpha=0.5)+
  facet_wrap(~factor(condition_label,levels=conditions_labels)+promoter,scale="free",ncol=6)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  expand_limits(y=0)+
  xlab("")+
  ylab("Normalized IR log-concentration")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/IR_norm_l_c_vs_parameters_facet_promoter_condition.pdf",width=12,height=8)


