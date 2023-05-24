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

mean_lambda_df <- myparameters %>% 
  filter(name=="mean_lambda") %>% 
  select(promoter,condition,date,final) %>% 
  rename(mean_lambda=final) %>% 
  left_join(doubling_times_df,by=c("promoter","condition","date")) %>% 
  select(condition,promoter,date,mean_lambda)

mean_lambda_err_df <- myerrors %>% 
  filter(epsilon==0.05) %>% 
  mutate(mean_lambda_err=sqrt(mean_lambda)) %>% 
  select(condition,promoter,date,mean_lambda_err)

toplot <- left_join(gamma_lambda_df,gamma_lambda_err_df,by=c("condition","date","promoter")) %>% 
  left_join(mean_lambda_df,by=c("condition","date","promoter")) %>% 
  left_join(mean_lambda_err_df,by=c("condition","date","promoter"))

toplot %>% 
  filter(1/gamma_lambda<1500) %>% 
  mutate(autocor_time=1/gamma_lambda,
         sd_autocor_time=abs(1/(gamma_lambda**2)*gamma_lambda_err)) %>% 
  ggplot()+
  geom_point(aes(mean_lambda,autocor_time,col=date))+
  geom_errorbar(aes(x=mean_lambda,ymin=autocor_time-sd_autocor_time,ymax=autocor_time+sd_autocor_time,col=date))+
  geom_errorbar(aes(y=autocor_time,xmin=mean_lambda-mean_lambda_err,xmax=mean_lambda+mean_lambda_err,col=date))+
  theme(text = element_text(size = 2),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))+
  theme_cowplot()+
  scale_colour_manual(values=c25)+
  xlab("Mean growth-rate (min-1)")+
  ylab("Autocorrelation time (min)")

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/fig2/correlation_time_vs_mean_growth_rate_Cdate.pdf",width=8,height=5)

toplot %>% 
  filter(1/gamma_lambda<1500) %>% 
  mutate(autocor_time=1/gamma_lambda,
         sd_autocor_time=abs(1/(gamma_lambda**2)*gamma_lambda_err)) %>% 
  ggplot()+
  geom_point(aes(mean_lambda,autocor_time,col=promoter))+
  geom_errorbar(aes(x=mean_lambda,ymin=autocor_time-sd_autocor_time,ymax=autocor_time+sd_autocor_time,col=promoter))+
  geom_errorbar(aes(y=autocor_time,xmin=mean_lambda-mean_lambda_err,xmax=mean_lambda+mean_lambda_err,col=promoter))+
  theme(text = element_text(size = 2),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))+
  theme_cowplot()+
  scale_colour_manual(values=c25)+
  xlab("Mean growth-rate (min-1)")+
  ylab("Autocorrelation time (min)")

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/fig2/correlation_time_vs_mean_growth_rate_Cpromoter.pdf",width=8,height=5)
