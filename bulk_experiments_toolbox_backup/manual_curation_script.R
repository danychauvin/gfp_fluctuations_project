#Start manual curation by running the following loop
#Beginning of the loop
for(w in unique(raw_data$well)){
  
  r <- raw_data %>% 
    filter(well==w)
  e <- mydata_exp_growth %>% 
    filter(well==w)
  fo <- raw_data_fluo_od %>% 
    filter(well==w)
  #---------------------------------------------------------------A
  p1a <- r %>% 
    ggplot() +
    geom_point(aes(time_min/60,log(corrected_od)),alpha=0.2,col="black")+
    theme_cowplot()
  
  #--grid
  p1 <- p1a+
    geom_point(data=e,aes(time_min/60,log(corrected_od)),col="green")+
    stat_smooth(data=e,aes(time_min/60,log(corrected_od)),method="lm_left", fullrange=TRUE,col='blue',alpha=0.2)+
    ggpubr::stat_cor(data=e,aes(time_min/60,log(corrected_od)),method = "pearson",size=4, label.y.npc = 0.6,label.x.npc = 0)+
    ggpubr::stat_regline_equation(data=e,aes(time_min/60,log(corrected_od)),size=4,label.y.npc = 0.8,label.x.npc = 0)+
    xlab("time (h)")+
    ylab("log(c_od)")
  #---------------------------------------------------------------B
  p2a <- r %>%
    filter(well==w) %>%
    ggplot() +
    geom_point(aes(time_min/60,log(corrected_fluo)),alpha=0.2,col="black")+
    theme_cowplot()
  
  p2 <- p2a+
    geom_point(data=e,aes(time_min/60,log(corrected_fluo)),col="green")+
    stat_smooth(data=e,aes(time_min/60,log(corrected_fluo)),method="lm_left", fullrange=TRUE,col='blue',alpha=0.2)+
    ggpubr::stat_cor(data=e,aes(time_min/60,log(corrected_fluo)),method = "pearson",size=4, label.y.npc = 0.6,label.x.npc = 0)+
    ggpubr::stat_regline_equation(data=e,aes(time_min/60,log(corrected_fluo)),size=4,label.y.npc = 0.8,label.x.npc = 0)+
    xlab("time (h)")+
    ylab("log(c_fluo)")
  
  #---------------------------------------------------------------C
  p3a <- r %>%
    filter(well==w) %>%
    ggplot() +
    geom_point(aes(corrected_od,corrected_fluo),alpha=0.2,col="black")+
    theme_cowplot()
  p3b <- e %>% ####
    filter(well==w)
  
  p3 <- p3a+
    geom_point(data=e,aes(corrected_od,corrected_fluo),col="green")+
    stat_smooth(data=e,aes(corrected_od,corrected_fluo),method="lm_left", fullrange=TRUE,col='blue',alpha=0.2)+
    ggpubr::stat_cor(data=e,aes(corrected_od,corrected_fluo),method = "pearson",size=4, label.y.npc = 0.6,label.x.npc = 0)+
    ggpubr::stat_regline_equation(data=e,aes(corrected_od,corrected_fluo),size=4,label.y.npc = 0.8,label.x.npc = 0)+
    xlab("log(c_od)")+
    ylab("fluo")
  #---------------------------------------------------------------D
  p4a <- r %>%
    filter(well==w) %>%
    ggplot() +
    geom_point(aes(corrected_od,corrected_fluo),alpha=0.2,col="black")+
    theme_cowplot()
  p4b <- fo %>% ####
    filter(well==w)
  
  p4 <- p4a +
    geom_point(data=fo,aes(corrected_od,corrected_fluo),col="red")+
    stat_smooth(data=fo,aes(corrected_od,corrected_fluo),method="lm_left", fullrange=TRUE,col='blue',alpha=0.2)+
    ggpubr::stat_cor(data=fo,aes(corrected_od,corrected_fluo),method = "pearson",size=4, label.y.npc = 0.6,label.x.npc = 0)+
    ggpubr::stat_regline_equation(data=fo,aes(corrected_od,corrected_fluo),size=4,label.y.npc = 0.8,label.x.npc = 0)+
    geom_point(aes(x=0, y=0), colour="red",shape="+",size=4)+
    xlab("log(c_od)")+
    ylab("fluo")
  
  title <- ggplot()+
    draw_label(paste(w,unique(r$condition),sep="_"),size=15)
  titletxt <- paste(w,unique(r$condition))
  
  pg <- plot_grid(p1,p2,p3,p4, labels =c("A","B","C","D"))
  pgg <- plot_grid(title,pg,nrow=2,rel_heights = c(1,10))
  print(pgg)
  
  good <- readline(prompt=sprintf("Is %s Valid?",titletxt))
  ggsave(sprintf("./curated_figures/%s_%s.png",good,titletxt))
  
  raw_data <- raw_data %>% 
    mutate(valid=ifelse(well==w,good,valid))
}
# End of the loop


