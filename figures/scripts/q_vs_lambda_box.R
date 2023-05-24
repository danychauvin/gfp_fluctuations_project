#We want to plot, a heat-map of actual growth-rate versus production (mean).For each condition for hi1.
#We assume here that data were already imported. Using load_simulation_data.R

datatypes <- c("allnoisesources","nodivisionnoise","nogrowthnoise","noprodnoise","nonoise")
replicas <- c("1","2")
#Correct for mistake with rrnB third replica
experimental_data <- experimental_data %>% 
  mutate(replica=ifelse(replica==3,1,replica))

toplot <- experimental_data %>% 
  filter(!((condition=="glucoseaa020")&(promoter=="med3"))) %>%
  filter(promoter!="6300") %>% 
  mutate(replica=as.character(replica))
  #filter(q-sqrt(abs(cov_qq))>0) %>% 
  #mutate(l_q=log(q))

boxPlotXY <- function(df,xvar,yvar,binWidth){
  xvar <- enquo(xvar)
  yvar <- enquo(yvar)
  
  p <- df %>% 
    ggplot(aes(x=!!xvar, y=!!yvar)) +
    geom_boxplot(fill="skyblue", aes(group = cut_number(!!xvar, n=binWidth)))
  return(p)
}

plot_box_condition <- function(.cond){

  toplot %>% 
    filter(condition==.cond) %>% 
    boxPlotXY(.,lambda,q,10)+ggpubr::stat_cor(method = "pearson",size=3)+
    facet_wrap(~factor(promoter,levels=promoters)+
                factor(condition,levels=conditions)+
                factor(replica,levels=c("1","2")),scales="free")+
    theme(aspect.ratio = 1)+
    theme(text = element_text(size = 2),axis.text = element_text(size = 5),
          legend.spacing.y = unit(0.05, 'cm'))+
    theme_cowplot()
  
  ggsave(sprintf("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/fig2/q_vs_lambda_box_%s.pdf",.cond),width=24,height=16)
}

plot_box_condition("acetate005")
plot_box_condition("glycerol040")
plot_box_condition("glucose020")
plot_box_condition("glucoseaa020")