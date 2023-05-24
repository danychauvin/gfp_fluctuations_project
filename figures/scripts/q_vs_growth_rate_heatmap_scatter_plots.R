#We want to plot, a heat-map of actual growth-rate versus production (mean).For each condition for hi1.
#We assume here that data were already imported. Using load_simulation_data.R

datatypes <- c("allnoisesources","nodivisionnoise","nogrowthnoise","noprodnoise","nonoise")
replicas <- c("1","2")
#Correct for mistake with rrnB third replica
experimental_data <- experimental_data %>% 
  mutate(replica=ifelse(replica==3,1,replica))

toplot <- experimental_data %>% 
  filter(!((condition=="glucoseaa020")&(promoter=="med3"))) %>%
  filter(promoter!="6300")
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

toplot %>% 
  boxPlotXY(.,lambda,q,10)+ggpubr::stat_cor(method = "pearson",size=3)+
  facet_wrap(~promoter+condition+replica,scales="free")+
  theme(aspect.ratio = 1)

ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/fig2/q_vs_lambda_box.pdf",width=24,height=16)