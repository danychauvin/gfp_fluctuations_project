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
  filter(q-sqrt(abs(cov_qq))>0) %>% 
  mutate(l_q=log(q))

plot_heatscatter <- function(.dataset){
  
  x <- toplot %>%
    filter(dataset==.dataset) %>% 
    .$lambda
  
  y <- toplot %>%
    filter(dataset==.dataset) %>%
    .$l_q
  
  nb_cell_cycles <- nrow(toplot %>%
                           filter(dataset==.dataset) %>%
                           distinct(cell))
  
  nb_measurement <- nrow(toplot %>%
                           filter(dataset==.dataset))
  
  nb_gl <- nrow(toplot %>%
                  filter(dataset==.dataset) %>%
                  distinct(gl_id))
  
  df <- data.frame(x,y)
  
  ## Use densCols() output to get density at each point
  x <- densCols(x,y, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  
  ## Map densities to colors
  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                              "#FCFF00", "#FF9400", "#FF3100"))(256)
  df$col <- cols[df$dens]
  
  ## Plot it, reordering rows so that densest points are plotted on top
  #pdf(file=paste("./figures/",cond,".pdf",sep=""))
  par(tck=-0.05)#,mar=c(2,2,2,2) + 0.1,mai=c(0.1,0.1,0.1,0.1))
  plot(y~x, data=df[order(df$dens),], pch=20, col=col, cex.lab=1, xlab="Exponential growth-rate (/min)",ylab="log GFP Volumic production (GFP/min/um)",yaxt="n",xaxt="n",main=sprintf("%s\n %s obs., %s cycles, %s gls",.dataset,nb_measurement,nb_cell_cycles,nb_gl),cex.main=1)
  #axis(2, mgp=c(3, .5, .5))
  #axis(1, mgp=c(3, .5, .5))
  axis(2, padj = 0)
  axis(1, padj = 0)
  #axis(tck=-0.01)
  #cex.main=0.50,main=sprintf("%s, %s observations, %s cycles, %s gls",.dataset,nb_measurement,nb_cell_cycles,nb_gl)
  #dev.off()}
}

#pdf( "./all_heatmap_scatter.pdf", width = 4, height = 4 )
#plot_heatscatter("hi1_acetate005_20221121")
#dev.off()

plot_empty <- function(){
  plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10))}

#all in one picture
#1 cm = 0.39 inch
pdf(file="/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/fig2/log_q_vs_growth_rate_heatmap_scatter_plots_stringent_filter_above_0.pdf",width=4*8,height=4*6)
par(mfrow=c(6,8))

for(.p in unique(toplot$promoter)){
  for(.c in conditions){
    for(.r in replicas){
      .d <- unique(toplot %>% filter(promoter==.p,condition==.c,replica==.r) %>% .$dataset)
      print(.d)
      if(length(.d)==0){
        plot_empty()
      }else{
        plot_heatscatter(.d)}
    }}}
  
dev.off()

