# We assume data are already loaded in raw_data dataframe.

raw_data %>% 
  left_join(condition_label_df,by=c("condition")) %>% 
  mutate(condition=factor(condition,levels=conditions)) %>% 
  mutate(promoter=factor(promoter,levels=promoters)) %>% 
  mutate(conditions_label=factor(condition_label,levels=conditions_labels)) %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  arrange(condition_label) %>% 
  ggplot(aes(x=alpha_oriented_bbox*60))+
  geom_freqpoly(aes(y=after_stat(density),col=interaction(date)),bins=40)+
  theme_cowplot()+
  scale_y_continuous(trans="log10")+
  facet_wrap(~fct_relevel(condition_label,levels=conditions_labels),scale="free")+
  scale_color_manual(values=c25)+
  ylab("Density")+
  xlab("Growth-rate (min-1)")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))+
  labs(col="Date") 
  
ggsave("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/figures/pdfs/growth_rate_distributions_full_cell_cycle_raw_data_Fcondition.pdf",width=12,height=8,units="cm")

# Plot all growth rate distributions (Stratified by promoter and condition, replica)

```{r}
library(RColorBrewer)
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown",
  "blue","red",
  "yellow","orange","black","red"
)

myframes %>% 
  left_join(label_df,by=c("condition")) %>% 
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
  
ggsave("./raw_data_sanity_checks/raw_data_growth_rate_distributions_replica.pdf",width=24,height=16,units="cm")
```

# Plot all growth rate distributions, facetted per promoter
```{r}
myframes %>% 
  left_join(label_df,by=c("condition")) %>% 
  mutate(condition=factor(condition,levels=conditions)) %>% 
  mutate(promoter=factor(promoter,levels=promoters)) %>% 
  mutate(conditions_label=factor(condition_label,levels=conditions_labels)) %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  arrange(condition_label) %>% 
  ggplot(aes(x=alpha_oriented_bbox*60))+
  geom_freqpoly(aes(y=after_stat(density),group=interaction(date,condition),col=as.character(replica)),binwidth=0.0003)+
  theme_cowplot()+
  #scale_y_continuous(trans="log10")+
  facet_wrap(~promoter,scale="free_y")+
  scale_color_manual(values=c25)+
  ylab("Density")+
  xlab("Growth-rate (min-1)")+
  labs(col="Promoter.date.condition") 

ggsave("./raw_data_sanity_checks/raw_data_growth_rate_facetted_per_promoter.pdf",width=50,height=25,units="cm")
```

# Plot all concentration distributions (raw-data)

```{r}
myframes %>% 
  mutate(concentration=gfp_nb/length_um_oriented_bbox) %>% 
  #filter(concentration>0) %>% 
  left_join(label_df,by=c("condition")) %>% 
  mutate(condition=factor(condition,levels=conditions)) %>% 
  mutate(promoter=factor(promoter,levels=promoters)) %>% 
  mutate(conditions_label=factor(condition_label,levels=conditions_labels)) %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  arrange(condition_label) %>% 
  ggplot(aes(x=concentration))+
  geom_freqpoly(aes(y=after_stat(density),col=as.character(replica)),bins=80)+
  theme_cowplot()+
  scale_y_continuous(trans="log10")+
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


ggsave("./raw_data_sanity_checks/raw_data_concentration_distributions.pdf",width=24,height=16,units="cm")
```

```{r}
myframes %>% 
  group_by(condition,promoter,date) %>% 
  mutate(m_concentration=median(gfp_nb/length_um_oriented_bbox)) %>% 
  ungroup() %>% 
  distinct(condition,promoter,date,m_concentration)
```

```{r}
mean_gr <- myframes %>% 
  #filter(replica==1) %>% 
  group_by(condition,promoter,date) %>% 
  mutate(mean_lambda=compute_weighted_mean(alpha_oriented_bbox,sd_alpha_oriented_bbox)) %>% 
  ungroup() %>% 
  distinct(condition,promoter,date,mean_lambda)

myframes %>% 
  mutate(x=gfp_nb/length_um_oriented_bbox) %>% 
  left_join(mean_gr, by=c("condition","promoter","date")) %>% 
  group_by(condition,date,promoter) %>% 
  mutate(q25=quantile(x,0.25)) %>% 
  mutate(q75=quantile(x,0.75)) %>%
  mutate(q50=quantile(x,0.50)) %>%
  ungroup() %>% 
  distinct(condition,promoter,date,.keep_all=TRUE) %>% 
  mutate(iq=(log(q75)-log(q25))) %>% 
  #mutate(iq_normalized=(q75-q25)/mean_x) %>% 
  ggplot()+
  geom_point(aes(mean_lambda,iq,col=as.character(date)))+
  #geom_line(aes(conditions_label,iq_normalized,col=promoter,group=interaction(promoter,replica)))+
  theme_cowplot()+
  scale_color_manual(values=c25)+
  facet_wrap(~promoter)+
  expand_limits(y=0)+
  ylab("Interquartile range")+
  xlab("Mean growth-rate") +
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))+
  labs(col="Date") 

#ggsave("./denoised_data_distributions/interquartile_range_log_all_dates.pdf",width=12,height=8,units="cm")
```



```{r}
myframes %>% 
  mutate(concentration=gfp_nb/length_um_oriented_bbox) %>% 
  #filter(concentration>0) %>% 
  left_join(label_df,by=c("condition")) %>% 
  mutate(condition=factor(condition,levels=conditions)) %>% 
  mutate(promoter=factor(promoter,levels=promoters)) %>% 
  mutate(conditions_label=factor(condition_label,levels=conditions_labels)) %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  arrange(condition_label) %>% 
  ggplot(aes(x=concentration))+
  geom_freqpoly(aes(y=after_stat(density),col=interaction(date,fct_relevel(condition_label,levels=conditions_labels))),bins=80)+
  theme_cowplot()+
  scale_y_continuous(trans="log10")+
  scale_x_continuous(trans="log10")+
  facet_wrap(~promoter,scale="free",ncol=2)+
  #facet_wrap(~promoter+fct_relevel(condition_label,levels=conditions_labels),scale="free",ncol=4)+
  #facet_grid(fct_relevel(condition_label,levels=conditions_labels)~promoter,scale="free")+
  scale_color_manual(values=c25)+
  ylab("Density")+
  xlab("Concentration (#GFP/um)")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))+
  labs(col='condition') 
  guides(fill = guide_legend(byrow = TRUE))


ggsave("./raw_data_sanity_checks/raw_data_concentration_distributions_logxy.pdf",width=24,height=16,units="cm")
```

# Plot all volumic production distributions (raw-data)

```{r}
myframes %>% 
  filter(mean_q>0) %>% 
  left_join(label_df,by=c("condition")) %>% 
  mutate(condition=factor(condition,levels=conditions)) %>% 
  mutate(promoter=factor(promoter,levels=promoters)) %>% 
  mutate(conditions_label=factor(condition_label,levels=conditions_labels)) %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  arrange(condition_label) %>% 
  ggplot(aes(x=mean_q*60))+
  geom_freqpoly(aes(y=after_stat(density),col=interaction(date,fct_relevel(condition_label,levels=conditions_labels))),bins=80)+
  theme_cowplot()+
  scale_y_continuous(trans="log10")+
  facet_wrap(~promoter,scale="free",ncol=2)+
  #facet_wrap(~promoter+fct_relevel(condition_label,levels=conditions_labels),scale="free",ncol=4)+
  #facet_grid(fct_relevel(condition_label,levels=conditions_labels)~promoter,scale="free")+
  scale_color_manual(values=c25)+
  ylab("Density")+
  xlab("Volumic production  (#GFP/min/um)")+
  theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        legend.spacing.y = unit(0.05, 'cm'))+
  labs(col='condition') 
  guides(fill = guide_legend(byrow = TRUE))


ggsave("./raw_data_sanity_checks/raw_data_production_distributions.pdf",width=24,height=16,units="cm")
```

# Inspect how growth rate, production and concentration change over time.

```{r}
myframes %>% 
  left_join(label_df,by=c("condition")) %>% 
  mutate(condition=factor(condition,levels=conditions)) %>% 
  mutate(promoter=factor(promoter,levels=promoters)) %>% 
  mutate(conditions_label=factor(condition_label,levels=conditions_labels)) %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  arrange(condition_label) %>% 
  ggplot()+
  geom_point(aes(frame,alpha_oriented_bbox*60))+
  theme_cowplot()+
  #scale_y_continuous(trans="log10")+
  facet_wrap(~date+fct_relevel(condition_label,levels=conditions_labels)+promoter,scale="free")+
  scale_color_manual(values=c25)+
  ylab("Growth-rate (min-1)")+
  xlab("Frame")
  #theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        #legend.spacing.y = unit(0.05, 'cm'))+
  #guides(fill = guide_legend(byrow = TRUE))

ggsave("./raw_data_sanity_checks/growth_rate_vs_time.pdf",width=50,height=50,limitsize=FALSE)
```
# Inspect how growth rate, problem with 6300

```{r}
myframes %>% 
  filter(condition=="glucoseaa020",date=="20230302",promoter=="6300") %>% 
  left_join(label_df,by=c("condition")) %>% 
  mutate(condition=factor(condition,levels=conditions)) %>% 
  mutate(promoter=factor(promoter,levels=promoters)) %>% 
  mutate(conditions_label=factor(condition_label,levels=conditions_labels)) %>% 
  distinct(cell,.keep_all=TRUE) %>% 
  arrange(condition_label) %>% 
  ggplot()+
  geom_point(aes(frame,alpha_oriented_bbox*60))+
  theme_cowplot()+
  #scale_y_continuous(trans="log10")+
  facet_wrap(~date+fct_relevel(condition_label,levels=conditions_labels)+promoter+gl_id,nrow=1)+
  scale_color_manual(values=c25)+
  ylab("Growth-rate (min-1)")+
  xlab("Frame")
  #theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        #legend.spacing.y = unit(0.05, 'cm'))+
  #guides(fill = guide_legend(byrow = TRUE))

ggsave("./figures/growth_rate_vs_time.pdf",width=50,height=10,limitsize=FALSE)
```
# Long tail for growth rate in glucoseaa, 20230302, why?

```{r}
myframes %>% 
  filter(condition=="glucoseaa020",date=="20230302") %>% 
  ggplot(aes(x=alpha_oriented_bbox*60,col=promoter))+
  geom_freqpoly(aes(y=after_stat(density)),bins=80)+
  theme_cowplot()+
  scale_y_continuous(trans="log10")
```
```{r}
myframes %>% 
  filter(condition=="glucoseaa020",date=="20230302") %>% 
  filter(alpha_oriented_bbox*60<0.01) %>% 
  ggplot(aes(x=alpha_oriented_bbox*60,col=promoter))+
  geom_freqpoly(aes(y=after_stat(density)),bins=80)+
  theme_cowplot()+
  scale_y_continuous(trans="log10")
```
```{r}
selected_cells <- myframes %>% 
  filter(condition=="glucoseaa020",date=="20230302") %>% 
  filter(alpha_oriented_bbox*60<0.01) %>%
  distinct(cell,promoter) 
  #filter(row_number()<10)

myframes %>% 
  semi_join(selected_cells,by=c("cell","promoter")) %>% 
  ggplot()+
  geom_point(aes(time_sec,length_um))+
  facet_wrap(~cell+promoter,scale="free")
```

```{r}
selected_cells <- myframes %>% 
  filter(condition=="glucoseaa020",date=="20230302") %>% 
  filter(alpha_oriented_bbox*60<0.01) %>%
  group_by(cell) %>% 
  filter(row_number()==1) %>% 
  ungroup() %>% 
  distinct(cell,promoter,mean_cell_rank,frame) %>% View()

myframes %>% 
  semi_join(selected_cells,by=c("cell","promoter"))

```

# Inspect traces

```{r}
myframes %>% 
  filter(condition=="glucoseaa020",date=="20230302",promoter=="6300") %>% 
  filter(alpha_oriented_bbox*60<0.01) %>% 
  left_join(label_df,by=c("condition")) %>% 
  mutate(condition=factor(condition,levels=conditions)) %>% 
  mutate(promoter=factor(promoter,levels=promoters)) %>% 
  mutate(conditions_label=factor(condition_label,levels=conditions_labels)) %>% 
  #distinct(cell,.keep_all=TRUE) %>% 
  arrange(condition_label) %>% 
  ggplot()+
  geom_point(aes(frame,length_um))+
  geom_line(aes(frame,length_um,group=cell))+
  theme_cowplot()+
  #scale_y_continuous(trans="log10")+
  facet_wrap(~date+fct_relevel(condition_label,levels=conditions_labels)+promoter+gl_id,ncol=1)+
  scale_color_manual(values=c25)+
  ylab("Growth-rate (min-1)")+
  xlab("Frame")
  #theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        #legend.spacing.y = unit(0.05, 'cm'))+
  #guides(fill = guide_legend(byrow = TRUE))

ggsave("./figures/growth_rate_vs_time.pdf",width=20,height=40,limitsize=FALSE)
```

# Checking offset is well set

```{r}
#Problem with #6300
myframes2 <- myframes %>% 
  filter(promoter=="rrnB",condition=="acetate005",date=="20230228") %>% 
  left_join(label_df,by=c("condition")) %>% 
  mutate(condition=factor(condition,levels=conditions)) %>% 
  mutate(promoter=factor(promoter,levels=promoters)) %>% 
  mutate(conditions_label=factor(condition_label,levels=conditions_labels)) %>% 
  arrange(condition_label) %>% 
  filter(alpha_oriented_bbox*60<0.01)

myframes %>% 
  filter(promoter=="rrnB",condition=="acetate005",date=="20230228") %>% 
  #filter(promoter=="6300",condition=="glucose020") %>% 
  left_join(label_df,by=c("condition")) %>% 
  mutate(condition=factor(condition,levels=conditions)) %>% 
  mutate(promoter=factor(promoter,levels=promoters)) %>% 
  mutate(conditions_label=factor(condition_label,levels=conditions_labels)) %>% 
  arrange(condition_label) %>% 
  #mutate(mean_cell_length=mean(length_um,na.rm=TRUE),
  #       sd_cell_length=sd(length_um,na.rm=TRUE),
  #       min_threshold=mean_cell_length-4*sd_cell_length,
  #       max_threshold=mean_cell_length+4*sd_cell_length) %>% 
  #mutate(filamenting=ifelse(length_um>max_threshold,TRUE,FALSE)) %>% 
  #group_by(cell) %>% 
  #filter(!any(filamenting)) %>% 
  #ungroup() %>% 
  ggplot()+
  geom_line(aes(frame,length_um,group=cell),alpha=0.5)+
  #geom_line(data=myframes2,aes(frame,length_um,group=cell),col="red")+
  theme_cowplot()+
  scale_y_continuous(trans="log10")+
  #facet_wrap(~gl_id,ncol=1)+
  #scale_color_manual(values=c25)+
  ylab("Length (um)")+
  xlab("Frame")
  #theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        #legend.spacing.y = unit(0.05, 'cm'))+
  #guides(fill = guide_legend(byrow = TRUE))

ggsave("./figures/6300_length_vs_time.pdf",width=50,height=10,limitsize=FALSE)
```

```{r}
myframes2 <- myframes %>% 
  filter(promoter=="rrnB",condition=="acetate005",date=="20230228") %>% 
  left_join(label_df,by=c("condition")) %>% 
  mutate(condition=factor(condition,levels=conditions)) %>% 
  mutate(promoter=factor(promoter,levels=promoters)) %>% 
  mutate(conditions_label=factor(condition_label,levels=conditions_labels)) %>% 
  arrange(condition_label) %>% 
  filter(alpha_oriented_bbox*60<0.01)


myframes %>% 
  filter(promoter=="rrnB",condition=="acetate005",date=="20230228") %>% 
  left_join(label_df,by=c("condition")) %>% 
  mutate(condition=factor(condition,levels=conditions)) %>% 
  mutate(promoter=factor(promoter,levels=promoters)) %>% 
  mutate(conditions_label=factor(condition_label,levels=conditions_labels)) %>% 
  arrange(condition_label) %>% 
  ggplot()+
  geom_line(aes(length_um,vertical_top,group=cell),alpha=0.5)+
  #geom_line(data=myframes2,aes(length_um,vertical_top,group=cell),col="red")+
  theme_cowplot()+
  #scale_y_continuous(trans="log10")+
  #facet_wrap(~gl_id,ncol=1)+
  #scale_color_manual(values=c25)+
  ylab("Vertical top")+
  xlab("Length (um)")
  #theme(text = element_text(size = 8),axis.text = element_text(size = 5),
        #legend.spacing.y = unit(0.05, 'cm'))+
  #guides(fill = guide_legend(byrow = TRUE))

ggsave("./figures/6300_length_vs_time.pdf",width=50,height=10,limitsize=FALSE)
```


```{r}
myframes_complete %>% 
  filter(promoter=="rrnB",condition=="acetate005",date=="20230228") %>% 
  filter(frame<250) %>% 
  group_by(frame) %>% 
  count()
  
  

  
  
```
# Compute inference parameters best guess

```{r}
myframes %>% 
  filter(condition %in% conditions) %>% 
  group_by(date,condition,promoter,cell) %>% 
  sample_n(1) %>% 
  ungroup %>% 
  group_by(date,condition,promoter) %>% 
  mutate(wtot=sum(1/sd_alpha_oriented_bbox**2,na.rm = TRUE)) %>%
  mutate(w=1/sd_alpha_oriented_bbox**2) %>% 
  mutate(mean_l=sum(alpha_oriented_bbox*w)*60/wtot) %>% 
  mutate(#mean_l=compute_weighted_mean(alpha_oriented_bbox,sd_alpha_oriented_bbox)*60,
         #mean_l=mean(alpha_oriented_bbox)*60,
         mean_doubling_time=log(2)/mean_l) %>% 
  ungroup() %>% 
  distinct(date,condition,promoter,mean_doubling_time,mean_l)
```

# In raw data, exponential fit errors bars are too small for slow growing cells, which have an important weight in the weighted mean, which distorts it.

```{r}
myframes %>% 
  group_by(date,condition,promoter) %>% 
  mutate(wtot=sum(1/sd_alpha_oriented_bbox**2,na.rm = TRUE)) %>%
  ungroup() %>% 
  ggplot()+
  geom_point(aes(alpha_oriented_bbox*60,(1/sd_alpha_oriented_bbox**2)/wtot))+
  facet_wrap(~promoter+condition,scale="free")

ggsave("./figures/alpha_sd_alpha.pdf",height=30,width=30,limitsize = FALSE)
```

What's happening is that, some cells in acetate, are characterized by very slow growth rate. And when they do, the uncertainty on their growth rate is very small. Which turns into a huge weight in the computation of the weighted mean. But... these cells actually do not really grow. Let's take a look at them in the data set.

```{r}
conditions <- c("acetate005")
promoters <- c("hi1")

myframes <- myframes %>% 
  mutate(date=as.character(date)) %>% 
  mutate(dataset=paste(promoter,condition,date,sep="_"))

dataset <- myframes %>% distinct(dataset) %>% .$dataset
dataset <- c("hi1_acetate005_20221121")

library(RColorBrewer)
c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

for(.d in dataset){
  
ncells <- length(data <- myframes %>% 
  filter(dataset==.d) %>% 
  distinct(cell) %>% 
    .$cell)

scale <- rep(c25,length.out=ncells)
  
plot_length <- myframes %>% 
  mutate(slow=ifelse(alpha_oriented_bbox*60<0.001,1,0.1)) %>% 
  filter(dataset==.d) %>% 
  ggplot()+
  geom_line(aes(time_sec/3600,length_um,group=cell,col=cell,alpha=slow))+
  #geom_line(aes(time_sec/3600,length_um,group=cell,col=cell),alpha=1)+
  theme_cowplot()+
  scale_y_continuous(trans="log10")+
  ylab("Length (um)")+
  xlab("Time (h)")+
  scale_color_manual(values=scale)+
  labs(subtitle=sprintf("%s",.d))+
  theme(legend.position = "none")

plot <- plot_grid(plot_lengtncol = 1)

dest <- as.character(sprintf("./figures/%s.pdf",.d))
ggsave(dest,plot_length,width=25,height=10,limitsize = FALSE)}
```
They look fine. They only grow slow. They look pretty straight too. They actually do speed up at the end.
Another way to look at this, is to question the r2.

```{r}
myframes %>% 
  distinct(cell,.keep_all = TRUE) %>% 
  mutate(mean_cell_rank=factor(as.character(mean_cell_rank),levels=as.character(c(0:10)))) %>% 
  ggplot()+
  geom_point(aes(alpha_oriented_bbox*60,logl_time_r2,col=mean_cell_rank),alpha=0.5)+
  facet_wrap(~promoter+condition,scale="free")

ggsave("./figures/alpha_logl_r2.pdf",height=30,width=30,limitsize = FALSE)
```
Cells that tend to grow slow, also behave less "exponential". They are also the mother cells.
Then I can run the analysis, not looking at them, but only at the cells above. A better way to compute the mean doubling time, is then to compute the mean of the division time.

# Computing mean division time

```{r}
conditions <- c("acetate005","glycerol040","glucose020","glucoseaa020")
conditions_labels <- c("acetate","glycerol","glucose","glucose+a.a")
promoters <- c("hi1","hi3","med2","med3","rrnB","rplN")
label_df <- tibble(condition=conditions,condition_label=conditions_labels)

myframes %>% 
  filter(condition %in% conditions) %>% 
  group_by(date,condition,promoter,cell) %>% 
  sample_n(1) %>% 
  ungroup %>% 
  group_by(condition,promoter) %>% 
  mutate(mean_division_time=mean(div_time)) %>% 
  ungroup() %>% 
  distinct(date,condition,promoter,mean_division_time)
```


