# Manual curation of 96 well plates data
# David Schar, Dany Chauvin
# 20210302

# Loading packages and functions
```{r message=FALSE, warning=FALSE}
source("./load_packages_functions.R")
source("./modeling_functions.R")
```

# Importing data
Here data of your choice can be imported, by providing correct paths. Here importing data which are stored in `path_to_data_list`.

```{r message=FALSE, warning=FALSE}
path_to_data_list <- "/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/constitexpr/20210208_constitexpr_testrun2/dataList.csv"
source("./import_data.R")
```
# Preparing the data for manual curation
We add `valid` which will be use to store our choice regarding whether or not the data is valid for inference.

```{r}
mydata <- mydata %>% 
  mutate(valid=NA) %>% 
  mutate(well=paste(plate,row,column,strain,sep="_"))
```

# Correct for blank

Ideally, the blank should be computed well per well, but this can be difficult in practice. Here using the value we obtained by computing the mean OD value of wells in which no cells were inoculating. This mean value should be obtained independently for each experiment. We showed it does not seem to be dependent on the condition used.

The same should be done for the fluorescence of empty wells. Here copy pasting the values that I got from testrun2.

```{r}
mean_blank_od <- 0.03840125
mydata <- mydata %>% 
  mutate(corrected_od=od-mean_blank_od)

mean_blank_fluo <- 35.3075
mydata <- mydata %>% 
  mutate(corrected_fluo=fluo-mean_blank_fluo)
```

# Plot log(corrected_od) & log(fluo) for all plates

In order to quickly check the data, we first plot all kinetics.

```{r}
plates <- unique(mydata$plate)
```

```{r message=FALSE,warning=FALSE}
for(i in plates){
  mydata_to_plot <- mydata %>% 
    filter(plate==i) %>% 
    filter(time_min/60<30)
  
  condition <- unique(mydata_to_plot$condition)
  
  print(
  mydata_to_plot %>%
  ggplot()+
  geom_point(aes(time_min/60,corrected_od))+
  facet_grid(row~column)+ 
  theme(strip.text.x = element_blank())+
  labs(subtitle=paste(i,condition,sep=" "))+
    xlab("time (h)")+
    ylab("corrected OD"))
  }
```
I do the same for the fluorescence (corrected fluorescence can be used to figure out the growth rate).

```{r message=FALSE,warning=FALSE}
for(i in plates){
  mydata_to_plot <- mydata %>% 
    filter(plate==i)%>% 
    filter(time_min/60<30)
  
  condition <- unique(mydata_to_plot$condition)
  
  print(
  mydata_to_plot %>%
  ggplot()+
  geom_point(aes(time_min,log(corrected_fluo)))+
  facet_grid(row~column)+ 
  theme(strip.text.x = element_blank())+
  labs(subtitle=paste(i,condition,sep=" "))+
    xlab("time (h)")+
    ylab("corrected OD"))
  }
```

# ACETATE

Based on the previous results, one can set parameters to focus on relevant part of the kinetics experimentally acquired.

Manual curation is performed below per condition.

```{r}
cond <- "acetate"
```

### Keeping all data

```{r message=FALSE,warning=FALSE}
raw_data <- mydata %>%
  filter(condition==cond) %>% #keeping only one relevant condition
  filter(!grepl("empty",strain)) %>% #filter out empty wells
  filter(max(od)>0.05) %>% #to get rid of wells in which there are no growth at all 
  filter(time_min/60<30) #clearly nothing intersting to see after 30 hours.
```

### Keeping only exponential growth regime

To identify exponential growth, we arbitrarily define a time threshold (to discard experimental noise probably due to condensation at the beginning of the kinetics), `time_min_threshold` and coefficient 'rel' so that only observations that are characterized by corrected_od<rel*maximum_corrected_od are considered.

```{r message=FALSE,warning=FALSE}
#parameters for exponential growth
rel <- 0.6
time_min_threshold <- 5*60 #min

mydata_exp_growth <- mydata %>% 
  filter(condition==cond) %>% # select only one condition
  filter(!grepl("empty",strain)) %>% #filter out empty wells
  filter(max(od)>0.05) %>% #to get rid of wells in which there are no growth at all 
  filter(time_min>time_min_threshold) %>% #discard beginning of the kinetics
  group_by(well) %>%
  arrange(time_min) %>% # the following is to discard all observations after corrected_od>rel*max_od
  mutate(global_max=max(corrected_od,na.rm=TRUE)) %>%
  mutate(relative_max=rel*global_max) %>% 
  mutate(above_max_time=ifelse(corrected_od>=relative_max,time_min,NA)) %>% 
  mutate(time_limit=min(above_max_time,na.rm=TRUE)) %>% 
  filter(time_min<time_limit) %>% 
  ungroup()
```

Checking results.

```{r message=FALSE,warning=FALSE}
plates <- unique(mydata_exp_growth$plate)

for(i in plates){
  mydata_to_plot <- mydata_exp_growth %>% 
    filter(plate==i)
  
  condition <- unique(mydata_to_plot$condition)
  
  print(
  mydata_to_plot %>%
  ggplot()+
  geom_point(aes(time_min/60,log(corrected_od)))+
  facet_grid(row~column)+ 
  theme(strip.text.x = element_blank())+
  labs(subtitle=paste(i,condition,sep=" ")))
  
  print(
  mydata_to_plot %>%
  ggplot()+
  geom_point(aes(time_min/60,log(corrected_fluo)))+
  facet_grid(row~column)+ 
  theme(strip.text.x = element_blank())+
  labs(subtitle=paste(i,condition,sep=" ")))
  }
```

### Keeping data meant to be used for fluo vs corrected_od inference

We noticed that the linear relationship between fluo and corrected_od seems to extent to early stationary phase. By iterating and looking carefully at the fits, we came up with arbitrary `min_corrected_od` and `max_corrected_od` values where a linear fit is acceptable between fluo and corrected_od. The value of these parameters depends on the condition.

```{r message=FALSE,warning=FALSE}
min_corrected_od <- 0.05
max_corrected_od <- 0.2

raw_data_fluo_od <- mydata %>% 
  filter(between(corrected_od,min_corrected_od,max_corrected_od)) %>% 
  filter(condition==cond) %>% 
  filter(!grepl("empty",strain)) %>% #filter out empty wells
  filter(max(od)>min_corrected_od) %>% #to get rid of wells in which there are no growth at all 
  filter(time_min/60<30)
```

### Performing the manual curation

Run the first loop in manual_curation_script.R to fill the `valid` criteria. Data within a well are considered valid is the linear fit on the red data points (figure D), is acceptable.

Once the manual curation has been done, raw_data can be saved into a csv file. For the inference only the observations for which valid==TRUE will be considered.

```{r}
# Saving data
raw_data %>% 
  readr::write_csv("./ace_curated.csv")
```

Now, the same can be done for other conditions.

# GLYCEROL

Based on the previous results, one can set parameters to focus on relevant part of the kinetics experimentally acquired.

Manual curation is performed below per condition.

```{r}
cond <- "glycerol"
```

### Keeping all data

```{r message=FALSE,warning=FALSE}
raw_data <- mydata %>%
  filter(condition==cond) %>% #keeping only one relevant condition
  filter(!grepl("empty",strain)) %>% #filter out empty wells
  filter(max(od)>0.05) %>% #to get rid of wells in which there are no growth at all 
  filter(time_min/60<30) #clearly nothing intersting to see after 30 hours.
```

### Keeping only exponential growth regime

To identify exponential growth, we arbitrarily define a time threshold (to discard experimental noise probably due to condensation at the beginning of the kinetics), `time_min_threshold` and coefficient 'rel' so that only observations that are characterized by corrected_od<rel*maximum_corrected_od are considered.

```{r message=FALSE,warning=FALSE}
#parameters for exponential growth
rel <- 0.6
time_min_threshold <- 2*60

mydata_exp_growth <- mydata %>% 
  filter(condition==cond) %>% # select only one condition
  filter(!grepl("empty",strain)) %>% #filter out empty wells
  filter(max(od)>0.05) %>% #to get rid of wells in which there are no growth at all 
  filter(time_min>time_min_threshold) %>% #discard beginning of the kinetics
  group_by(well) %>%
  arrange(time_min) %>% # the following is to discard all observations after corrected_od>rel*max_od
  mutate(global_max=max(corrected_od,na.rm=TRUE)) %>%
  mutate(relative_max=rel*global_max) %>% 
  mutate(above_max_time=ifelse(corrected_od>=relative_max,time_min,NA)) %>% 
  mutate(time_limit=min(above_max_time,na.rm=TRUE)) %>% 
  filter(time_min<time_limit) %>% 
  ungroup()
```

Checking results.

```{r message=FALSE,warning=FALSE}
plates <- unique(mydata_exp_growth$plate)

for(i in plates){
  mydata_to_plot <- mydata_exp_growth %>% 
    filter(plate==i)
  
  condition <- unique(mydata_to_plot$condition)
  
  print(
  mydata_to_plot %>%
  ggplot()+
  geom_point(aes(time_min/60,log(corrected_od)))+
  facet_grid(row~column)+ 
  theme(strip.text.x = element_blank())+
  labs(subtitle=paste(i,condition,sep=" ")))
  
  print(
  mydata_to_plot %>%
  ggplot()+
  geom_point(aes(time_min/60,log(corrected_fluo)))+
  facet_grid(row~column)+ 
  theme(strip.text.x = element_blank())+
  labs(subtitle=paste(i,condition,sep=" ")))
  }
```

### Keeping data meant to be used for fluo vs corrected_od inference

We noticed that the linear relationship between fluo and corrected_od seems to extent to early stationary phase. By iterating and looking carefully at the fits, we came up with arbitrary `min_corrected_od` and `max_corrected_od` values where a linear fit is acceptable between fluo and corrected_od. The value of these parameters depends on the condition.

```{r message=FALSE,warning=FALSE}
min_corrected_od <- 0.05
max_corrected_od <- 0.4

raw_data_fluo_od <- mydata %>% 
  filter(between(corrected_od,min_corrected_od,max_corrected_od)) %>% 
  filter(condition==cond) %>% 
  filter(!grepl("empty",strain)) %>% #filter out empty wells
  filter(max(od)>min_corrected_od) %>% #to get rid of wells in which there are no growth at all 
  filter(time_min/60<30)
```

### Performing the manual curation

Run the first loop in manual_curation_script.R to fill the `valid` criteria. Data within a well are considered valid is the linear fit on the red data points (figure D), is acceptable.

Once the manual curation has been done, raw_data can be saved into a csv file. For the inference only the observations for which valid==TRUE will be considered.

```{r}
# Saving data
raw_data %>% 
  readr::write_csv("./gly_curated.csv")
```

Now, the same can be done for other conditions.

# GLUCOSE

Based on the previous results, one can set parameters to focus on relevant part of the kinetics experimentally acquired.

Manual curation is performed below per condition.

```{r}
cond <- "glucose"
```

### Keeping all data

```{r message=FALSE,warning=FALSE}
raw_data <- mydata %>%
  filter(condition==cond) %>% #keeping only one relevant condition
  filter(!grepl("empty",strain)) %>% #filter out empty wells
  filter(max(od)>0.05) %>% #to get rid of wells in which there are no growth at all 
  filter(time_min/60<30) #clearly nothing intersting to see after 30 hours.
```

### Keeping only exponential growth regime

To identify exponential growth, we arbitrarily define a time threshold (to discard experimental noise probably due to condensation at the beginning of the kinetics), `time_min_threshold` and coefficient 'rel' so that only observations that are characterized by corrected_od<rel*maximum_corrected_od are considered.

```{r message=FALSE,warning=FALSE}
#parameters for exponential growth
rel <- 0.6
time_min_threshold <- 60

mydata_exp_growth <- mydata %>% 
  filter(condition==cond) %>% # select only one condition
  filter(!grepl("empty",strain)) %>% #filter out empty wells
  filter(max(od)>0.05) %>% #to get rid of wells in which there are no growth at all 
  filter(time_min>time_min_threshold) %>% #discard beginning of the kinetics
  group_by(well) %>%
  arrange(time_min) %>% # the following is to discard all observations after corrected_od>rel*max_od
  mutate(global_max=max(corrected_od,na.rm=TRUE)) %>%
  mutate(relative_max=rel*global_max) %>% 
  mutate(above_max_time=ifelse(corrected_od>=relative_max,time_min,NA)) %>% 
  mutate(time_limit=min(above_max_time,na.rm=TRUE)) %>% 
  filter(time_min<time_limit) %>% 
  ungroup()
```

Checking results.

```{r message=FALSE,warning=FALSE}
plates <- unique(mydata_exp_growth$plate)

for(i in plates){
  mydata_to_plot <- mydata_exp_growth %>% 
    filter(plate==i)
  
  condition <- unique(mydata_to_plot$condition)
  
  print(
  mydata_to_plot %>%
  ggplot()+
  geom_point(aes(time_min/60,log(corrected_od)))+
  facet_grid(row~column)+ 
  theme(strip.text.x = element_blank())+
  labs(subtitle=paste(i,condition,sep=" ")))
  
  print(
  mydata_to_plot %>%
  ggplot()+
  geom_point(aes(time_min/60,log(corrected_fluo)))+
  facet_grid(row~column)+ 
  theme(strip.text.x = element_blank())+
  labs(subtitle=paste(i,condition,sep=" ")))
  }
```

### Keeping data meant to be used for fluo vs corrected_od inference

We noticed that the linear relationship between fluo and corrected_od seems to extent to early stationary phase. By iterating and looking carefully at the fits, we came up with arbitrary `min_corrected_od` and `max_corrected_od` values where a linear fit is acceptable between fluo and corrected_od. The value of these parameters depends on the condition.

```{r message=FALSE,warning=FALSE}
min_corrected_od <- 0.05
max_corrected_od <- 0.5

raw_data_fluo_od <- mydata %>% 
  filter(between(corrected_od,min_corrected_od,max_corrected_od)) %>% 
  filter(condition==cond) %>% 
  filter(!grepl("empty",strain)) %>% #filter out empty wells
  filter(max(od)>min_corrected_od) %>% #to get rid of wells in which there are no growth at all 
  filter(time_min/60<30)
```

### Performing the manual curation

Run the first loop in manual_curation_script.R to fill the `valid` criteria. Data within a well are considered valid is the linear fit on the red data points (figure D), is acceptable.

Once the manual curation has been done, raw_data can be saved into a csv file. For the inference only the observations for which valid==TRUE will be considered.

```{r}
# Saving data
raw_data %>% 
  readr::write_csv("./glu_curated.csv")
```

Now, the same can be done for other conditions.

# GLUCOSEAA

Based on the previous results, one can set parameters to focus on relevant part of the kinetics experimentally acquired.

Manual curation is performed below per condition.

```{r}
cond <- "glucoseaa"
```

### Keeping all data

```{r message=FALSE,warning=FALSE}
raw_data <- mydata %>%
  filter(condition==cond) %>% #keeping only one relevant condition
  filter(!grepl("empty",strain)) %>% #filter out empty wells
  filter(max(od)>0.05) %>% #to get rid of wells in which there are no growth at all 
  filter(time_min/60<30) #clearly nothing intersting to see after 30 hours.
```

### Keeping only exponential growth regime

To identify exponential growth, we arbitrarily define a time threshold (to discard experimental noise probably due to condensation at the beginning of the kinetics), `time_min_threshold` and coefficient 'rel' so that only observations that are characterized by corrected_od<rel*maximum_corrected_od are considered.

```{r message=FALSE,warning=FALSE}
#parameters for exponential growth
rel <- 0.6
time_min_threshold <- 60

mydata_exp_growth <- mydata %>% 
  filter(condition==cond) %>% # select only one condition
  filter(!grepl("empty",strain)) %>% #filter out empty wells
  filter(max(od)>0.05) %>% #to get rid of wells in which there are no growth at all 
  filter(time_min>time_min_threshold) %>% #discard beginning of the kinetics
  group_by(well) %>%
  arrange(time_min) %>% # the following is to discard all observations after corrected_od>rel*max_od
  mutate(global_max=max(corrected_od,na.rm=TRUE)) %>%
  mutate(relative_max=rel*global_max) %>% 
  mutate(above_max_time=ifelse(corrected_od>=relative_max,time_min,NA)) %>% 
  mutate(time_limit=min(above_max_time,na.rm=TRUE)) %>% 
  filter(time_min<time_limit) %>% 
  ungroup()
```

Checking results.

```{r message=FALSE,warning=FALSE}
plates <- unique(mydata_exp_growth$plate)

for(i in plates){
  mydata_to_plot <- mydata_exp_growth %>% 
    filter(plate==i)
  
  condition <- unique(mydata_to_plot$condition)
  
  print(
  mydata_to_plot %>%
  ggplot()+
  geom_point(aes(time_min/60,log(corrected_od)))+
  facet_grid(row~column)+ 
  theme(strip.text.x = element_blank())+
  labs(subtitle=paste(i,condition,sep=" ")))
  
  print(
  mydata_to_plot %>%
  ggplot()+
  geom_point(aes(time_min/60,log(corrected_fluo)))+
  facet_grid(row~column)+ 
  theme(strip.text.x = element_blank())+
  labs(subtitle=paste(i,condition,sep=" ")))
  }
```

### Keeping data meant to be used for fluo vs corrected_od inference

We noticed that the linear relationship between fluo and corrected_od seems to extent to early stationary phase. By iterating and looking carefully at the fits, we came up with arbitrary `min_corrected_od` and `max_corrected_od` values where a linear fit is acceptable between fluo and corrected_od. The value of these parameters depends on the condition.

```{r message=FALSE,warning=FALSE}
min_corrected_od <- 0.05
max_corrected_od <- 0.75

raw_data_fluo_od <- mydata %>% 
  filter(between(corrected_od,min_corrected_od,max_corrected_od)) %>% 
  filter(condition==cond) %>% 
  filter(!grepl("empty",strain)) %>% #filter out empty wells
  filter(max(od)>min_corrected_od) %>% #to get rid of wells in which there are no growth at all 
  filter(time_min/60<30)
```

### Performing the manual curation

Run the first loop in manual_curation_script.R to fill the `valid` criteria. Data within a well are considered valid is the linear fit on the red data points (figure D), is acceptable.

Once the manual curation has been done, raw_data can be saved into a csv file. For the inference only the observations for which valid==TRUE will be considered.

```{r}
# Saving data
raw_data %>% 
  readr::write_csv("./gluaa_curated.csv")
```

End of the manual curation!
