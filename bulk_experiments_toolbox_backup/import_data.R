# Import necessary packages and functions
# Dany 20210204


myconditions <- readr::read_csv(path_to_data_list,
                                col_types = cols(
                                  date=col_character(),
                                  description=col_character(),
                                  data_path=col_character(),
                                  plate_map_path=col_character()))

# Import the data

mydata <- myconditions %>%
  group_by(date,description,plate) %>% #date here is as per datalist, this is not the date of the measurement, description and plate should refer, together, to a unique experiment
    do((function(.df){
    plate_ <- as.character(unique(.df$plate))
    date_ <- as.character(unique(.df$date))
    description_ <- as.character(unique(.df$description))
    path_ <- as.character(unique(.df$data_path))
    new_df <- list.files(path_,sprintf("%s_%s_[0-9]{8}_[0-9]{6}_FE_.txt$",description_,plate_),full.names = TRUE) %>% 
      unlist() %>% 
      lapply(read_Biotek_Synergy2_kinetic) %>% 
      #extract from the name the time
      bind_rows() %>% 
      mutate(date=date_,
             description=description_,
             plate=plate_)
    return(new_df)})(.)) %>% 
  select(date,description,plate,row,column,od,fluo,measurement_date,measurement_time,fluo_channel,od_channel) %>% 
  ungroup() %>% 
  rename(exp_start_date=date)
  
# Import plate layout
  
mylayouts <- myconditions %>%
  group_by(date,description,plate) %>% 
  do((function(.df){
    plate_ <- as.character(unique(.df$plate))
    date_ <- as.character(unique(.df$date))
    description_ <- as.character(unique(.df$description))
    path_ <- as.character(unique(.df$plate_map_path))
    new_df <- read_plate_layout(path_) %>% 
      bind_rows() %>% 
      mutate(date=date_,
             description=description_,
             plate=plate_)
    return(new_df)})(.)) %>% 
  ungroup() %>% 
  rename(exp_start_date=date)

# Completed data

mydata <- mydata %>% 
  left_join(mylayouts,by=c("exp_start_date","plate","description","row","column"))

# Convert measurement date and time into universal time, using posix format
# "date" is the date of the measurement
mydata <- mydata %>% 
  mutate(measurement_date_posix=paste(substr(measurement_date,1,4),substr(measurement_date,5,6),substr(measurement_date,7,8),sep="/")) %>% 
  mutate(measurement_time_posix=paste(substr(measurement_time,1,2),substr(measurement_time,3,4),substr(measurement_time,5,6),sep=":")) %>%
  mutate(dateTime=as.POSIXct(paste(measurement_date_posix,measurement_time_posix,sep=" "))) %>% 
  group_by(description,plate) %>% 
  mutate(time_min=difftime(dateTime,min(dateTime),units="min"),
         time_min=as.double(time_min)) %>% 
  ungroup()

# Convert OD and fluorescence to numbers, add a fluo_sat column
mydata <- mydata %>% 
  mutate(fluo_sat=FALSE) %>% 
  mutate(fluo_sat=ifelse(fluo_channel=="OVRFLW",TRUE,fluo_sat),
         fluo=ifelse(fluo=="OVRFLW",NA,fluo)) %>% 
  mutate(od=as.double(od)) %>% 
  mutate(fluo=as.double(fluo))

# Data should be ready there

