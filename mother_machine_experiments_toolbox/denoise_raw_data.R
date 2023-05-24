# Load necessary functions and packages
source("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/mother_machine_experiments_toolbox/load_functions_and_packages.R")
# Define where raw data are stored
data_folder <- c("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/raw_curated_data_complete_cycles")
# Output folder of the denoising procedure
denoising_folder <- c("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoising_raw_data")
#Path to the folder where parameters initial parameters are stored
parameter_folder<-c("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoising_raw_data/parameters_with_bleaching")
#Label for the denoising
label <- "denoising"

# Generating time stamp for the inference procedure
time_stamp <- c(year(now()),month(now()),day(now()),hour(now()),minute(now()),round(second(now()),0))
time_stamp <- unlist(lapply(time_stamp,function(.l){
  .new_l <- as.character(.l)
  if(str_length(.new_l)==1){
    .new_l <- paste("0",.new_l,sep="")}
  return(.new_l)}))
time_stamp <- paste(time_stamp,sep="",collapse="")
inference_id <- sprintf("%s_%s",label,time_stamp)

# Define necessary variables and start denoising
folder_name <- 'denoising_raw_data' #folder in which results are stored
root_folder <- denoising_folder

codePath <- "/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/gfp_gaussian_process/bin/gfp_gaussian"
tolerance <- "1e-40"

inference_input_path_all <- paste(root_folder,"/inference_input_all",sep="")
csv_config_path <- paste(inference_input_path_all,"/csv_config.txt",sep="")
system(sprintf("mkdir %s/%s 2> /dev/null",root_folder,inference_id))

##Take the folder content as the input to run the inference.
output_path <- paste(root_folder,inference_id,sep="/")

##parameters are saved in .txt files, and data are saved in .csv files
list_of_data_files <- lapply(data_folder,function(.l){
  list.files(path=.l,pattern=".csv$",recursive = TRUE,full.names=TRUE)}) %>% 
  unlist()

df_of_data_files <- tibble(path=list_of_data_files) %>% 
    tidyr::extract(path,"promoter","[0-9a-zA-Z]{1,}_([0-9a-zA-Z]{1,})_[0-9]{8}_rawdata.csv$",remove=FALSE,convert=FALSE) %>% 
    tidyr::extract(path,"condition","([0-9a-zA-Z]{1,})_[0-9a-zA-Z]{1,}_[0-9]{8}_rawdata.csv$",remove=FALSE,convert=FALSE) %>%
    tidyr::extract(path,"date","[0-9a-zA-Z]{1,}_[0-9a-zA-Z]{1,}_([0-9]{8})_rawdata.csv$",remove=FALSE,convert=FALSE) %>% 
  rename(data_path=path)

list_of_parameter_files<- lapply(parameter_folder,function(.l){
  list.files(path=.l,pattern=".txt$",recursive = TRUE,full.names=TRUE)}) %>% 
  unlist()

df_of_parameter_files <- tibble(path=list_of_parameter_files) %>% 
    tidyr::extract(path,"promoter","[a-z0-9A-Z]{1,}_([0-9a-zA-Z]{1,})[_filtered]{0,}_parameters.txt$",remove=FALSE,convert=FALSE) %>% 
    tidyr::extract(path,"condition","([a-z0-9A-Z]{1,})_[0-9a-zA-Z]{1,}[_filtered]{0,}_parameters.txt$",remove=FALSE,convert=FALSE) %>% 
  filter(!(grepl("csv_config",path))) %>% 
  rename(parameter_path=path)

df_of_experimental_files <- left_join(df_of_data_files,df_of_parameter_files,by=c("condition","promoter")) %>% 
  select(condition,promoter,data_path,parameter_path) %>% 
  mutate(cmd=sprintf("%s -i %s -b %s -c %s -o %s/ -m -p -l 1 -t %s -space log -j -div binomial",codePath,data_path,parameter_path,csv_config_path,output_path,as.character(tolerance)))

command_path <- sprintf("%s/commands.cmd",output_path) 
results_path <- sprintf("%s/results_%%a",output_path)
errors_path <- sprintf("%s/errors_%%a",output_path)
system(sprintf("cp %s/jobs_template.sh %s/jobs.sh",paste(root_folder,"/inference_input_all",sep=""),output_path))
jobs_path <- sprintf("%s/jobs.sh",output_path)
system(sprintf("sed -i \"s+path_to_output+%s+g\" %s",results_path,jobs_path))
system(sprintf("sed -i \"s+path_to_error+%s+g\" %s",errors_path,jobs_path))
system(sprintf("sed -i \"s+path_to_commands+%s+g\" %s",command_path,jobs_path))

##Writting commands.cmd
write_csv(df_of_experimental_files %>% 
            select(cmd),command_path,col_names = FALSE)

#The following line has be sent to login20...
runAllJobsLine <- sprintf("sbatch --array=1-$(cat %s|wc -l):1 jobs.sh",command_path)
print(runAllJobsLine)

system(sprintf("cd /scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoising_raw_data/%s_%s && %s",label,time_stamp,runAllJobsLine))