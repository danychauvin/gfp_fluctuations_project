# Importing necessary packages
source("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/mother_machine_experiments_toolbox/load_functions_and_packages.R")
#Denoising output path
denoising_folder <- "/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoising_raw_data/denoising_20230512104746"
correlation_folder <- "correlations_results"


# Preparing and running computation on the cluster
# Parameters are saved in .txt files, and data are saved in .csv files
list_of_data_files <- lapply(denoising_folder,function(.l){
  list.files(path=.l,pattern="joints.csv$",recursive = TRUE,full.names=TRUE)}) %>% 
  unlist()

df_of_data_files <- tibble(path=list_of_data_files) %>% 
    tidyr::extract(path,"promoter","[0-9a-zA-Z]{1,}_([0-9a-zA-Z]{1,})_[0-9]{8}_rawdata_f01234578910_b_joints.csv$",remove=FALSE,convert=FALSE) %>% 
    tidyr::extract(path,"condition","([0-9a-zA-Z]{1,})_[0-9a-zA-Z]{1,}_[0-9]{8}_rawdata_f01234578910_b_joints.csv$",remove=FALSE,convert=FALSE) %>%
    tidyr::extract(path,"date","[0-9a-zA-Z]{1,}_[0-9a-zA-Z]{1,}_([0-9]{8})_rawdata_f01234578910_b_joints.csv$",remove=FALSE,convert=FALSE) %>% 
  rename(data_path=path)

dt_df <- tibble(condition=c("acetate005","glycerol040","glucose020","glucoseaa020"),dt=c("12","6","3","1.5"))

code_path <- "/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/gfp_gaussian_process/python_src/correlation_from_joint.py"

df_of_experimental_files <- df_of_data_files %>% 
  left_join(dt_df,by=c("condition")) %>% 
  distinct(data_path,condition,dt) %>% 
  mutate(cmd=sprintf("python %s -d %s -k %s -dt %s -delimiter _",code_path,data_path,condition,dt))

output_path <- sprintf("%s/%s",denoising_folder,correlation_folder)
system(sprintf("mkdir %s/%s &> /dev/null",denoising_folder,correlation_folder))
job_template <- "/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/gfp_fluctuations_project/denoising_raw_data/inference_input_all/jobs_template.sh"

command_path <- paste(output_path,"/compute_correlation_functions.cmd",sep="")
results_path <- sprintf("%s/results_%%a",output_path)
errors_path <- sprintf("%s/errors_%%a",output_path)
system(sprintf("cp %s %s/jobs.sh",job_template,output_path))
jobs_path <- sprintf("%s/jobs.sh",output_path)
system(sprintf("sed -i \"s+path_to_output+%s+g\" %s",results_path,jobs_path))
system(sprintf("sed -i \"s+path_to_error+%s+g\" %s",errors_path,jobs_path))
system(sprintf("sed -i \"s+path_to_commands+%s+g\" %s",command_path,jobs_path))

#Writting commands.cmd
write_csv(df_of_experimental_files %>% 
            select(cmd),command_path,col_names = FALSE)

#The following line has be sent to login20...
runAllJobsLine <- sprintf("sbatch --array=1-$(cat %s|wc -l):1 jobs.sh",command_path)
print(runAllJobsLine)
system(sprintf("cd %s && %s",output_path,runAllJobsLine))