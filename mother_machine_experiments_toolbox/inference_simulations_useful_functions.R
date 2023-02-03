generate_subsampled_data <- function(.myframes_to_mycells,exported_data_folder_var,step_var,.type){
  myframes_to_mycells_sub <- .myframes_to_mycells %>% 
    filter(type==.type) %>% 
    semi_join(mydataset,by=c("date","promoter","vector","condition","type")) %>% 
    mutate(exp_reference=paste(date,condition,promoter,.type,sep=".")) %>% 
    group_by(exp_reference) %>% 
    arrange(frame) %>%
    mutate(sub_cell=NA) %>%
    do((function(.df){
      new_df <- .df %>% 
        mutate(sub_cell=NA) %>% 
        mutate(sub_parent=NA) 
      for(i in 0:(step_var-1)){
        new_df <- new_df %>% 
          mutate(sub_cell=ifelse(frame%%step_var==i,paste(cell,step_var,i,sep="."),sub_cell)) %>% 
          mutate(sub_parent=ifelse(frame%%step_var==i,paste(parent,step_var,i,sep="."),sub_parent))}
      return(new_df)})(.)) %>%
    ungroup()
  
  # Export the data in the current form
  myframes_to_mycells_sub %>% 
    ungroup() %>% 
    group_by(condition,promoter,type) %>% 
    arrange(sub_cell,time_sec) %>% 
    do((function(.l){
      readr::write_csv(.l,sprintf("%s/%s_%s.csv",exported_data_folder_var,unique(.l$condition),unique(.l$promoter)))
      return(data_frame())})(.)) %>%
    ungroup()}

generate_initial_parameter_df <- function(){
  lambda_df <- mycells %>%
    semi_join(mydataset,by=c("date","vector","condition","promoter")) %>% 
    group_by(condition,promoter) %>% 
    distinct(cell,.keep_all = TRUE) %>% 
    mutate(mean_lambda_full=mean(alpha*60)) %>% #mean_lambda should be /min
    mutate(var_lambda_full=var(alpha*60)) %>% #var_lambda should be in /min2
    mutate(mean_q_full=mean(q*60,na.rm=TRUE)) %>% 
    mutate(var_q_full=var(q*60,na.rm=TRUE)) %>% #per min
    mutate(mean_div_time=mean(div_time/60)) %>% 
    mutate(gamma_lambda=1/mean_div_time) %>% 
    mutate(mean_dl=mean(dl,na.rm=TRUE)) %>% 
    ungroup %>% 
    distinct(condition,promoter,.keep_all=TRUE) %>% 
    select(condition,promoter,mean_lambda_full,var_lambda_full,gamma_lambda,var_q_full,mean_q_full,mean_dl) %>% 
    arrange(condition,promoter)
  
  mycells_alpha <- mycells %>% 
    semi_join(mydataset,by=c("date","vector","condition","promoter")) %>% 
    distinct(cell,alpha)
  
  q_x_df <- myframes_to_mycells %>% 
    semi_join(mydataset,by=c("date","vector","condition","promoter")) %>% 
    left_join(mycells_alpha,by=c("cell")) %>%
    group_by(condition,promoter) %>% 
    mutate(res=log(length_predict)-log(length_um)) %>% 
    mutate(var_x_res=var(res,na.rm=TRUE)) %>% 
    mutate(mean_c=mean(gfp_nb/length_um)) %>% 
    mutate(var_x=(mean(log(length_um))/100)**2) %>% 
    mutate(var_g=(mean(gfp_nb)/100)**2) %>% 
    ungroup() %>% 
    group_by(cell) %>% 
    mutate(lag_gfp_nb=lag(gfp_nb),
           lag_length_um=lag(length_um),
           dg=gfp_nb-lag_gfp_nb,
           dgdt=dg/t_interval,
           q_cell=dgdt*1/length_um,
           mean_q_cell=mean(q_cell,na.rm=TRUE)) %>% #per min
    ungroup() %>% 
    group_by(promoter,condition) %>% 
    mutate(var_q_all=var(mean_q_cell,na.rm=TRUE)) %>% 
    mutate(mean_q_all=mean(mean_q_cell,na.rm=TRUE),
           var_q_square=mean_q_all**2) %>%
    ungroup() %>% 
    distinct(condition,promoter,.keep_all=TRUE) %>% 
    select(condition,promoter,mean_q_all,var_q_all,var_q_square,var_g,var_x,mean_c,var_x_res,t_interval)
  
  # Adding to the parameter, information about asymmetry at division.
  daughter_mother_df <- myframes_to_mycells %>% 
    semi_join(mydataset,by=c("date","vector","condition","promoter")) %>% 
    distinct(cell,parent)
  
  parent_df <- mycells %>% 
    semi_join(mydataset,by=c("date","vector","condition","promoter")) %>% 
    left_join(daughter_mother_df,by=c("cell")) %>% 
    select(condition,promoter,cell,l_div,g_div) %>% 
    rename(l_div_mother=l_div,
           g_div_mother=g_div,
           parent=cell)
  
  asymetry_df <- mycells %>% 
    semi_join(mydataset,by=c("date","vector","condition","promoter")) %>% 
    left_join(daughter_mother_df,by=c("cell")) %>% 
    left_join(parent_df,by=c("parent","condition","promoter")) %>% #some parent in mycells might not be in cell. Hence some NA's
    select(condition,promoter,cell,parent,l_div_mother,g_div_mother,l_birth,g_birth) %>% 
    drop_na() %>% 
    group_by(condition,promoter) %>% 
    mutate(var_dx=var(-log(l_div_mother)+log(l_birth)+2,na.rm=TRUE)) %>% 
    mutate(var_dg=var(g_birth-g_div_mother/2),na.rm=TRUE) %>% 
    ungroup() %>% 
    distinct(condition,promoter,.keep_all = TRUE) %>% 
    select(condition,promoter,var_dx,var_dg)
  
  final_df <- left_join(lambda_df,q_x_df,by=c("condition","promoter")) %>% 
    left_join(asymetry_df,by=c("condition","promoter")) %>% 
    mutate(beta=0) %>% # Beta = 0 is hardcoded.
    mutate(gamma_q=gamma_lambda)
  
  return(final_df)}

generating_parameter_files <- function(replica,fuzzy_parameters){
  
  to_vector <- function(m,n){
    return(signif(c(m,m,m/n,m*n),2))
  }
  
  to_vector_beta <- function(m){
    return(signif(c(m,m,1e-10,1e-2),2))
  }
  
  to_vector_free <- function(m){
    return(signif(c(m,m),2))
  }
  
  to_vector_fixed <- function(m){
    return(signif(c(m),2))
  }
  
  param_df %>% 
    #filter(promoter %in% promoters) %>% 
    #filter(condition %in% conditions) %>% 
    group_by(condition,promoter) %>% 
    do((function(.df){
      if(fuzzy_parameters==TRUE){
        fuzzy_parameters_val <- runif(11,min=-0.5,max=0.5)
      }else{fuzzy_parameters_val <- rep(c(0),11)}
      condition <- .df$condition
      promoter <- .df$promoter
      mean_lambda <- to_vector_free(.df$mean_lambda_full) #to_vector(.df$mean_lambda_full*(1+fuzzy_parameters_val[1]),1.5) #.df$mean_lambda_full
      mean_q <- to_vector_free(.df$mean_q_full) #to_vector(.df$mean_q_full*(1+fuzzy_parameters_val[2]),1.5) #.df$mean_q_full
      #if(condition=="ace"){
        #gamma_lambda <- to_vector_free(.df$gamma_lambda*0.5*(1+fuzzy_parameters_val[3])) #In acetate the autocorrelation time is larger
      #}else{
      gamma_lambda <- to_vector_free(.df$gamma_lambda*1*(1+fuzzy_parameters_val[3]))#}
      gamma_q <- to_vector_free(.df$gamma_q*5*(1+fuzzy_parameters_val[4]))
      
      #if(condition=="ace"){
      #  var_lambda_val <- 2*(.df$var_lambda_full)*(.df$gamma_lambda*0.5)
      #}else{
      var_lambda_val <- 2*(.df$var_lambda_full)*(.df$gamma_lambda*1)#}
      
      var_q_val <- 2*(.df$var_q_all)*(.df$gamma_q)*5
      var_lambda <- to_vector_free(var_lambda_val*(1+fuzzy_parameters_val[5]))
      var_q <- to_vector_free(var_q_val*(1+fuzzy_parameters_val[6]))
      var_x <- to_vector_free(.df$var_x_res*(1+fuzzy_parameters_val[7]))
      var_g <- to_vector_free(.df$var_g*(1+fuzzy_parameters_val[8]))
      beta <- to_vector_fixed(.df$beta)#*(1+fuzzy_parameters_val[9]))
      var_dx <- to_vector_free(.df$var_dx*(1+fuzzy_parameters_val[10]))
      var_dg <- to_vector_free(.df$var_dg*(1+fuzzy_parameters_val[11]))
      
      input_path_rep <- inference_input_path #paste(inference_input_path,replica,sep="_")
      output_path_rep <- inference_output_path #paste(inference_output_path,replica,sep="_")
      #system(sprintf('mkdir %s',input_path_rep))
      #system(sprintf('mkdir %s',output_path_rep))
      system(sprintf('touch %s/%s_%s_parameters.txt',input_path_rep,condition,promoter))
      system(sprintf('touch %s/%s_%s_parameters.txt',input_path_rep,condition,promoter))
      cat("mean_lambda = ",paste(mean_lambda,collapse=", "),"\n",file=sprintf('%s/%s_%s_parameters.txt',input_path_rep,condition,promoter),sep=' ',append=TRUE)
      cat("gamma_lambda = ",paste(gamma_lambda,collapse=", "),"\n",file=sprintf('%s/%s_%s_parameters.txt',input_path_rep,condition,promoter),sep='',append=TRUE)
      cat("var_lambda = ",paste(var_lambda,collapse=", "),"\n",file=sprintf('%s/%s_%s_parameters.txt',input_path_rep,condition,promoter),sep='',append=TRUE)
      cat("mean_q = ",paste(mean_q,collapse=", "),"\n",file=sprintf('%s/%s_%s_parameters.txt',input_path_rep,condition,promoter),sep='',append=TRUE)
      cat("gamma_q = ",paste(gamma_q,collapse=", "),"\n",file=sprintf('%s/%s_%s_parameters.txt',input_path_rep,condition,promoter),sep='',append=TRUE)
      cat("var_q = ",paste(var_q,collapse=", "),"\n",file=sprintf('%s/%s_%s_parameters.txt',input_path_rep,condition,promoter),sep='',append=TRUE)
      cat("var_x = ",paste(var_x,collapse=", "),"\n",file=sprintf('%s/%s_%s_parameters.txt',input_path_rep,condition,promoter),sep='',append=TRUE)
      cat("var_g = ",paste(var_g,collapse=", "),"\n",file=sprintf('%s/%s_%s_parameters.txt',input_path_rep,condition,promoter),sep='',append=TRUE)
      cat("beta = ",paste(beta,collapse=", "),"\n",file=sprintf('%s/%s_%s_parameters.txt',input_path_rep,condition,promoter),sep='',append=TRUE)
      cat("var_dx = ",paste(var_dx,collapse=", "),"\n",file=sprintf('%s/%s_%s_parameters.txt',input_path_rep,condition,promoter),sep='',append=TRUE)
      cat("var_dg = ",paste(var_dg,collapse=", "),"\n",file=sprintf('%s/%s_%s_parameters.txt',input_path_rep,condition,promoter),sep='',append=TRUE)
      return(data_frame())
    })(.)) %>% 
    ungroup()}

generating_parameter_files_simulation_from_tidied <- function(.myparameters,.replica,.type){
  
  names <- c("mean_lambda","gamma_lambda","var_lambda","mean_q","gamma_q","var_q","beta","var_x","var_g","var_dx","var_dg")
  
  .myparameters_new <- .myparameters %>%
    tidyr::extract(output_directory,"type","/outputs_([a-z]{1,})$",remove=FALSE,convert=FALSE) %>%
    filter(type==.type) %>%
    select(-c(type,output_directory)) %>% 
    pivot_longer(cols=starts_with(c("mean","gamma","var","beta","dt_m","dv"))) %>%
    filter(name %in% names) %>% 
    rename(final=value)
  
  .myparameters_new$name <- factor(.myparameters_new$name,levels=names)
  .myparameters_new <- arrange(.myparameters_new,by=name) %>% 
    group_by(condition,promoter) %>% 
    do((function(.df){
      
      condition <- unique(.df$condition)
      promoter <- unique(.df$promoter)
      
      .new_df <- .df %>% 
        mutate(equal="=") %>% 
        select(name,equal,final) %>% 
        mutate(step=final/10) %>%
        rename(init=final) %>% 
        mutate(command=sprintf("%s = %s, %s",name,as.character(init),as.character(step)),
               command=ifelse(grepl("= 0, 0",command),sprintf("%s = %s",name,as.character(init)),command)) %>% 
        select(command)
      
      input_path_rep <- inference_input_path #paste(inference_input_path,.replica,sep="_")
      output_path_rep <- inference_output_path #paste(inference_output_path,.replica,sep="_")
      #system(sprintf('mkdir %s',input_path_rep))
      #system(sprintf('mkdir %s',output_path_rep))
      system(sprintf('touch %s/%s_%s_parameters.txt',input_path_rep,condition,promoter))
      #system(sprintf('touch %s/%s_%s_parameters.txt',input_path_rep,condition,promoter))
      readr::write_delim(.new_df,sprintf('%s/%s_%s_parameters.txt',input_path_rep,condition,promoter),delim="",col_names=FALSE)
      
      return(data_frame())
    })(.)) %>% 
    ungroup()}

generating_parameter_files_simulation <- function(.myparameters,.replica,.type){
  .myparameters %>% 
    group_by(condition,promoter) %>% 
    do((function(.df){
      
      condition <- unique(.df$condition)
      promoter <- unique(.df$promoter)
      
      .new_df <- .df %>% 
        mutate(equal="=") %>% 
        select(name,equal,final,step) %>% 
        mutate(step=final/10) %>% 
        rename(init=final) %>% 
        mutate(step=ifelse(name=="beta","",as.character(step))) %>% 
        mutate(command=sprintf("%s = %s, %s",name,as.character(init),as.character(step)),
               command=ifelse(name=="beta",sprintf("%s = %s",name,as.character(init)),command)) %>% 
        select(command)
        
      input_path_rep <- inference_input_path #paste(inference_input_path,.replica,sep="_")
      output_path_rep <- inference_output_path #paste(inference_output_path,.replica,sep="_")
      #system(sprintf('mkdir %s',input_path_rep))
      #system(sprintf('mkdir %s',output_path_rep))
      system(sprintf('touch %s/%s_%s_parameters.txt',input_path_rep,condition,promoter))
      system(sprintf('touch %s/%s_%s_parameters.txt',input_path_rep,condition,promoter))
      readr::write_delim(.new_df,sprintf('%s/%s_%s_parameters.txt',input_path_rep,condition,promoter),delim="",col_names=FALSE)
      
      return(data_frame())
    })(.)) %>% 
    ungroup()}

prepare_for_slurm <- function(replica,tolerance){
  
  #temporary and ugly
  conditions <- as.character(unique(mydataset$condition))
  promoters <- as.character(unique(mydataset$promoter))
  
  input_path_rep <- inference_input_path #paste(inference_input_path,replica,sep="_")
  output_path_rep <- inference_output_path #paste(inference_output_path,replica,sep="_")
  ref_tibble <- tibble()
  k <- 0
  for(i in 1:length(conditions)){
    for(j in 1:length(promoters)){
      #Append jobs to .input/commands.cmd
      lineToWrite <- sprintf("%s -i %s/%s -b %s -c %s -o %s/ -m -p -l 1 -t %s -space log -j",
                             codePath,
                             exported_data_folder,
                             sprintf("%s_%s.csv",conditions[i],promoters[j]),
                             sprintf("%s/%s_%s_parameters.txt",input_path_rep,conditions[i],promoters[j]),
                             csv_config_path,
                             output_path_rep,
                             tolerance)
      
      system(sprintf('touch %s/commands.cmd',input_path_rep))
      write(lineToWrite,file=sprintf('%s/commands.cmd',input_path_rep),append=TRUE)
      k <- k+1
      cond_vec <- list(conditions[i],promoters[j],k)
      ref_tibble <- rbind(ref_tibble,cond_vec)
    }}
  #Modify jobs.sh
  command_path <- sprintf("%s/commands.cmd",input_path_rep) 
  results_path <- sprintf("%s/results_%%a",output_path_rep)
  errors_path <- sprintf("%s/errors_%%a",output_path_rep)
  system(sprintf("cp %s/jobs_template.sh %s/jobs.sh",inference_input_path_all,input_path_rep))
  jobs_path <- sprintf("%s/jobs.sh",input_path_rep)
  
  system(sprintf("sed -i \"s+path_to_output+%s+g\" %s",results_path,jobs_path))
  system(sprintf("sed -i \"s+path_to_error+%s+g\" %s",errors_path,jobs_path))
  system(sprintf("sed -i \"s+path_to_commands+%s+g\" %s",command_path,jobs_path))
  
  #The following line has be sent to login20...
  runAllJobsLine <- sprintf("sbatch --array=1-$(cat %s|wc -l):1 jobs.sh",command_path)
  print(runAllJobsLine)}
# Here one can open a terminal within R studio, ssh login to login20
# ssh rocasu25@login20.scicore.unibas.ch
# cd folderPath/optimization_code
# copy paste runAllJobsLine to the terminal
# jobs should be submitted to the slurm scheduler


recover_nb_points_single <- function(fileName){

  nfileName <- as.character(unique(fileName))
  nlines <- length(readr::read_lines(nfileName))
  if(nlines<20){return(NA)}else{
  nb_points <- readr::read_lines(nfileName) %>% 
             .[[20]] %>% 
             str_split(pattern=",") %>% 
             unlist() %>% 
             .[[2]] %>% 
             as.double() %>% rep(times=length(fileName))
  return(nb_points)}}

recover_log_likelihood_single <- function(fileName){

  nfileName <- as.character(unique(fileName))
  nlines <- length(readr::read_lines(nfileName))
  if(nlines<20){return(NA)}else{
  log_likelihood <- readr::read_lines(nfileName) %>% 
             .[[21]] %>% 
             str_split(pattern=",") %>% 
             unlist() %>% 
             .[[2]] %>% 
             as.double() %>% rep(times=length(fileName))
  return(log_likelihood)}}


recover_final_parameters <- function(folderName){
  fileList <- list.files(folderName,pattern='*_final.csv',full.names = TRUE)
  #nameList <- list.files(folderName,pattern='*_final.csv',full.names = FALSE)
  #conditions <- lapply(nameList,function(.l) str_split(.l,"_")[[1]][1])
  #promoters <- lapply(nameList,function(.l) str_split(.l,"_")[[1]][2])
  if(length(fileList)==0){
    return(tibble())
  }else{
  df <- fileList %>% 
    lapply(function(.l) readr::read_csv(.l,col_names=TRUE,n_max=11) %>% 
             mutate(filepath=.l) %>%
             mutate(nb_points= recover_nb_points_single(filepath)) %>%
             mutate(log_likelihood= recover_log_likelihood_single(filepath))) %>%
    do.call(rbind,.) %>% 
    tidyr::extract(filepath,"replica",".*/inference_output_([1-9]{1})/[a-z]{1,}_[a-z0-9A-Z]{1,}_[a-z0-9]{1,}_[a-z0-9]{1,}_final.csv$",convert=TRUE,remove=FALSE) %>% 
    tidyr::extract(filepath,"condition",".*/inference_output_[1-9]{1}/([a-z]{1,})_[a-z0-9A-Z]{1,}_[a-z0-9]{1,}_[a-z0-9]{1,}_final.csv$",convert=TRUE,remove=FALSE) %>% 
    tidyr::extract(filepath,"promoter",".*/inference_output_[1-9]{1}/[a-z]{1,}_([a-z0-9A-Z]{1,})_[a-z0-9]{1,}_[a-z0-9]{1,}_final.csv$",convert=TRUE,remove=FALSE)
  return(df)}
}

#DEBUG
# Multiple and identical initial conditions in the inference results df. Why is that?
#inference_results_df %>% filter(condition=="glycerol",promoter=="rrnB",name=="gamma_q",replica==5)

#DEBUG
#list_of_replica_df <- lapply(folderNames, recover_final_parameters)
#recover_final_parameters("/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/constitexpr/20210224_constitexpr_inference/inference_output_5")
#t <- "/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/constitexpr/20210224_constitexpr_inference/inference_output_5"
#fileList <- list.files(t,pattern='*_final.csv',full.names = TRUE)
#for(e in fileList){
#  print(e)
#  recover_nb_points_single(e)}
# For a reason that I do not get, rpmb glycerol has final, but no last output... which is due, I think to the fact that it was building that after 6 hours.
###

recover_init_parameters <- function(folderName){
  fileList <- list.files(folderName,pattern='*parameters.txt',full.names = TRUE)
  #nameList <- list.files(folderName,pattern='*_final.csv',full.names = FALSE)
  #conditions <- lapply(nameList,function(.l) str_split(.l,"_")[[1]][1])
  #promoters <- lapply(nameList,function(.l) str_split(.l,"_")[[1]][2])
  if(length(fileList)==0){
    return(tibble())
  }else{
    df <- fileList %>% 
      lapply(function(.l) convert_text(.l)) %>% 
      bind_rows() %>% 
        tidyr::extract(filepath,"replica",".*/inference_input_([1-9]{1})/[a-z]{1,}_[a-z0-9A-Z]{1,}_parameters.txt$",convert=TRUE,remove=FALSE) %>% 
        tidyr::extract(filepath,"condition",".*/inference_input_[1-9]{1}/([a-z]{1,})_[a-z0-9A-Z]{1,}_parameters.txt$",convert=TRUE,remove=FALSE) %>% 
        tidyr::extract(filepath,"promoter",".*/inference_input_[1-9]{1}/[a-z]{1,}_([a-z0-9A-Z]{1,})_parameters.txt$",convert=TRUE,remove=FALSE)
    return(df)}
}

recover_init_parameters_from_final <- function(folderName){
  #fileList <- list.files(folderName,pattern='*parameters.txt',full.names = TRUE)
  #nameList <- list.files(folderName,pattern='*_final.csv',full.names = FALSE)
  #conditions <- lapply(nameList,function(.l) str_split(.l,"_")[[1]][1])
  #promoters <- lapply(nameList,function(.l) str_split(.l,"_")[[1]][2])
  if(length(fileList)==0){
    return(tibble())
  }else{
    df <- fileList %>% 
      lapply(function(.l) convert_text(.l)) %>% 
      bind_rows() %>% 
      tidyr::extract(filepath,"replica",".*/inference_input_([1-9]{1})/[a-z]{1,}_[a-z0-9A-Z]{1,}_parameters.txt$",convert=TRUE,remove=FALSE) %>% 
      tidyr::extract(filepath,"condition",".*/inference_input_[1-9]{1}/([a-z]{1,})_[a-z0-9A-Z]{1,}_parameters.txt$",convert=TRUE,remove=FALSE) %>% 
      tidyr::extract(filepath,"promoter",".*/inference_input_[1-9]{1}/[a-z]{1,}_([a-z0-9A-Z]{1,})_parameters.txt$",convert=TRUE,remove=FALSE)
    return(df)}
}


extract_val <- function(.l,.path_to_text){ 
  .new_l <- unlist(str_split(str_replace_all(.l," +",""),pattern="="))
  #print(.new_l)
  .name <- .new_l[1]
  #print(.name)
  .val <- as.double(unlist(str_split(.new_l[2],pattern=","))[1])
  .df <- tibble(name=.name,init=as.double(.val),filepath=.path_to_text)
  return(.df)}

convert_text <- function(.path_to_text){
  .t <- readr::read_lines(.path_to_text) %>%
    lapply(function(.l) extract_val(.l,.path_to_text)) %>% 
    bind_rows()
  return(.t)}

#DEBUG
#l <- '/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/constitexpr/20210216_constitexpr_inference/inference_input_1'
#g <- '/scicore/home/nimwegen/rocasu25/Documents/Projects/biozentrum/constitexpr/20210216_constitexpr_inference/inference_input_1/glucoseaa_hi1_parameters.txt'
#t <- readr::read_lines(l)

#fileList <- list.files(l,pattern='*parameters.txt',full.names = TRUE)
#df <- fileList %>% 
#  lapply(function(.l) convert_text(.l))

#convert_text(g)
#



recover_final_parameters_errors <- function(folderName){
  fileList <- list.files(folderName,pattern='*_final.csv',full.names = TRUE)
  #nameList <- list.files(folderName,pattern='*_final.csv',full.names = FALSE)
  #conditions <- lapply(nameList,function(.l) str_split(.l,"_")[[1]][1])
  #promoters <- lapply(nameList,function(.l) str_split(.l,"_")[[1]][2])
  if(length(fileList)==0){
    return(tibble())
  }else{
  df <- fileList %>% 
    lapply(function(.l) readr::read_csv(.l,skip=14,n_max=3) %>% mutate(filepath=.l)) %>% 
    do.call(rbind,.) %>% 
    tidyr::extract(filepath,"replica",".*/inference_output_([1-9]{1})/[a-z]{1,}_[a-z0-9A-Z]{1,}_[a-z0-9]{1,}_[a-z0-9]{1,}_final.csv$",convert=TRUE,remove=FALSE) %>% 
    tidyr::extract(filepath,"condition",".*/inference_output_[1-9]{1}/([a-z]{1,})_[a-z0-9A-Z]{1,}_[a-z0-9]{1,}_[a-z0-9]{1,}_final.csv$",convert=TRUE,remove=FALSE) %>% 
    tidyr::extract(filepath,"promoter",".*/inference_output_[1-9]{1}/[a-z]{1,}_([a-z0-9A-Z]{1,})_[a-z0-9]{1,}_[a-z0-9]{1,}_final.csv$",convert=TRUE,remove=FALSE)
  return(df)}
}



recover_log_likelihood_old <- function(folderName){
  fileList <- list.files(folderName,pattern='*_final.csv',full.names = TRUE)
  toLikelihoodPath <- function(.l){
    .i <- unlist(str_locate_all(.l,"_final.csv$"))[[1]]
    .new_l <- substr(.l,1,.i-1)
    return(paste(.new_l,".csv",sep=""))}
  likelihoodFileList <- lapply(fileList,toLikelihoodPath)
  if(length(likelihoodFileList)==0){
    return(tibble())
  }else{
    df <- likelihoodFileList %>% 
      lapply(function(.l) readr::read_csv(.l,skip=14) %>%
               tail(1) %>% 
               select(iteration,likelihood) %>% 
               mutate(filepath=.l)) %>%
      do.call(rbind,.) %>% 
      extract(filepath,"replica",".*/inference_output_([1-9]{1})/[a-z]{1,}_[a-z0-9A-Z]{1,}_[a-z0-9]{1,}_[a-z0-9]{1,}.csv$",convert=TRUE,remove=FALSE) %>% 
      extract(filepath,"condition",".*/inference_output_[1-9]{1}/([a-z]{1,})_[a-z0-9A-Z]{1,}_[a-z0-9]{1,}_[a-z0-9]{1,}.csv$",convert=TRUE,remove=FALSE) %>% 
      extract(filepath,"promoter",".*/inference_output_[1-9]{1}/[a-z]{1,}_([a-z0-9A-Z]{1,})_[a-z0-9]{1,}_[a-z0-9]{1,}.csv$",convert=TRUE,remove=FALSE)
    return(df)
      }
}

generate_conditions_df <- function(){
  input_path_rep <- paste(inference_input_path,replica,sep="_")
  output_path_rep <- paste(inference_output_path,replica,sep="_")
  ref_tibble <- tibble()
  k <- 0
  for(i in 1:4){
    for(j in 1:8){
      k <- k+1
      cond_vec <- list(conditions[i],promoters[j],k)
      ref_tibble <- rbind(ref_tibble,cond_vec)
    }}
  colnames(ref_tibble) <- list("condition","promoter","reference")
  return(ref_tibble)}

plot_inferred_final <- function(prom){
  
  avoid_null_plot <- function(.df,cond){
    if(nrow(.df)==0){
      plot <- ggplot()+
        theme_void()
      return(plot)
    }else{
      plot <- .df %>% 
        ggplot()+
        geom_point(aes(replica,init),col="blue")+
        geom_point(aes(replica,final),col="red")+
        facet_wrap(~name,scales="free")+
        theme_cowplot()+
        labs(subtitle=sprintf("%s,%s",prom,cond))+
        coord_cartesian(xlim=c(0,7))+
        ylab("parameter value")}
    return(plot)}
  
  p1 <- avoid_null_plot(inference_results_df %>% 
                          filter(condition=="acetate",promoter==prom),"acetate")
  
  p2 <- avoid_null_plot(inference_results_df %>% 
                          filter(condition=="glycerol",promoter==prom),"glycerol")
  
  p3 <- avoid_null_plot(inference_results_df %>% 
                          filter(condition=="glucose",promoter==prom),"glucose")
  
  p4 <- avoid_null_plot(inference_results_df %>% 
                          filter(condition=="glucoseaa",promoter==prom),"glucoseaa")
  
  ptot <- plot_grid(p1,p2,p3,p4,nrow=2,ncol=2)
  return(ptot)
}

plot_inferred_final_single_condition <- function(prom,cond){
  
  avoid_null_plot <- function(.df,cond){
    if(nrow(.df)==0){
      plot <- ggplot()+
        theme_void()
      return(plot)
    }else{
      
      .dummy <- .df %>% 
        group_by(name) %>% 
        mutate(max_param_final=2*max(final,na.rm=TRUE)) %>% 
        ungroup() %>% 
        distinct(name,.keep_all = TRUE)
      
      plot <- .df %>%
        ungroup() %>% 
        ggplot()+
        geom_point(aes(replica,init),col="blue",alpha=0.5)+
        geom_point(aes(replica,final),col="red",alpha=0.5)+
        geom_blank(data=.dummy,aes(x=1,y=max_param_final))+
        facet_wrap(~name,scales="free")+
        theme_cowplot()+
        labs(subtitle=sprintf("%s,%s",prom,cond))+
        coord_cartesian(xlim=c(0,7))+
        expand_limits(y=0)+
        ylab("parameter value")}
    return(plot)}
  
  p1 <- avoid_null_plot(inference_results_df %>% 
                          filter(condition==cond,promoter==prom),cond)
  return(p1)
}

plot_inferred_final_single_condition_w_error <- function(prom,cond){
  
  avoid_null_plot <- function(.df,cond){
    if(nrow(.df)==0){
      plot <- ggplot()+
        theme_void()
      return(plot)
    }else{
      
      .dummy <- .df %>% 
        group_by(name) %>% 
        mutate(max_param_final=2*max(final,na.rm=TRUE)) %>% 
        ungroup() %>% 
        distinct(name,.keep_all = TRUE)
      
      plot <- .df %>%
        ungroup() %>% 
        ggplot()+
        geom_point(aes(replica,init),col="blue",alpha=0.5)+
        geom_point(aes(replica,final),col="red",alpha=0.5)+
        geom_errorbar(aes(x=replica,ymin=final-sqrt(err),ymax=final+sqrt(err)),col="red",alpha=0.5)+
        geom_blank(data=.dummy,aes(x=1,y=max_param_final))+
        facet_wrap(~name,scales="free")+
        theme_cowplot()+
        labs(subtitle=sprintf("%s,%s",prom,cond))+
        coord_cartesian(xlim=c(0,7))+
        expand_limits(y=0)+
        ylab("parameter value")}
    return(plot)}
  
  p1 <- avoid_null_plot(inference_results_df %>% 
                          filter(condition==cond,promoter==prom),cond)
  return(p1)
}

plot_likelihood_final <- function(prom){
  
  avoid_null_plot <- function(.df,cond){
    if(nrow(.df)==0){
      plot <- ggplot()+
        theme_void()
      return(plot)
    }else{
      plot <- .df %>% 
        ggplot()+
        geom_point(aes(replica,log_likelihood),col="blue")+
        theme_cowplot()+
        labs(subtitle=sprintf("%s,%s",prom,cond))+
        coord_cartesian(xlim=c(0,7))+
        ylab("total log likelihood")}
    return(plot)}
  
  p1 <- avoid_null_plot(inference_results_df %>% 
                          filter(condition=="acetate",promoter==prom),"acetate")
  
  p2 <- avoid_null_plot(inference_results_df %>% 
                          filter(condition=="glycerol",promoter==prom),"glycerol")
  
  p3 <- avoid_null_plot(inference_results_df %>% 
                          filter(condition=="glucose",promoter==prom),"glucose")
  
  p4 <- avoid_null_plot(inference_results_df %>% 
                          filter(condition=="glucoseaa",promoter==prom),"glucoseaa")
  
  ptot <- plot_grid(p1,p2,p3,p4,nrow=2,ncol=2)
  return(ptot)
}

plot_likelihood_final_nogluaa <- function(prom){
  
  avoid_null_plot <- function(.df,cond){
    if(nrow(.df)==0){
      plot <- ggplot()+
        theme_void()
      return(plot)
    }else{
      plot <- .df %>% 
        ggplot()+
        geom_point(aes(replica,log_likelihood),col="blue")+
        theme_cowplot()+
        labs(subtitle=sprintf("%s,%s",prom,cond))+
        coord_cartesian(xlim=c(0,7))+
        ylab("total log likelihood")}
    return(plot)}
  
  p1 <- avoid_null_plot(inference_results_df %>% 
                          filter(condition=="acetate",promoter==prom),"acetate")
  
  p2 <- avoid_null_plot(inference_results_df %>% 
                          filter(condition=="glycerol",promoter==prom),"glycerol")
  
  p3 <- avoid_null_plot(inference_results_df %>% 
                          filter(condition=="glucose",promoter==prom),"glucose")
  
  ptot <- plot_grid(p1,p2,p3,nrow=2,ncol=2)
  return(ptot)
}



  
