run_inference_clean <- function(cond,reps,.mydata_curated){
  
    mydata_curated_cond <- .mydata_curated %>% 
      filter(!grepl("empty",strain)) %>% #filter out empty wells
      filter(condition==cond) %>% 
      filter(replica %in% reps) %>% 
      #filter(valid=="y") %>% 
      #filter(between(corrected_od,min_corrected_od,max_corrected_od)) %>% 
      filter(max(od)>min_corrected_od) %>% #to get rid of wells in which there are no growth at all
      ungroup() %>% 
      group_by(description,plate,well) %>%
      arrange(time_min) %>% 
      mutate(above_max_time=ifelse(corrected_od>=max_corrected_od,time_min,NA)) %>% 
      mutate(time_limit=min(c(above_max_time,Inf),na.rm=TRUE)) %>% 
      filter(time_min<time_limit) %>% 
      mutate(below_min_time=ifelse(corrected_od<=min_corrected_od,time_min,NA)) %>% 
      mutate(time_limit=max(c(below_min_time,-Inf),na.rm=TRUE)) %>% 
      filter(time_min>time_limit) %>% 
      ungroup() %>% 
      select(-c(above_max_time,time_limit,below_min_time))

    mydata_inference <- mydata_curated_cond %>% 
      group_by(description,well) %>% 
      mutate(mean_xy_p=mean(corrected_od*fluo)) %>% 
      mutate(mean_x2_p=mean((corrected_od)**2)) %>% 
      mutate(mean_y2_p=mean((fluo)**2)) %>% 
      mutate(mean_x_p=mean(corrected_od)) %>%
      mutate(mean_y_p=mean(fluo)) %>%
      mutate(var_x_p=var(corrected_od)) %>% 
      mutate(Bp=(mean_y_p*mean_x2_p-mean_xy_p*mean_x_p)/(var_x_p)) %>% 
      mutate(Qp2=(mean_y2_p*mean_x2_p-mean_xy_p**2)/(var_x_p)-Bp**2) %>% 
      mutate(Np=n()) %>% 
      ungroup() %>% 
      mutate(aip=fluo-corrected_od*mean_xy_p/mean_x2_p) %>% 
      mutate(bip=1-corrected_od*mean_x_p/mean_x2_p) %>% 
      group_by(description,well) %>% 
      mutate(mean_apbp=mean(aip*bip),
             mean_ap2=mean(aip**2),
             mean_bp2=mean(bip**2)) %>% 
      mutate(Np=n()) %>% 
      ungroup() %>% 
      distinct(description,well,.keep_all = TRUE)

    #approx_beta <- mydata_curated %>%
    #  group_by(strain) %>%
    #  mutate(intercept=linear_mod_intercept(fluo,corrected_od)) %>% 
    #  ungroup() %>% 
    #  distinct(promoter,replicate,.keep_all=TRUE) %>% 
    #  .$intercept %>% 
    #  mean
    
    approx_beta <- 500
    
    opt_beta_iterative <- iterative_search(approx_beta,mydata_inference)
    opt_beta_dichotomic <- dichotomic_search(approx_beta,mydata_inference)
    err_beta_iterative <- compute_error_beta(opt_beta_iterative,mydata_inference)
    err_beta_dichotomic <- compute_error_beta(opt_beta_dichotomic,mydata_inference)

    #Setting some error bars on beta
    
    mydata_inference <- mydata_inference %>% 
      mutate(beta_predict_dic=opt_beta_dichotomic) %>% 
      mutate(beta_predict_ite=opt_beta_iterative) %>% 
      mutate(beta_var_dic=err_beta_dichotomic) %>%
      mutate(beta_var_ite=err_beta_iterative) %>% 
      mutate(beta_sd_ite=sqrt(beta_var_ite)) %>% 
      mutate(beta_sd_dic=sqrt(beta_var_dic)) %>%
      mutate(alpha_tot_predict=(mean_xy_p-beta_predict_ite*mean_x_p)/(mean_x2_p)) %>% 
      mutate(var_alpha_tot=var_x_p/(mean_x2_p**2)*((beta_predict_ite-Bp)**2+Qp2)/Np) %>% 
      mutate(sd_alpha_tot=(sqrt(var_alpha_tot)))
    
    #mydata_inference %>% 
    #  select(beta_predict_ite,beta_sd_ite,beta_predict_dic,beta_sd_dic) %>% 
    #  filter(row_number()==1)
    
    mydata_inference <- mydata_inference %>% 
      select(description,well,alpha_tot_predict,sd_alpha_tot,beta_predict_ite,beta_sd_ite,replica,condition)
    
    return(mydata_inference)
}
