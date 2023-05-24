generate_traces <- function(control,rep){
  
  if(control==TRUE){
    alpha_p <- 0
  }else{
    alpha_p <- runif(1,min_alpha_p,max_alpha_p)}
  
  generate_single_trace <- function(r){
    gr <- rnorm(1,mean_x,std_x)*log(2)/60
    lod <- lod_ini + gr*time_min
    residuals <- rnorm(length(lod),0,exp_err_od)
    #residuals <- rnorm(length(lod),0,0)
    od_noisy <- exp(lod)+residuals
    lod_noisy <- log(od_noisy)
    f <- exp(lod)*(alpha_0+alpha_p)+beta
    #f_noise <- rnorm(length(f),0,0)
    #f_noise <- rnorm(length(f),0,sqrt(mean(f)))
    #f_noise <- rnorm(length(f),0,100)
    f_noise <- c(rnorm(1,0,exp_err_f*f[1]))
    for(i in c(2:length(f))){
      f_noise <- c(f_noise,rnorm(1,0,exp_err_f*f[i]))
      #f_noise <- c(f_noise,rnorm(1,0,0.001*f[i]))
    }
    f_noisy <- f+f_noise
    new_df <- tibble(time_min=time_min,corrected_od=od_noisy,fluo=f_noisy,replicate=r,alpha_p=alpha_p,alpha_0=alpha_0,beta=beta)
    return(new_df)}
  
  .final_df <- lapply(c(1:rep), generate_single_trace)
  .final_df <- do.call(rbind,.final_df)
  return(.final_df)
}

compute_dL_dB <- function(.beta,.mydata_inference){
  compute_dL_dB_p <- function(map2,mbp2,mapbp,np,.beta){
    val_p <- np*(mapbp-mbp2*.beta)/(map2+mbp2*.beta**2-2*.beta*mapbp)}
  
  new_df <- .mydata_inference %>% 
    mutate(dL_dB_p=compute_dL_dB_p(mean_ap2,mean_bp2,mean_apbp,Np,.beta))
  dL_dB <- sum(new_df$dL_dB_p)
  return(dL_dB)
}

dichotomic_search <- function(.beta_max_ini,.mydata_inference){
  beta <- 0
  dL_dB <- compute_dL_dB(beta,.mydata_inference)
  if(dL_dB<=0){
    print(sprintf("beta,dL_dB = %s,%s",as.character(beta),as.character(dL_dB)))
    return(0)
  }else if(dL_dB>=0){
    beta_min <- 0
    beta_max <- .beta_max_ini}
  
  dL_dB <- compute_dL_dB(beta_max,.mydata_inference)
  
  while(dL_dB>0){
    beta_max <- beta_max*2
    dL_dB <- compute_dL_dB(beta_max,.mydata_inference)
  }
  print("Negative dL_dB, beginning dichotomic search")
  print(sprintf("With beta_max,dL_dB = %s,%s",as.character(beta_max),as.character(dL_dB)))
  
  while((2*abs(beta_min-beta_max)/(beta_max+beta_min))>1e-3){
    print(sprintf("beta_max,beta_min = %s,%s",as.character(beta_max),as.character(beta_min)))
    beta <- (beta_max+beta_min)/2
    dL_dB <- compute_dL_dB(beta,.mydata_inference)
    print(sprintf("beta,dL_dB = %s,%s",as.character(beta),as.character(dL_dB)))
    if(dL_dB>=0){
      beta_min <- beta
    }else{
      beta_max <- beta}}
  
  return(beta)
}


compute_wp <- function(.Bp,.Qp2,np,.beta){
  wp_val <- (np-1)/((.beta-.Bp)**2 + .Qp2)
  return(wp_val)}

compute_beta <- function(.df){
  .new_df <- .df %>% 
    mutate(num=wp*Bp)
  val <- (sum(.new_df$num))/(sum(.new_df$wp))
  return(val)}

compute_error_beta <- function(.beta,.mydata_inference){
  .new_df <- .mydata_inference %>% 
    mutate(dL2_dB2=(Np-1)*(Qp2-(.beta-Bp)**2)/(Qp2+(.beta-Bp)**2)**2)
  val <- 1/sum(.new_df$dL2_dB2)
  return(val)}

iterative_search <- function(.beta_ini,.mydata_inference){
  
  .old_beta <- .beta_ini
  
  .df <- .mydata_inference %>% 
    mutate(wp=compute_wp(Bp,Qp2,Np,.old_beta))
  
  .new_beta <- compute_beta(.df)
  
  while(abs(.new_beta-.old_beta)>0.01){
    print(sprintf("old_beta,new_beta=%s,%s",as.character(.old_beta),as.character(.new_beta)))
    .old_beta <- .new_beta
    .df <- .mydata_inference %>% 
      mutate(wp=compute_wp(Bp,Qp2,Np,.old_beta))
    .new_beta <- compute_beta(.df)}
  
  return(.new_beta)
  
}


compute_error <- function(a,varx,vary,covarxy){
  v <- varx*vary-covarxy**2
  dlP2_da2 <- 4*(1-a**2)/(a**2+1)**2+2*varx**2*((varx*a-covarxy)**2-v)/(((varx*a-covarxy)**2+v)**2)+2*vary**2*((vary*a+covarxy)**2-v)/(((vary*a+covarxy)**2+v)**2)
  sigma <- sqrt(1/(abs(dlP2_da2)))
}


