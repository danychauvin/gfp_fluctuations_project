# 20201211, Dany Chauvin
# INSTALLING AND LOADING NECESSARY PACKAGES 
# Following packages are necessary to import data from deepMoma and analyze these data.

# Uncomment below when using for the first time.
#options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/__linux__/centos7/latest"))
#install.packages("tidyverse")
#install.packages("cowplot")
#install.packages("devtools")
#install.packages("ggcorrplot")
#install.packages("ggpubr")
#install.packages("comprehenr")
#devtools::install_github(c('hadley/multidplyr'))
#devtools::install_github(c('julou/ggCustomTJ'))
#devtools::install_github('vanNimwegenLab/vngMoM',auth_token="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")


library(tidyverse)
library(cowplot)
library(devtools)
library(multidplyr)
library(ggCustomTJ)
library(vngMoM)
library(ggcorrplot)
library(ggpubr)
library(nloptr)
library(RColorBrewer)
library(RcppArmadillo)
library(comprehenr)
# Notes about Rcpp Armadillo. It is already installed together with Lapack in
# /scicore/soft/apps/R/4.1.0-foss-2018b/lib64/R
# But for whathever reason, vngMoM installs it too, on top of the already incorporated install.
# So use: remove.packages("RcppArmadillo").
# WARNING: libstdc++.so.6 library, with proper GLIBCXX versions is necessary to use RcppArmadillo.
# Be sure that GCCcore/8.3.0 is loaded/installed.
# To do so while using Rstudio on the scicore server "service06", do the following, add the following to your ~/.bashrc file:
# `if [[ "$HOSTNAME" = *service06* ]]; then
#     ml GCCcore/8.3.0
#  fi'

# Then use renv::
#install.packages("renv")
#renv::init()
#renv::snapshot()

# ------------------------------------ Functions necessary for the import and transformation of MM Data -----------------------------------------------------------------

# Important functions that are used to import the data
# SET NECESSARY FUNCTIONS TO GENERATE PATHS TO PREPROCESSED DATA FILES ####
data2preproc_file <- function(.f)
    basename(.f) %>% sub("ExportedCellStats_", "", .) %>% 
    file_path_sans_ext %>% paste0("_frames.txt")
data2preproc <- function(.f)
    file.path(data2preproc_dir(.f), data2preproc_file(.f))

# EXPONENTIAL FIT FOR VOLUME AND LENGTH
# Here only considering error in y axis

fit_exp_elongation_predict <- function(.x_l,.y_l){
  .x <- .x_l
  .y <- log(.y_l)
  
  mxy <- mean(.x*.y,na.rm=TRUE) 
  mx2 <- mean(.x**2,na.rm=TRUE)
  my2 <- mean(.y**2,na.rm=TRUE)
  mx <-  mean(.x,na.rm=TRUE)
  my <- mean(.y,na.rm=TRUE)

  .Vxx <- mx2-mx**2
  .Vyy <- my2-my**2
  .Vxy <- mxy-mx*my
  
  .s <- .Vxy/.Vxx
  .i <- mean(.y,na.rm=TRUE)-.Vxy/.Vxx*mean(.x,na.rm=TRUE)
  .output <- exp(.s*.x+.i)
  return(.output)}

fit_exp_elongation_slope <- function(.x_l,.y_l){
  .x <- .x_l
  .y <- log(.y_l)
  
  mxy <- mean(.x*.y,na.rm=TRUE) 
  mx2 <- mean(.x**2,na.rm=TRUE)
  my2 <- mean(.y**2,na.rm=TRUE)
  mx <-  mean(.x,na.rm=TRUE)
  my <- mean(.y,na.rm=TRUE)
  
  .Vxx <- mx2-mx**2
  .Vyy <- my2-my**2
  .Vxy <- mxy-mx*my
  
  .s <- .Vxy/.Vxx
  .i <- mean(.y,na.rm=TRUE)-.Vxy/.Vxx*mean(.x,na.rm=TRUE)
  
  .output <- rep(c(.s),times=length(.x_l))
  return(.output)}

fit_exp_elongation_sd_slope <- function(.x_l,.y_l){
  .x <- .x_l
  .y <- log(.y_l)
  
  mxy <- mean(.x*.y,na.rm=TRUE) 
  mx2 <- mean(.x**2,na.rm=TRUE)
  my2 <- mean(.y**2,na.rm=TRUE)
  mx <-  mean(.x,na.rm=TRUE)
  my <- mean(.y,na.rm=TRUE)
  
  .Vxx <- mx2-mx**2
  .Vyy <- my2-my**2
  .Vxy <- mxy-mx*my
  
  .s <- .Vxy/.Vxx
  .i <- mean(.y,na.rm=TRUE)-.Vxy/.Vxx*mean(.x,na.rm=TRUE)
  
  .p <- length(.y)
  
  .y_pred <- .x*.s+.i
  .sse <- sum((.y-.y_pred)**2)/(.p-2)
  .sd <- sqrt(.sse/sum((.x-mx)**2))
    
  .output <- .sd
  .output <- rep(c(.output),times=length(.x_l))
  return(.output)}
  
fit_exp_elongation_intercept <- function(.x_l,.y_l){
  .x <- .x_l
  .y <- log(.y_l)
  
  mxy <- mean(.x*.y,na.rm=TRUE) 
  mx2 <- mean(.x**2,na.rm=TRUE)
  my2 <- mean(.y**2,na.rm=TRUE)
  mx <-  mean(.x,na.rm=TRUE)
  my <- mean(.y,na.rm=TRUE)
  
  .Vxx <- mx2-mx**2
  .Vyy <- my2-my**2
  .Vxy <- mxy-mx*my
  
  .s <- .Vxy/.Vxx
  .i <- mean(.y,na.rm=TRUE)-.Vxy/.Vxx*mean(.x,na.rm=TRUE)
  
  .p <- length(.y)
  
  .y_pred <- .x*.s+.i
  .sse <- sum((.y-.y_pred)**2)/(.p-2)
  .sd <- sqrt(.sse/sum((.x-mx)**2))
  
  .output <- .i
  .output <- rep(c(.output),times=length(.x_l))
  return(.output)}

# COMPUTE CELL VOLUME

compute_vol <- function(.l,.w){
    .vol_um <- ((.l-.w)*(.w/2)**2*3.14)+4/3*3.14*(.w/2)**3
    return(.vol_um)}

#-------------------------------------------------- Function necessary for the import and transformation of Bulk Data --------------------------------------------------

read_Biotek_Synergy2_kinetic <- function(.path) {
  #Extract the folder name where the spec data are stored: this is the expId. Must be of the sort: date_description. Ex: 20191115_test
  #  read_spec_kinetic() is written such as to call read.table only once 
  # (much faster than calling it for each timepoint and channel)
  #extract_date_time
  .measurement_date <- str_match(.path,"[a-zA-Z0-9]{1,}_[a-zA-Z0-9]{1,}_([0-9]{8})_[0-9]{6}_FE_.txt$")[[2]]
  .measurement_time <- str_match(.path,"[a-zA-Z0-9]{1,}_[a-zA-Z0-9]{1,}_[0-9]{8}_([0-9]{6})_FE_.txt$")[[2]]
  #print(measurement_date)
  #print(measurement_time)
  
  .lines <- readLines(.path)
  #print(.lines)
  .od_ch <- .lines[[1]]
  .fluo_ch <- .lines[[12]]
  .od_values <- unlist(.lines[c(3:10)])
  .fluo_values <- unlist(.lines[c(14:21)])
  
  noFirstnoLast <- function(.l){
    new_l <- str_split(.l,",")[[1]]
    new_l <- new_l[2:(length(new_l)-1)]
    return(new_l)}
  
  noFirstnoSecondLast <- function(.l){
    new_l <- str_split(.l,",")[[1]]
    new_l <- new_l[2:(length(new_l)-2)]
    return(new_l)}
  
  .formatted_od_values <- lapply(.od_values,noFirstnoLast) %>% unlist()
  .formatted_fluo_values <- lapply(.fluo_values,noFirstnoSecondLast) %>% unlist()
  .rows <- rep(LETTERS[c(1:8)],each=12)
  .cols <- rep(c(1:12),times=8)
  new_df <- data_frame(row=.rows,column=.cols,od=.formatted_od_values,fluo=.formatted_fluo_values,measurement_date=.measurement_date,measurement_time=.measurement_time,od_channel=.od_ch,fluo_channel=.fluo_ch)
}

read_plate_layout <- function(.path) {
  #Extract the folder name where the spec data are stored: this is the expId. Must be of the sort: date_description. Ex: 20191115_test
  #  read_spec_kinetic() is written such as to call read.table only once 
  # (much faster than calling it for each timepoint and channel)
  .lines <- readLines(.path)
  #print(.lines)
  .l_idx <- stringr::str_detect(.lines, "Type:") %>% which %>% (function(.x) .x)
  
  noFirst <- function(.l){
    new_l <- str_split(.l,",")[[1]]
    new_l <- new_l[2:length(new_l)]
    return(new_l)}
  
  return_plate <- function(.index){
    .col_title <- str_match(.lines[[.index]],"Type: ([a-zA-Z0-9]{1,}),,,,,,,,,,,,$")[[2]]
    .data <- .lines[c((.index+2):(.index+9))]
    .values <- lapply(.data,noFirst) %>% unlist()
    .rows <- rep(LETTERS[c(1:8)],each=12)
    .cols <- rep(c(1:12),times=8)
    new_df <- data_frame(row=.rows,column=.cols,values=.values,type=rep(c(.col_title),each=96))
    return(new_df)}
  
  new_df <- .l_idx %>% lapply(return_plate) %>% 
    bind_rows() %>% 
    pivot_wider(id_cols=c(row,column),names_from=type,values_from=values)
  
  return(new_df)}


# Here a few modeling functions that have written in order to simplify predictions efforts with dplyr.

#predict_od
#input: od, time_min
#produces a linear fit between log(od) and time_min
#output: predicted_od = OD_0 * exp(alpha * time_min) where OD_0 and alpha are inferred

predict_od <- function(.c_od,.t_min){
  if (length(.c_od) != length(.t_min)) 
    stop("variables of different lengths.")
  if (length(.c_od) < 2) 
    return(NA)
  stats::lm(log(.c_od) ~ .t_min) %>% stats::predict(se.fit = TRUE) %>% 
    .[["fit"]] %>% exp
}

#predict_alpha
#input: od, time_min
#produces a linear fit between log(od) and time_min
#output: alpha, the exponential growth rate, in min-1

predict_alpha <- function(.p_od,.t_min){
  if (length(.p_od) != length(.t_min)) 
    stop("variables of different lengths.")
  if (length(.p_od) < 2) 
    return(NA)
  alpha <- (log(last(.p_od))-log(first(.p_od)))/(last(.t_min)-first(.t_min))
  return(alpha)
}

predict_fluo_od <- function(.c_f,.c_od){
  if (length(.c_f) != length(.c_od)) 
    stop("variables of different lengths.")
  if (length(.c_f) <= 2) 
    return(NA)
  stats::lm(log(.c_f) ~ log(.c_od)) %>% stats::predict(se.fit = TRUE) %>% 
    .[["fit"]] %>% exp
}


predict_lfluo_od <- function(.c_f,.c_od){
  if (length(.c_f) != length(.c_od)) 
    stop("variables of different lengths.")
  if (length(.c_f) <= 2) 
    return(NA)
  stats::lm(.c_f ~ .c_od) %>% stats::predict(se.fit = TRUE) %>% 
    .[["fit"]] %>% exp
}

#########################################################################################



#predict_intercept_fluo_od
#input: corrected_fluo,corrected_od
#produces a linear fit between log(fluo) and log(od) (power law of the type fluo = A * od ^ B)
#output: coefficient A

predict_intercept_fluo_od <- function(.c_f,.c_od){
  if (length(.c_f) != length(.c_od)) 
    stop("variables of different lengths.")
  if (length(.c_f) <= 2) 
    return(NA)
  stats::lm(log(.c_f) ~ log(.c_od)) %>% stats::coef(se.fit = TRUE) %>% 
    .[[1]] %>% exp %>% rep(times=length(.c_f))
}

#predict_power_fluo_od
#input: corrected_fluo,corrected_od
#produces a linear fit between log(fluo) and log(od) (power law of the type fluo = A * od ^ B)
#output: coefficient B

predict_power_fluo_od <- function(.c_f,.c_od){
  if (length(.c_f) != length(.c_od)) 
    stop("variables of different lengths.")
  if (length(.c_f) <= 2) 
    return(NA)
  stats::lm(log(.c_f) ~ log(.c_od)) %>% stats::coef(se.fit = TRUE) %>% 
    .[[2]] %>% rep(times=length(.c_f))
}


#predict_intercept_fluo_od
#input: corrected_fluo,corrected_od
#produces a linear fit between log(fluo) and log(od) (power law of the type fluo = A * od ^ B)
#output: coefficient A

predict_lintercept_fluo_od <- function(.c_f,.c_od){
  if (length(.c_f) != length(.c_od)) 
    stop("variables of different lengths.")
  if (length(.c_f) <= 2) 
    return(NA)
  stats::lm(.c_f ~ .c_od) %>% stats::coef(se.fit = TRUE) %>% 
    .[[1]] %>% rep(times=length(.c_f))
}

#predict_power_fluo_od
#input: corrected_fluo,corrected_od
#produces a linear fit between log(fluo) and log(od) (power law of the type fluo = A * od ^ B)
#output: coefficient B

predict_lslope_fluo_od <- function(.c_f,.c_od){
  if (length(.c_f) != length(.c_od)) 
    stop("variables of different lengths.")
  if (length(.c_f) <= 2) 
    return(NA)
  stats::lm(.c_f ~ .c_od) %>% stats::coef(se.fit = TRUE) %>% 
    .[[2]] %>% rep(times=length(.c_f))
}

#Example of application
#mydata %>% 
#  group_by(condition) %>% 
#  mutate(predicted_fluo=predict_fluo_od(corrected_fluo,corrected_od)) %>% #Predict fluo
#  mutate(predicted_intercept_fluo_od=predict_intercept_fluo_od(corrected_fluo,corrected_od)) %>% #Predict intercept
#  mutate(predicted_power_fluo_od=predict_power_fluo_od(corrected_fluo,corrected_od)) %>% #Predict power
#  ungroup() %>% 
#  ggplot()+
#  geom_point(aes(corrected_od,corrected_fluo,col=interaction(plate,row,column)),alpha=0.7)+ # Compare predicted data and experimental data
#  geom_line(aes(corrected_od,predicted_fluo,col=interaction(plate,row,column)),alpha=0.7,col="red")+
#  geom_line(aes(corrected_od,predicted_intercept_fluo_od*corrected_od**predicted_power_fluo_od,col=interaction(plate,row,column)),alpha=0.7,col="green")+
#  facet_wrap(~condition,scales="free")+
#  scale_color_manual(values = getPalette(colourCount))+
#  theme_cowplot()

# To assess how good a fit is, we use the Pearson determination coefficient, between prediction and experimental values
compute_r2 <- function(.e,.p){
  if (length(.e) != length(.p)) 
    stop("variables of different lengths.")
  if (length(.e) < 2) 
    return(NA)
  r2 <- 1-sum((.e-.p)**2)/sum((.e-mean(.e))**2) 
  return(r2)
}


linear_mod_predict <- function(.c_f,.c_od){
  if (length(.c_f) != length(.c_od)) 
    stop("variables of different lengths.")
  if (length(.c_f) <= 2) 
    return(NA)
  stats::lm(.c_f ~ .c_od) %>% stats::predict(se.fit = TRUE) %>% 
    .[["fit"]]
}

linear_mod_intercept <- function(.c_f,.c_od){
  if (length(.c_f) != length(.c_od)) 
    stop("variables of different lengths.")
  if (length(.c_f) <= 2) 
    return(NA)
  stats::lm(.c_f ~ .c_od) %>% stats::coef(se.fit = TRUE) %>% 
    .[[1]] %>% rep(times=length(.c_f))
}

linear_mod_slope <- function(.c_f,.c_od){
  if (length(.c_f) != length(.c_od)) 
    stop("variables of different lengths.")
  if (length(.c_f) <= 2) 
    return(NA)
  stats::lm(.c_f ~ .c_od) %>% stats::coef(se.fit = TRUE) %>% 
    .[[2]] %>% rep(times=length(.c_f))
}

compute_slope_se <- function(.x,.y,.p){
  num <- 1/(length(.x)-2)*sum((.p-.y)**2)
  den <- sum((.x-mean(.x))**2)
  se <- sqrt(num/den)
  return(se)
}

predictdf.lm_right <- 
  function(model, xseq, se, level){
    ## here the main code: truncate to x values at the right
    init_range = range(model$model$x)
    xseq <- xseq[xseq >=init_range[1]]
    ggplot2:::predictdf.default(model, xseq[-length(xseq)], se, level)
  }

predictdf.lm_left <- 
  function(model, xseq, se, level){
    init_range = range(model$model$x)
    ## here the main code: truncate to x values at the left
    xseq <- xseq[xseq <=init_range[2]]
    ggplot2:::predictdf.default(model, xseq[-length(xseq)], se, level)
  }

lm_right <- function(formula,data,...){
  mod <- lm(formula,data)
  class(mod) <- c('lm_right',class(mod))
  mod
}

## decorate lm object with a new class lm_left
lm_left <- function(formula,data,...){
  mod <- lm(formula,data)
  class(mod) <- c('lm_left',class(mod))
  mod
}

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


#------------------------------------------------------------------- Statistical functions-----------------------------------------------------------------------
compute_slope_errxy_robust <- function(.x,.y){
  p <- length(.x)
  .Vxx <- 1/p*sum((.x-1/p*sum(.x))**2)
  .Vyy <- 1/p*sum((.y-1/p*sum(.y))**2)
  .Vxy <- 1/p*sum(.x*.y-1/(p**2)*sum(.x)*sum(.y))
  
  .slope_plus <- (.Vyy-.Vxx)/(2*.Vxy)+sqrt(1+((.Vyy-.Vxx)/(2*.Vxy))**2)#assuming that this is the best slope
  .slope_minus <- (.Vyy-.Vxx)/(2*.Vxy)-sqrt(1+((.Vyy-.Vxx)/(2*.Vxy))**2)#assuming that this is the best slope
  
  .beta_plus <- 1/p*sum(.y-.slope_plus*.x)
  .beta_minus <- 1/p*sum(.y-.slope_minus*.x)
  
  delta <- compute_rmsd(.slope_plus,.beta_plus,.x,.y)-compute_rmsd(.slope_minus,.beta_minus,.x,.y)
  
  if(delta>0){
    .slope <- .slope_minus 
    .beta <- .beta_minus 
    #.sd_slope <- compute_error(.slope,.Vxx,.Vyy,.Vxy)
  }else if(delta<=0){
    .slope <- .slope_plus 
    .beta <- .beta_plus
    #.sd_slope <- compute_error(.slope,.Vxx,.Vyy,.Vxy)
    }
  else{
    .slope <- NA
    .beta <- NA}
  
  return(rep(c(.slope),times=length(.x)))
}

compute_intercept_errxy_robust <- function(.x,.y){
  p <- length(.x)
  
  .Vxx <- 1/p*sum((.x-1/p*sum(.x))**2)
  .Vyy <- 1/p*sum((.y-1/p*sum(.y))**2)
  .Vxy <- 1/p*sum(.x*.y-1/(p**2)*sum(.x)*sum(.y))
  
  .slope_plus <- (.Vyy-.Vxx)/(2*.Vxy)+sqrt(1+((.Vyy-.Vxx)/(2*.Vxy))**2)#assuming that this is the best slope
  .slope_minus <- (.Vyy-.Vxx)/(2*.Vxy)-sqrt(1+((.Vyy-.Vxx)/(2*.Vxy))**2)#assuming that this is the best slope
  
  .beta_plus <- 1/p*sum(.y-.slope_plus*.x)
  .beta_minus <- 1/p*sum(.y-.slope_minus*.x)
  
  delta <- compute_rmsd(.slope_plus,.beta_plus,.x,.y)-compute_rmsd(.slope_minus,.beta_minus,.x,.y)
  
  if(delta>0){
    .slope <- .slope_minus 
    .beta <- .beta_minus 
    #.sd_slope <- compute_error(.slope,.Vxx,.Vyy,.Vxy)
  }else if(delta<=0){
    .slope <- .slope_plus 
    .beta <- .beta_plus
    #.sd_slope <- compute_error(.slope,.Vxx,.Vyy,.Vxy)
    }
  else{
    .slope <- NA
    .beta <- NA}
  
  return(rep(c(.beta),times=length(.x)))
}

compute_sdslope_errxy_robust <- function(.x,.y){
  p <- length(.x)
  
  .Vxx <- 1/p*sum((.x-1/p*sum(.x))**2)
  .Vyy <- 1/p*sum((.y-1/p*sum(.y))**2)
  .Vxy <- 1/p*sum(.x*.y-1/(p**2)*sum(.x)*sum(.y))
  
  .slope_plus <- (.Vyy-.Vxx)/(2*.Vxy)+sqrt(1+((.Vyy-.Vxx)/(2*.Vxy))**2)#assuming that this is the best slope
  .slope_minus <- (.Vyy-.Vxx)/(2*.Vxy)-sqrt(1+((.Vyy-.Vxx)/(2*.Vxy))**2)#assuming that this is the best slope
  
  .beta_plus <- 1/p*sum(.y-.slope_plus*.x)
  .beta_minus <- 1/p*sum(.y-.slope_minus*.x)
  
  delta <- compute_rmsd(.slope_plus,.beta_plus,.x,.y)-compute_rmsd(.slope_minus,.beta_minus,.x,.y)
  
  if(delta>0){
    .slope <- .slope_minus 
    .beta <- .beta_minus 
    .sd_slope <- compute_error(.slope,.Vxx,.Vyy,.Vxy,p)
  }else if(delta<=0){
    .slope <- .slope_plus 
    .beta <- .beta_plus
    .sd_slope <- compute_error(.slope,.Vxx,.Vyy,.Vxy,p)}
  else{
    .slope <- NA
    .beta <- NA}
  
  return(rep(c(.sd_slope),times=length(.x)))
}

compute_error <- function(a,varx,vary,covarxy,m){
  # With Erik's corrections
  #Computes the error on the best slope a (computed with compute_slope_errxy_robust), given the variances of x and y, s well as the covariance.
  v <- varx*vary-covarxy**2
  ax <- covarxy/varx
  ay <- covarxy/vary
  dlP2_da2 <- (m-1)*(2*(1-a**2)/(a**2+1)**2+
                       ((a-ax)**2-v/(varx**2))/(((a-ax)**2+v/(varx**2))**2)+
                       ((a+ay)**2-v/(vary**2))/(((a+ay)**2+v/(vary**2))**2))
  sigma <- sqrt(1/(abs(dlP2_da2)))
}

compute_error_2 <- function(a,varx,vary,covarxy){
  #Computes the error on the best slope a (computed with compute_slope_errxy_robust), given the variances of x and y, s well as the covariance.
  v <- varx*vary-covarxy**2
  dlP2_da2 <- 4*(1-a**2)/(a**2+1)**2+2*varx**2*((varx*a-covarxy)**2-v)/(((varx*a-covarxy)**2+v)**2)+2*vary**2*((vary*a+covarxy)**2-v)/(((vary*a+covarxy)**2+v)**2)
  sigma <- sqrt(1/(abs(dlP2_da2)))
}


compute_rmsd <- function(.s,.b,.x,.y){
  #Computes the root mean square deviation, between a fit characterized by slope .s, intercept .b, and some data (.x,.y)
  .y_pred <- .s*.x+.b
  return(sqrt(sum((.y-.y_pred)**2,na.rm=TRUE)))
}


#----------------------------------------------------------------Functions to manipulate cell phylogeny------------------------------------------------------------

concatenate_traces <- function(df){
  # Input: should have a growth_rate field, depending on the one that you want to consider for the computation.
  # Output: concatenated cell cycles
  cell_and_parent_1 <- df %>% 
    distinct(cell,.keep_all=TRUE) %>% 
    select(cell,parent,medium) %>% 
    rename(degree_1=parent)
  
  cell_and_parent_2 <- df %>%
    distinct(cell,.keep_all = TRUE) %>% 
    select(cell,parent,medium) %>% 
    rename(degree_1=cell) %>% 
    rename(degree_2=parent)
  
  cell_and_parent_3 <- df %>%
    distinct(cell,.keep_all = TRUE) %>% 
    select(cell,parent,medium) %>% 
    rename(degree_2=cell) %>% 
    rename(degree_3=parent)
  
  cell_and_parent_4 <- df %>%
    distinct(cell,.keep_all = TRUE) %>% 
    select(cell,parent,medium) %>% 
    rename(degree_3=cell) %>% 
    rename(degree_4=parent)
  
  phylogenies <- left_join(cell_and_parent_1,cell_and_parent_2,by=c("degree_1","medium")) %>% 
    left_join(cell_and_parent_3,by=c("degree_2","medium")) %>% 
    left_join(cell_and_parent_4,by=c("degree_3","medium")) %>% 
    mutate(degree_0=cell) %>%
    rename(trace_ref=cell) %>% 
    select(trace_ref,degree_0,degree_1,degree_2,degree_3,degree_4,medium) %>% #here limited to 3 consecutive growth
    gather(degree,cell_id,c(degree_0,degree_1,degree_2,degree_3,degree_4)) %>%  
    mutate(degree=as.double(substr(degree,nchar(degree),nchar(degree)))) %>% 
    arrange(trace_ref,desc(degree)) %>% 
    rename(cell=cell_id)
  
  # Now join phylogenies to another df
  concatenatedTraces <- phylogenies %>% 
    left_join(df,by=c("cell","medium")) %>% 
    drop_na() %>% #Dropping cells that do not exist in the dataset
    select(trace_ref,degree,cell,time_sec,medium,growth_rate) %>% 
    arrange(trace_ref,time_sec)
  
  return(concatenatedTraces)
}


#--------------------------------------------------------------Functions to compute autocorrelation functions-------------------------------------------------------------------

autocorrelationFunction <- function(full_df,mean_div_time,dt,nCellCycle,cond){
  results <- c()
  M <- round(mean_div_time/dt*nCellCycle,0) #Autocorrelation is computed over nCellCycle * mean number of data point per cell cycle
  for (i in 1:M){
    results[i] <- autocorrelationFunction_lag(full_df,i)}
  results_df <- tibble(condition=cond,lag_step=c(0:(M-1)),lag_min=c(0:(M-1))*dt,autocorrelation=results)
  return(results_df)
}

autocorrelationFunction_lag <- function(df,lag){
  #input: df with trace_ref,cell,time_sec and growth_rate, lag as an integer
  #output: autocorrelation value for a given lag.
  
  new_df <- df %>% 
    mutate(uid=paste(cell,time_sec,sep=".")) %>% 
    group_by(trace_ref) %>% 
    arrange(time_sec) %>%
    select(trace_ref,uid,growth_rate) %>% 
    do((function(.df){
      N <- nrow(.df)
      if (lag>N){
        return(tibble())
      }else{
        concat <- cbind(
          .df[c(lag:N),],
          .df[c(1:(N-lag+1)),])
        colnames(concat) <- list('trace_ref','uid_1','growth_rate_1','trace_ref_2','uid_2','growth_rate_2')
        return(concat)}
    })(.)) %>% 
    ungroup() %>%
    select(-c(trace_ref,trace_ref_2)) %>% 
    unique()
  corvalue <- (mean(new_df$growth_rate_1*new_df$growth_rate_2)-mean(new_df$growth_rate_1)*mean(new_df$growth_rate_2))/(sd(new_df$growth_rate_1)*sd(new_df$growth_rate_2))
  return(corvalue)}

#-------------------------------------------------------------- Inference functions --------------------------------------------------------------------------------

# FUNCTION TO SET OPTIMIZATION READY

set_optimization_ready <- function(data,cPath,fPath,fName,me,va,ste,num,N_cell){
  
  #Make the necessary folders
  dir.create(fPath,showWarnings=FALSE)
  dir.create(paste(fPath,"/","input",sep=""),showWarnings=FALSE)
  dir.create(paste(fPath,"/","output",sep=""),showWarnings=FALSE)
  #Copy paste the code: the code is always overwritten.
  system(paste("\\cp -r",cPath,fPath,sep=' '))
  
  # First, I will assume that in a given conditions, all growth lanes, all strains are alike.
  # Therefore, I will randomly sample 100 cell cycles.
  cell_sample_for_optimization <- data %>% 
    filter(medium==me) %>%
    group_by(cell) %>% 
    filter(n()>=num) %>% 
    ungroup() %>% 
    distinct(cell) %>% 
    sample_n(N_cell)
  
  sample_for_optimization <- data %>% 
    semi_join(cell_sample_for_optimization,by=c("cell"))
  
  #Adapting the format to the script
  sample_for_optimization <- sample_for_optimization %>%
    rename(cell_ori=cell) %>% 
    #parent_ID=parent_id) %>%
    mutate(pos=paste(pos,orientation,sep=".")) %>%  #NECESSARY CAUSE ATHOS CODE RECOMPUTES WRONGLY THE LANE_ID...
    mutate(pos=paste(pos,"_",sep="")) %>% #NECESSARY CAUSE ATHOS CODE RECOMPUTES WRONGLY THE LANE_ID...
    mutate(gl=paste(gl,"_",sep="")) %>% #NECESSARY CAUSE ATHOS CODE RECOMPUTES WRONGLY THE LANE_ID...
    mutate(cell=paste(lane_ID,id,sep="_"))
  #select(cell,date,time_sec,lane_ID,parent_ID,length_um)
  
  #Write sample
  readr::write_csv(sample_for_optimization,sprintf(paste(fPath,"/input/",fName,sep="")))
  
  #Append jobs to .input/commands.cmd
  lineToWrite <- sprintf("python %s/optimization_code/parameters_find.py %s/input/%s %s %s %s",fPath,fPath,fName,va,ste,num)
  system(sprintf('touch %s/input/commands.cmd',fPath))
  write(lineToWrite,file=sprintf('%s/input/commands.cmd',fPath),append=TRUE,N_cell)
  
  #Modify jobs.sh
  command_path <- sprintf("%s/input/commands.cmd",fPath) 
  results_path <- sprintf("%s/output/results_%%a",fPath)
  errors_path <- sprintf("%s/output/errors_%%a",fPath)
  jobs_path <- sprintf("%s/optimization_code/jobs.sh",fPath)
  system(sprintf("sed -i \"s+path_to_output+%s+g\" %s",results_path,jobs_path))
  system(sprintf("sed -i \"s+path_to_error+%s+g\" %s",errors_path,jobs_path))
  system(sprintf("sed -i \"s+path_to_commands+%s+g\" %s",command_path,jobs_path))
  
}

set_prediction_ready <- function(fPath,results_df,data){
  
  #Create necessary folders
  dir.create(paste(fPath,"/","prediction_input",sep=""),showWarnings=FALSE)
  dir.create(paste(fPath,"/","prediction_output",sep=""),showWarnings=FALSE)
  
  #Preparing the dataframe of conditions
  prediction_parameters <- results_df %>% 
    mutate(va="length_um") %>% 
    group_by(medium) %>% 
    arrange(log_likelihood) %>% 
    filter(row_number()==1) %>% #Here keeping only the parameters giving the best log likelihood.
    ungroup() %>% 
    arrange(medium) %>% 
    select(medium,va,ml,gamma,sl2,sm2,sd2) %>% 
    mutate(index=row_number()) %>% 
    mutate(input_path=sprintf("%s/prediction_input/%s_%s.csv",fPath,medium,as.character(index)),
           output_path=sprintf("%s/prediction_output/%s_%s.csv",fPath,medium,as.character(index)),
           lineToWrite=sprintf("python %s/optimization_code/path_predictions.py %s %s %s %s %s %s %s %s",fPath,input_path,output_path,va,ml,gamma,sl2,sm2,sd2))
  
  index_df <- prediction_parameters %>% 
    select(medium,index)
  
  #Copy paste data
  data %>%
    left_join(index_df,by=c("medium")) %>% 
    rename(cell_ori=cell) %>% 
    mutate(pos=paste(pos,orientation,sep=".")) %>%  #NECESSARY CAUSE ATHOS CODE RECOMPUTES WRONGLY THE LANE_ID...
    mutate(pos=paste(pos,"_",sep="")) %>% #NECESSARY CAUSE ATHOS CODE RECOMPUTES WRONGLY THE LANE_ID...
    mutate(gl=paste(gl,"_",sep="")) %>% #NECESSARY CAUSE ATHOS CODE RECOMPUTES WRONGLY THE LANE_ID...
    mutate(cell=paste(lane_ID,id,sep="_")) %>% 
    group_by(medium) %>% 
    do((function(.df){
      medium <- unique(.df$medium)
      readr::write_csv(.df,sprintf("%s/prediction_input/%s_%s.csv",fPath,medium,unique(.df$index)))
      return(data.frame())})(.))
  
  
  command_path <- sprintf("%s/prediction_input/commands.cmd",fPath)
  results_path <- sprintf("%s/prediction_output/results_%%a",fPath)
  errors_path <- sprintf("%s/prediction_output/errors_%%a",fPath)
  jobs_path <- sprintf("%s/optimization_code/prediction_jobs.sh",fPath)
  #system(sprintf('rm %s',command_path))
  system(sprintf('touch %s',command_path))
  
  #Adding jobs to command
  prediction_parameters %>% 
    arrange(index) %>% 
    group_by(index) %>% 
    do((function(.df){
      write(unique(.df$lineToWrite),file=sprintf('%s/prediction_input/commands.cmd',fPath),append=TRUE)
      return(data_frame())})(.)) %>% 
    ungroup()
  
  #Writing the prediction_jobs.sh
  command_path <- sprintf("%s/prediction_input/commands.cmd",fPath) 
  results_path <- sprintf("%s/prediction_output/results_%%a",fPath)
  errors_path <- sprintf("%s/prediction_output/errors_%%a",fPath)
  jobs_path <- sprintf("%s/optimization_code/prediction_jobs.sh",fPath)
  system(sprintf("sed -i \"s+path_to_output+%s+g\" %s",results_path,jobs_path))
  system(sprintf("sed -i \"s+path_to_error+%s+g\" %s",errors_path,jobs_path))
  system(sprintf("sed -i \"s+path_to_commands+%s+g\" %s",command_path,jobs_path))
  
  runAllJobsLine <- sprintf("sbatch --array=1-$(cat %s|wc -l):1 prediction_jobs.sh",command_path)
  print(runAllJobsLine)
  return(prediction_parameters)
}

dichotomic_search_growth_rate <- function(.var,.option){
  if(.option=="div_time"){
    .doubling_time <- .var # As a doubling times I am using first, the div times.
    .gr <- log(2)/mean(.var,na.rm=TRUE)
  }else if(.option=="alpha"){
    .doubling_time <- log(2)/.var # As a doubling times I am using first, the div times.
    .gr <- mean(.var,na.rm=TRUE)
  }
  p <- length(.var) #number of observations
  # As a first approximation for the population exponential growth rate, I am using mean(log(2)/.div_time), but I could also use .gr <- mean(alpha). 
  .approx_inf <- 1/2*mean(.gr,na.rm=TRUE) #low initial value
  .approx_sup <- 3*mean(.gr,na.rm=TRUE) #high initial value
  .delta_inf <- sum(exp(-.approx_inf*.doubling_time))-p/2 #Corresponding the function value
  .delta_sup <- sum(exp(-.approx_sup*.doubling_time))-p/2 #Corresponding the function value
  while(.delta_inf<0){
    .approx_inf <- 1/2*.approx_inf
    .delta_inf <- sum(exp(-.approx_inf*.doubling_time))-p/2} #This value should be POSITIVE to start with
  while(.delta_sup>0){
    .approx_sup <- 2*.approx_sup
    .delta_sup <- sum(exp(-.approx_sup*.doubling_time))-p/2} #This value should be NEGATIVE to start with
  # Simple dichotomic search for root
  .approx_new <- 1/2*(.approx_sup-.approx_inf)+.approx_inf #New growth rate value
  .delta_new <- sum(exp(-.approx_new*.doubling_time))-p/2 #New function value
  while((.approx_sup-.approx_inf)>0.00000001){ #Search should stop when difference between .approx_sup and .approx_inf becomes too small
    #print(.approx_sup)
    #print(.approx_inf)
    if(.delta_new<0){ # If function is negative, the intermediate value becomes the upper value, a new intermediate value is computed.
      .approx_sup <- .approx_new
      .delta_sup <- .delta_new
      .approx_new <- 1/2*(.approx_sup-.approx_inf)+.approx_inf
      .delta_new <- sum(exp(-.approx_new*.doubling_time))-p/2}
    else{ # Otherwise, the intermediate value becomes the inf value. New intermediate is computed.
      .approx_inf <- .approx_new
      .delta_inf <- sum(exp(-.approx_inf*.doubling_time))-p/2
      .approx_new <- 1/2*(.approx_sup-.approx_inf)+.approx_inf
      .delta_new <- sum(exp(-.approx_new*.doubling_time))-p/2}}
  #print(c(.approx_inf,.delta_inf,.approx_sup,.delta_sup))
  return(rep(c(.approx_sup+.approx_inf)/2,times=p))
}

