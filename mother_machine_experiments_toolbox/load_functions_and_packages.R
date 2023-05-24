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
#scale_color_manual(values=c25)+

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

library(pracma)
library(shiny)
library(grid)
library(plotly)
library(latex2exp)
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
library(RColorBrewer)
library(lubridate)
library(knitr)
library(parallel)
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


# ------------------------------------ Function to call cluster -----------------------------------------------------------------

load_cluster <- function(){
  if(!exists("mycluster")){
    mycluster <- (parallel::detectCores()-1) %>%  # do not use more than 30 cores
      default_cluster() %>% cluster_library( # load currently loaded packages on each core
        names(sessionInfo()$otherPkgs))
  }
  return(mycluster)}

copy_cluster <- function(cl,.v){
  cluster_copy(cl,.v)
}

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


######### Simulations and inference useful functions #########

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

compute_weighted_mean <- function(.m,.sd){
  p <- length(.m)
  .w <- 1/(.sd)**2
  .w_tot <- sum(.w)
  .w_m <- sum(.m*.w)/.w_tot
  return(rep(c(.w_m),times=p))
}

compute_weighted_sd <- function(.m,.sd){
  p <- length(.m)
  .w <- 1/(.sd)**2
  .w_tot <- sum(.w)
  .w_m <- sum(.m*.w)/.w_tot
  .w_sd <- sqrt(sum(.w*(.m-.w_m)**2)/.w_tot)
  return(rep(c(.w_sd),times=p))
}

#-------------------------------------------------- Plotting functions --------------------------------------------------


plot_dataset_statistics <- function(.condition,.promoter,.full_cell_cycle_only){
  
  if(.full_cell_cycle_only==TRUE){
    data <- myframes_complete
  }else{
    data <- myframes}
  
  myframes_toplot <- data %>% 
    filter(condition==.condition) %>% 
    filter(promoter %in% c(.promoter)) %>% 
    extract(gl_id,"gl_number","^20[0-9]{6}.[0-9]{1,}.([0-9]{1,})$",remove="FALSE",convert="FALSE") %>% 
    mutate(gl_number=as.double(gl_number)) %>% 
    #mutate(orientation=ifelse(gl_number>30,"right","left")) %>% 
    distinct(cell,.keep_all=TRUE)
  
  myframes_toplot$gl_number <- factor(myframes_toplot$gl_number,levels=sort(unique(myframes_toplot$gl_number)))
  myframes_toplot$mean_cell_rank <- factor(myframes_toplot$mean_cell_rank,levels=as.character(sort(unique(as.integer(myframes_toplot$mean_cell_rank)))))
  
  print(myframes_toplot %>% 
          ggplot() +
          geom_bar(aes(gl_number,fill=mean_cell_rank),position="dodge")+
          theme_cowplot()+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
          #coord_cartesian(ylim=c(1,3))+
          #labs(subtitle="Right growth lanes (bottom to top)")+
          facet_wrap(~paste("Position",pos,sep=" "),scale="free_x",nrow=1)+
          scale_fill_manual(values=c25)+
          guides(fill = guide_legend(title = "Cell rank"))+
          #theme(legend.position = "none")+
          xlab("Growth lane ID")+
          ylab("Number of cells"))}

plot_dataset_statistics_curated <- function(.condition,.promoter,.full_cell_cycle_only){
  
  if(.full_cell_cycle_only==TRUE){
    data <- myframes_complete_goodcells
  }else{
    data <- myframes %>% semi_join(mydf,by=c("cell"))}
  
  myframes_toplot <- data %>% 
    filter(condition==.condition) %>% 
    filter(promoter %in% c(.promoter)) %>% 
    extract(gl_id,"gl_number","^20[0-9]{6}.[0-9]{1,}.([0-9]{1,})$",remove="FALSE",convert="FALSE") %>% 
    mutate(gl_number=as.double(gl_number)) %>% 
    #mutate(orientation=ifelse(gl_number>30,"right","left")) %>% 
    distinct(cell,.keep_all=TRUE)
  
  myframes_toplot$gl_number <- factor(myframes_toplot$gl_number,levels=sort(unique(myframes_toplot$gl_number)))
  myframes_toplot$mean_cell_rank <- factor(myframes_toplot$mean_cell_rank,levels=as.character(sort(unique(as.integer(myframes_toplot$mean_cell_rank)))))
  
  print(myframes_toplot %>% 
          ggplot() +
          geom_bar(aes(gl_number,fill=mean_cell_rank),position="dodge")+
          theme_cowplot()+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
          #coord_cartesian(ylim=c(1,3))+
          #labs(subtitle="Right growth lanes (bottom to top)")+
          facet_wrap(~paste("Position",pos,sep=" "),scale="free_x",nrow=1)+
          scale_fill_manual(values=c25)+
          guides(fill = guide_legend(title = "Cell rank"))+
          #theme(legend.position = "none")+
          xlab("Growth lane ID")+
          ylab("Number of cells"))}

plot_concentration_distribution <- function(.date,.gl_min,.gl_max){
  
  .n <- nrow(myframes_complete_goodcells %>%
               semi_join(good_cells,by=c("cell")) %>% 
               ungroup() %>% 
               filter(date==.date) %>% 
               distinct(condition,promoter,pos))
  
  p <- myframes_complete_goodcells %>% 
    filter(date==.date) %>% 
    mutate(concentration=gfp_nb/length_um) %>% 
    group_by(cell) %>% 
    mutate(mean_concentration=mean(concentration,na.rm=TRUE)) %>% 
    ungroup() %>% 
    distinct(cell,.keep_all=TRUE) %>% 
    ggplot()+
    geom_boxplot(aes(gl_number,mean_concentration,group=gl_id))+
    facet_wrap(~paste(condition,promoter,pos,sep=" "),ncol=1,scales="free_y")+
    coord_cartesian(xlim=c(.gl_min,.gl_max))
  
  #ggsave(sprintf("./concentration_distributions_%s.pdf",.date),device="pdf",width=15,height=.n*3,limitsize=FALSE,plot=p)
  #ggsave("./concentration_distributions.pdf",device="pdf",width=15,height=.n*3,limitsize=FALSE,plot=p)
  #ggsave(sprintf("./concentration_distributions_%s.pdf",.date),device="pdf",width=15,height=.h,limitsize=FALSE,plot=p)
  
  print(p)
  
}

# --------------------------------------------------------- Random sampling -------------------------------------------------------------

random_selection_no_replacement <- function(.d,.N,.N0){
  #.d: input dataframe
  #.N: Number of observations to pick
  #.N0 number of 'bags'
  
  .df <- .d %>% 
    ungroup() %>% 
    #group_by(cell) %>%  #already grouped by condition, promoter, type, simulationid
    #sample_n(1) %>% 
    #mutate(lambda_cell=mean(lambda*60/log(2),na.rm=TRUE)) %>% 
    #distinct(cell,.keep_all = TRUE) %>% 
    sample_n(.N,replace=FALSE) %>% #no replacement
    mutate(sampling=1)
  
  for(i in 2:.N0){
    .df_new <- .d %>% 
      ungroup() %>% 
      #group_by(cell) %>%  #already grouped by condition, promoter, type, simulationid
      #sample_n(1) %>% 
      #mutate(lambda_cell=mean(lambda*60/log(2),na.rm=TRUE)) %>% 
      #ungroup() %>% 
      #distinct(cell,.keep_all = TRUE) %>% 
      sample_n(.N,replace=TRUE) %>% 
      mutate(sampling=i)
    .df <- rbind(.df,.df_new)}
  
  return(.df)}

random_selection_per_time_bin <- function(.d,.N0){
  .df <- .d %>% 
    ungroup() %>% 
    mutate(data_group=1) %>% 
    #already grouped by condition, promoter, date,time_bin
    group_by(gl_id) %>%
    sample_n(1) %>% 
    ungroup()
  
  for(i in 2:.N0){
    .df_new <- .d %>% 
      ungroup() %>% 
      mutate(data_group=i) %>% 
      #group_by(cell) %>%  #already grouped by condition, promoter, type, simulationid
      group_by(gl_id) %>% 
      sample_n(1) %>% 
      ungroup()
    .df <- rbind(.df,.df_new)}
  return(.df)}

# --------------------------------------------------------- Custom distribution routine -------------------------------------------------------------

erf_custom <- function(x){ 
  2 * pnorm(x * sqrt(2)) - 1}

compute_breaks <- function(.means,.sds,.N_breaks){
  .min <- min(.means-15*.sds)
  .max <- max(.means+4*.sds)
  .breaks <- seq(.min,.max,length.out=.N_breaks)
  #.means-3sigma*sd should ALL be strictly positive, filter should be implemented there
  .means_mins <- .means-3*.sds
  .means_log <- .means[.means_mins>0]
  .sds_log <- .sds[.means_mins>0]
  .min <- min(.means_log-3*.sds_log)
  .max <- max(.means_log+3*.sds_log)
  k <- (.max/.min)**(1/.N_breaks)
  .breaks_log <- unlist(lapply(seq(0,.N_breaks), function(.n) .min*k^(.n)))
  return(list(.breaks,.breaks_log))}

compute_hist_single_mean <- function(mean,sd,.breaks){
  .l_breaks <- length(.breaks)
  .left_breaks <- .breaks[1:(.l_breaks-1)]
  .right_breaks <- .breaks[2:.l_breaks]
  hist_count <- mapply(compute_error_function,mean,sd,.left_breaks,.right_breaks)
  return(hist_count)
}

compute_error_function <- function(mean,sd,.left_break,.right_break){
  if(sd==0){
    return(1*(data.table::between(mean,.left_break,.right_break)))
  }else{
    left_threshold <- mean-3*sd
    right_threshold <- mean+3*sd
    a <- (.left_break < left_threshold & .right_break < left_threshold)
    if(!(a %in% c(FALSE,TRUE))){
      print(.left_break)
      print(left_threshold)
      print(.right_break)
      print(right_threshold)}
    if(.left_break < left_threshold & .right_break < left_threshold){
      return(0)}
    else if(.left_break > right_threshold & .right_break > right_threshold){
      return(0)}
    else{
      .left_break <- max(.left_break,left_threshold) 
      .right_break <- min(.right_break,right_threshold)
      result <- 1/2*(erf_custom((.right_break-mean)/(sqrt(2)*sd))-erf_custom((.left_break-mean)/(sqrt(2)*sd)))
      return(result)}}}

compute_hist <- function(.means,.sds,.breaks){
  .l_breaks <- length(.breaks)
  .centers <- (.breaks[2:.l_breaks]-.breaks[1:(.l_breaks-1)])/2 + .breaks[1:(.l_breaks-1)]
  hist_count <- mapply(compute_hist_single_mean,.means,.sds,list(.breaks),SIMPLIFY = FALSE)
  #count
  hist_count_final <- Reduce("+", hist_count)
  #density
  .s <- sum(hist_count_final)*(.breaks[2]-.breaks[1])
  hist_density_final <- hist_count_final/.s
  #logx - count
  .breaks_non_zeros <- .breaks[.breaks>0]
  .centers_non_zeros <- .centers[.breaks[1:(.l_breaks-1)]>0]
  hist_count_logx <- hist_count_final[.breaks[1:(.l_breaks-1)]>0]
  .l_breaks <- length(.breaks_non_zeros)
  .l_w <- log(.breaks_non_zeros[2:.l_breaks])-log(.breaks_non_zeros[1:(.l_breaks-1)])
  hist_count_logx_final <- hist_count_logx/.l_w
  #logx - density
  .s <- sum(hist_count_logx)
  hist_density_logx_final <- hist_count_logx_final/.s
  return(list(.breaks,
              .centers,
              hist_count_final,
              hist_density_final))
}

##### Plotting functions testing routine

#testing erf in the context of a gaussian with mean 0 and sd 1
# a <- 1/2*(erf(1/sqrt(2))-erf(-1/sqrt(2))) = 68.2
# a <- 1/2*(erf_custom(1/sqrt(2))-erf(-1/sqrt(2)))
# a <- 1/2*(erf_custom(2/sqrt(2))-erf(-2/sqrt(2)))
# compute_gaussian <- function(mean,sd,center,width){
#  if(abs((mean-center))>3*sd){return(0)
#  }else{return((width)*1/(sd*sqrt(2*pi))*exp(-(mean-center)**2/(2*(sd**2))))}}
# 
# .means <- noise_contribution_data %>% filter(condition=="acetate005",promoter=="hi1",date=="20230228",datatype=="allnoisesources") %>% 
#   mutate(means=length_um,
#          sds=1) %>% 
#   .$means
# 
# .sds <- noise_contribution_data %>% filter(condition=="acetate005",promoter=="hi1",date=="20230228",datatype=="allnoisesources") %>% 
#   mutate(means=length_um,
#          sds=1) %>% 
#   .$sds
# 
# N_breaks <- 200
# .breaks <- compute_breaks(.means,.sds,N_breaks)
# .breaks_linear <- unlist(.breaks[1])
# .hist <- compute_hist(.means,.sds,.breaks_linear)
# 
# .breaks_log <- unlist(.breaks[2])
# .means_mins <- .means-3*.sds
# .means_log <- .means[.means_mins>0]
# .sds_log <- .sds[.means_mins>0]
# .hist_log <- compute_hist(.means_log,.sds_log,.breaks_log)
# 
# plot(unlist(.hist[2]),unlist(.hist[3]))
# plot(unlist(.hist[2]),unlist(.hist[4]))
# plot(unlist(.hist_log[2]),unlist(.hist_log[3]),log='x')
# plot(unlist(.hist_log[2]),unlist(.hist_log[4]),log='x')

# --------------------------------------------------------- Import data routine -------------------------------------------------------------


prepare_simulation_csv <- function(.l){
  
  promoter <- str_match(.l,"/[0-9a-zA-Z]{1,}_([0-9a-zA-Z]{1,})_[0-9]{1,}_[a-z]{1,}.csv$")[2]
  condition <- str_match(.l,"/([0-9a-zA-Z]{1,})_[0-9a-zA-Z]{1,}_[0-9]{1,}_[a-z]{1,}.csv$")[2]
  carbon_source_concentration <- str_match(condition,"^[a-z]{1,}0([0-9]{2})")[2]
  carbon_source_concentration=as.double(carbon_source_concentration)/100
  carbon_source <- str_match(condition,"^([a-z]{1,})0[0-9]{2}")[2]
  datatype <- str_match(.l,"/[0-9a-z]{1,}_[0-9a-zA-Z]{1,}_[0-9]{1,}_([a-z]{1,}).csv$")[2]
  date <- str_match(.l,"/[0-9a-z]{1,}_[0-9a-zA-Z]{1,}_([0-9]{1,})_[a-z]{1,}.csv$")[2]
  
  .df <- readr::read_csv(.l) %>% 
    mutate(path=.l) %>%
    rename(cell=cell_id,parent=parent_id,gfp_nb=gfp,lambda=lt,q=qt) %>% 
    mutate(time_sec=time_min*60,length_um=exp(log_length)) %>% 
    mutate(promoter=promoter,condition=condition,carbon_source_concentration=carbon_source_concentration,
           carbon_source=carbon_source, datatype=datatype, date=date)
  
  return(.df)
}


