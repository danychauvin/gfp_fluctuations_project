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

# Example of applications
#mydata %>%
#  group_by(plate,row,column) %>% 
#  arrange(time_min) %>% 
#  mutate(predicted_od=predict_od(corrected_od,time_min)) %>%
#  mutate(alpha=predict_alpha(predicted_od,time_min)) %>% # 
#  ungroup() %>% 
#  filter(condition=="glucose") %>% 
#  ggplot()+
#  geom_point(aes(time_min,corrected_od))+
#  geom_point(aes(time_min,predicted_od),col="red")+
#  facet_wrap(~interaction(plate,row,column),scales="free")+
#  theme_cowplot()


#predict_fluo_od #######################################################################
#input: corrected_fluo,corrected_od
#produces a linear fit between log(fluo) and log(od) (power law of the type fluo = A * od ^ B)
#output: fluo=A * od ^ B

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


