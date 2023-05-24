# First we want to get rid of control wells that are potentially contaminated
predict_od <- function(.c_od,.t_min){
  if (length(.c_od) != length(.t_min)) 
    stop("variables of different lengths.")
  if (length(.c_od) <= 2) 
    return(NA)
  stats::lm(log(.c_od) ~ .t_min) %>% stats::predict(se.fit = TRUE) %>% 
    .[["fit"]] %>% exp
}

predict_alpha <- function(.p_od,.t_min){
  if (length(.p_od) != length(.t_min)) 
    stop("variables of different lengths.")
  if (length(.p_od) <= 2) 
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

predict_intercept_fluo_od <- function(.c_f,.c_od){
  if (length(.c_f) != length(.c_od)) 
    stop("variables of different lengths.")
  if (length(.c_f) <= 2) 
    return(NA)
  stats::lm(log(.c_f) ~ log(.c_od)) %>% stats::coef(se.fit = TRUE) %>% 
    .[[1]] %>% exp %>% rep(times=length(.c_f))
}

predict_power_fluo_od <- function(.c_f,.c_od){
  if (length(.c_f) != length(.c_od)) 
    stop("variables of different lengths.")
  if (length(.c_f) <= 2) 
    return(NA)
  stats::lm(log(.c_f) ~ log(.c_od)) %>% stats::coef(se.fit = TRUE) %>% 
    .[[2]] %>% rep(times=length(.c_f))
}

.c_f <- c(1,2,3)
.c_od <- c(exp(0.5)*1,exp(0.5)*4,exp(0.5)*9)
test <- stats::lm(log(.c_od) ~ log(.c_f)) %>% stats::predict(se.fit = TRUE)
test <- stats::lm(log(.c_od) ~ log(.c_f)) %>% stats::coef(se.fit = TRUE) %>% .[[1]] %>% exp %>% rep(times=length(.c_f))

