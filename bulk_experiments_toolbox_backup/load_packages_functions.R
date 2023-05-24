# 20201211, Dany Chauvin

# INSTALLING AND LOADING NECESSARY PACKAGES 

# Following packages are necessary to import data from deepMoma and analyze these data.

# WARNING: libstdc++.so.6 library, with proper GLIBCXX versions is necessary to use some packages imported here. 
# Be sure that GCCcore/8.3.0 is loaded/installed.
# To do so while using Rstudio on the scicore server "service06", do the following, add the following to your ~/.bashrc file:
# `if [[ "$HOSTNAME" = *service06* ]]; then
#     ml GCCcore/8.3.0
#  fi'

# Currently, packages are installed in my HOME folder, from centOS7 binaries at packagemanager.rstudio.com. I am not yet using renv, planned in the future.
# On Rstudio, set R version to 4.0.3.
# If this is the first time you use these codes, uncomment to following lines to install the corresponding packages.

# Uncomment below when using for the first time.
#options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/__linux__/centos7/latest"))
#install.packages("tidyverse")
#install.packages(c("here","cowplot"))
#install.packages("devtools")
#remotes::install_github(c('hadley/multidplyr'))
#install.packages("ccaPP")
#install.packages("~/R/vngMoM.tar", repos = NULL)
#remotes::install_github(c('julou/ggCustomTJ'))
#renv::init()
#install.packages("reticulate")
#install.packages("ggcorrplot")
#install.packages("lemon")
#install.packages("LSD")
#install.packages("parallel")
#install.packages("ggpubr")
# Comment above if the packages are already installed.

# Necessary libraries
library(tidyverse)
library(RcppArmadillo)
library(tools)
library(here)
library(cowplot)
library(devtools)
library(multidplyr)
library(vngMoM)
library(ggCustomTJ)
library(renv)
library(svglite)
#Sys.setenv(RETICULATE_PYTHON = '/scicore/home/nimwegen/rocasu25/.local/share/r-miniconda/envs/r-reticulate/bin/python')
library(reticulate)
library(ggcorrplot)
library(lemon)
library(parallel)
library(broom)
library(stats)
library(ggpubr)
library(comprehenr)

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



