# # Initial setup of the folder
setwd("~/Documents/GitHub/NhaTrang_flu_rsv_interaction/")
library(here)

library(parallel)
library(doParallel)
library(reshape2)
library(MASS)
library(deSolve)
library(ggplot2)
library(data.table)
library(tmvtnorm)
library(Rcpp)
library(RColorBrewer)
library(gridExtra)
library(zoo)
library(forcats)
library(Momocs)
library(GGally)
library(tmvtnorm)
library(lubridate)

# # These needed always
source(here::here("fit", "compiler_function.R"))
source(here::here("fit","NT_fitting_mcmc.R"))
source(here::here("fit","NT_functions_dualmult.R"))
source(here::here("fit","NT_together_parameters.R"))
source(here::here("fit","NT_together_Rrats.R"))
# compile and load the cpp model
compileModel(here::here("fit","NT_model_dual.cpp"), "fit/build")
dyn.load("fit/build/NT_model_dual.so")
#install.packages("RcppPHE_0.1.0.tar.gz", repos=NULL, type="source")
#library(RcppPHE)

population_contacts <- as.matrix(read.csv(here::here("fit","Mixing_Matrix_model.csv"), header = T, row.names = 1))
population_numbers <- c(5230, 9070, 35654, 148255, 12530 )

#NOTE - this data is not publically available
multi_data <- readRDS(here::here("AWS","multi_data_dual_JULY.RDS"))

# open the fitting main script
file.edit(here::here("fit","AWS_main_run.R"))


########## visualise the fit ########
source(here::here("figures","visualise_funcs.R"))
source(here::here("figures","quantile_plot_functions.R"))
source(here::here("figures","quantile_plots_dual.R"))
source(here::here("figures","vaccine_testting.R"))
source(here::here("figures", "figure_generation_functions.R"))
file.edit(here::here("figures", "figure_generation.R"))



