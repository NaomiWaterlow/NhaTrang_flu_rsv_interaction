# NT together core transformations

# Core script for paralllel tempering of NT data.

array_num <- 16

library(doParallel)
library(reshape2)
library(MASS)
library(deSolve)
library(data.table)
library(tmvtnorm)
library(Rcpp)

########## Load in params and functions and initialise paralisation

# Specify the characteristics of the cluster
cores_to_use <- 3
cl <- parallel::makeCluster(cores_to_use, setup_strategy = "sequential"#, outfile="tester.txt"
)
registerDoParallel(cl)

###### Load in the other scripts ########

source("NT_fitting_mcmc.R")
source("NT_together_parameters.R")
source("NT_together_Rrats.R")
dyn.load("build/NT_model_dual.so")
population_contacts <- as.matrix(read.csv("Mixing_Matrix_model.csv", header = T, row.names = 1))
population_numbers <- c(5230, 9070, 35654, 148255, 12530 )
multi_data <- readRDS("multi_data_dual_JULY.RDS")

source("NT_functions_dualmult.R")

############ CONDITIONS AND RUN ############

input_params_theta <- read.csv("input_params_table_transform.csv")
input_params_proposal <- read.csv("input_proposal_table_transform.csv")
theta_bounds <- read.csv("theta_bounds.csv",row.names = 1)
theta_bounds[which(theta_bounds$Lower == Inf), "Lower"] <- - Inf
lower_bounds = theta_bounds[,"Lower"]
upper_bounds = theta_bounds[,"Upper"]
names(lower_bounds) =names(upper_bounds) = row.names(theta_bounds)
lower_bounds <- c(lower_bounds, overdispersion = 0)
upper_bounds <- c(upper_bounds, overdispersion = Inf)
lower_bounds <- c(lower_bounds, Dual_mult = 1)
upper_bounds <- c(upper_bounds, Dual_mult = Inf)

init_theta <- input_params_theta[,array_num]
names(init_theta) <- input_params_theta[,1]
init_theta <- c(init_theta, overdispersion = 1)
init_theta <- c(init_theta, Dual_mult = 10)

proposal_sd <- input_params_proposal[,array_num]/20
names(proposal_sd) <- input_params_proposal[,1]
proposal_sd <- c(proposal_sd, overdispersion = 0.01)
proposal_sd <- c(proposal_sd, Dual_mult = 0.1)

covmat <- matrix(0,ncol=length(init_theta), nrow=length(init_theta))
diag(covmat) <- proposal_sd
colnames(covmat) = rownames(covmat) = names(init_theta)

covmat <- as.matrix(nearPD(readRDS("20210708_dualmult.RDS"))[[1]])/10
propsal_sd <- diag(covmat)

R0ratI <- calculate_R0_rat_I(init_theta, population_numbers)
R0ratR <- calculate_R0_rat_R(init_theta, population_numbers)

Temperatures <- c(1,1.4,1.9,2.5,4,8,14,25,50,150,500,1000)
n_chains <- 12

load("Pop_dualmult_cont_16_2021-08-08.Rdata")

# Export parameters to the cores
clusterExport(cl, list("n_chains","population_numbers",
                       "multi_data", "R0ratR", "R0ratI", "covmat", "Temperatures",
                       "population_contacts", "gammaR", "gammaI", "rsv_sus", "season_params", 
                       "lower_bounds", "upper_bounds"))
# Export required functions to the cores
clusterExport(cl, list( "calculate_init_state",
                        "calc_priors", "calc_likelihood_overall",
                        "simulate_SIPR","summary_stats",
                        "temp_func", "step_fn", "swapping_rates", "subset_parameters", 
                        "adjust_unsymmetrical", "pars_transformation", 
                        "calc_likelihood_overall",
                        "calc_priors"))
# Export relevant packages to the cores
clusterEvalQ(cl, library("Rcpp"))
clusterEvalQ(cl, library("data.table"))
clusterEvalQ(cl, library("deSolve"))
clusterEvalQ(cl, library("tmvtnorm"))
clusterEvalQ(cl, dyn.load("build/NT_model_dual.so"))

time_start <- Sys.time()
name <- "Pop_dualmult_cont2"
total_trace <- trace_wrapper_cont(trace_previous = total_trace,
                               additional_run = 100000,
                               old_length_run = 100000,
                               each_run = 500,
                               n_chains = n_chains,
                               data_to_fit = multi_data,
                               population_numbers = population_numbers,
                               indep = 5,
                               covmat = covmat,
                               adapt_rate = 0.5,
                               adaption_starter = 10090000,
                               proposal_sd = proposal_sd,
                               Temperatures = Temperatures,
                               name = name,
                               save = T
)

time_end <- Sys.time()

print(time_end - time_start)
save(total_trace, file=paste0(name,"_", array_num,"_", Sys.Date() , ".Rdata"))
stopCluster(cl)


