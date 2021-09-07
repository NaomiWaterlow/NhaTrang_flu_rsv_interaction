# Functions for fitting the NT model, with dual cases over the season.
pars_transformation <- function(parameters, n_chains){
  
  if (parameters[["sig"]] == 0 ){
    parameters[["sig"]] = 1
  } else {
    parameters[["sig"]] = (1 - parameters[["sig"]])
  }
  
  parameters["rho"] <- exp(as.numeric(parameters["rho"])) 
  
  parameters["propR_1"] <- as.numeric(parameters["propR_1"])/10000
  parameters["propR_2"] <- as.numeric(parameters["propR_2"])/10000
  parameters["propR_3"] <- as.numeric(parameters["propR_3"])/10000
  parameters["propR_4"] <- as.numeric(parameters["propR_4"])/10000
  parameters["propR_5"] <- as.numeric(parameters["propR_5"])/10000
  parameters["propR_6"] <- as.numeric(parameters["propR_6"])/10000
  parameters["propR_7"] <- as.numeric(parameters["propR_7"])/10000
  parameters["propR_8"] <- as.numeric(parameters["propR_8"])/10000
  parameters["propR_9"] <- as.numeric(parameters["propR_9"])/10000
  parameters["propR_10"] <- as.numeric(parameters["propR_10"])/10000
  parameters["propR_11"] <- as.numeric(parameters["propR_11"])/10000
  
  parameters["propI_1"] <- as.numeric(parameters["propI_1"])/10000
  parameters["propI_2"] <- as.numeric(parameters["propI_2"])/10000
  parameters["propI_3"] <- as.numeric(parameters["propI_3"])/10000
  parameters["propI_4"] <- as.numeric(parameters["propI_4"])/10000
  parameters["propI_5"] <- as.numeric(parameters["propI_5"])/10000
  parameters["propI_6"] <- as.numeric(parameters["propI_6"])/10000
  parameters["propI_7"] <- as.numeric(parameters["propI_7"])/10000
  parameters["propI_8"] <- as.numeric(parameters["propI_8"])/10000
  parameters["propI_9"] <- as.numeric(parameters["propI_9"])/10000
  parameters["propI_10"] <- as.numeric(parameters["propI_10"])/10000
  parameters["propI_11"] <- as.numeric(parameters["propI_11"])/10000
  
  parameters["RSV_default"] <- 1
  
  return(parameters)
}


calc_likelihood_overall <- function (theta_season, init_state_season, season_data, log = FALSE)
{
  
  # times to run the model over
  theta_season[["gammaR"]] = gammaR
  theta_season[["gammaI"]] = gammaI
  theta_season[["Contact_Structure"]] = population_contacts
  theta_season[["num_grps"]] = 5
  theta_season[["RSV_Sus"]] = rsv_sus
  theta_season[["seedR"]] = 0
  theta_season[["seedI"]] = 0
  theta_season[["pop_sizes"]] = unname(population_numbers)
  theta_season[["trickleI"]] = 0
  theta_season[["trickleR"]] = 0
  
  times <- seq(from = 1, to = nrow(season_data)*7, by = 1)
  
  trajectory <- simulate_SIPR(theta_season, unname(init_state_season), times)
  # set up density calculators
  lik_table <- cbind(season_data,trajectory)
  # reporting rates
  
  lik_table$INF_modelcases_1 <- lik_table$INF_incidence_1*theta_season[["Idetect"]]
  lik_table$INF_modelcases_2 <- lik_table$INF_incidence_2*theta_season[["Idetect"]]*theta_season[["detect_multiplier"]]
  lik_table$RSV_modelcases_1 <- lik_table$RSV_incidence_1*theta_season[["Rdetect2"]]*theta_season[["RSV_mult"]]
  lik_table$RSV_modelcases_2 <- lik_table$RSV_incidence_2*theta_season[["Rdetect5"]]*theta_season[["RSV_mult"]]
  lik_table$Dual_modelcases_1 <- lik_table$Dual_incidence_1*theta_season[["Rdetect2"]]*theta_season[["RSV_mult"]]*theta_season[["Dual_mult"]]
  lik_table$Dual_modelcases_2 <- lik_table$Dual_incidence_2*theta_season[["Rdetect5"]]*theta_season[["RSV_mult"]]*theta_season[["Dual_mult"]]
  lik_table$all_0 <- lik_table$Hos_INF_0 + lik_table$Hos_RSV_0 + lik_table$both_0
  lik_table$all_1 <- lik_table$Hos_INF_1 + lik_table$Hos_RSV_1 + lik_table$both_1
  lik_table$all_modelcases_0 <-  lik_table$INF_modelcases_1+ lik_table$RSV_modelcases_1+lik_table$Dual_modelcases_1
  lik_table$all_modelcases_1 <-  lik_table$INF_modelcases_2+ lik_table$RSV_modelcases_2+lik_table$Dual_modelcases_2
  
  
  dens <- 0
  dens_multinom <- 0
  dens_negbinom <- 0
  
  for(j in 1:nrow(lik_table)){
    
    young_multinom <-  dmultinom(x = unlist(lik_table[j,c("Hos_INF_0", "Hos_RSV_0", "both_0")]), 
                                 prob = unlist(lik_table[j,c("INF_modelcases_1","RSV_modelcases_1","Dual_modelcases_1")]), 
                                 log = T) 
    
    old_multinom <-   dmultinom(x = unlist(lik_table[j,c("Hos_INF_1", "Hos_RSV_1", "both_1")]), 
                                prob = unlist(lik_table[j,c("INF_modelcases_2","RSV_modelcases_2","Dual_modelcases_2")]), 
                                log = T)
    
    young_negbinom <- dnbinom(x = lik_table[[j,"all_0"]],
                              mu = lik_table[[j,"all_modelcases_0"]],
                              size = theta_season[["overdispersion"]]  # Infinity is a poisson
                              , log=T)
    old_negbinom <- dnbinom(x = lik_table[[j,"all_1"]],
                            mu = lik_table[[j,"all_modelcases_1"]],
                            size = theta_season[["overdispersion"]]  # Infinity is a poisson
                            , log=T)
    
    dens_multinom <- dens_multinom + young_multinom + old_multinom
    dens_negbinom <- dens_negbinom + young_negbinom + old_negbinom 
    
  }
  dens <- dens_multinom + dens_negbinom
  # print(paste0("likelihood mulitnomial is ", dens_multinom))
  # print(paste0("likelihood negbinom is ",dens_negbinom))
  
  return(ifelse(log, dens, exp(dens)))
}


## Calculate the priors
calc_priors <- function (theta_season, R0ratI, R0ratR, log = FALSE)
{ # Calculate the prior for each parameters
  
  log_prior_betaR <- dunif((as.numeric(theta_season["betaR"])*R0ratR), min = 0.5, max=8, log=T)
  log_prior_propR1<- dunif(as.numeric(theta_season["propR_1"]), min=0, max=1, log=T)
  log_prior_propI1<- dunif(as.numeric(theta_season["propI_1"]), min=0, max=1, log=T)
  log_prior_propR2<- dunif(as.numeric(theta_season["propR_2"]), min=0, max=1, log=T)
  log_prior_propI2<- dunif(as.numeric(theta_season["propI_2"]), min=0, max=1, log=T)
  log_prior_propR3<- dunif(as.numeric(theta_season["propR_3"]), min=0, max=1, log=T)
  log_prior_propI3<- dunif(as.numeric(theta_season["propI_3"]), min=0, max=1, log=T)
  log_prior_propR4<- dunif(as.numeric(theta_season["propR_4"]), min=0, max=1, log=T)
  log_prior_propI4<- dunif(as.numeric(theta_season["propI_4"]), min=0, max=1, log=T)
  log_prior_propR5<- dunif(as.numeric(theta_season["propR_5"]), min=0, max=1, log=T)
  log_prior_propI5<- dunif(as.numeric(theta_season["propI_5"]), min=0, max=1, log=T)
  log_prior_propR6<- dunif(as.numeric(theta_season["propR_6"]), min=0, max=1, log=T)
  log_prior_propI6<- dunif(as.numeric(theta_season["propI_6"]), min=0, max=1, log=T)
  log_prior_propR7<- dunif(as.numeric(theta_season["propR_7"]), min=0, max=1, log=T)
  log_prior_propI7<- dunif(as.numeric(theta_season["propI_7"]), min=0, max=1, log=T)
  log_prior_propR8<- dunif(as.numeric(theta_season["propR_8"]), min=0, max=1, log=T)
  log_prior_propI8<- dunif(as.numeric(theta_season["propI_8"]), min=0, max=1, log=T)
  log_prior_propR9<- dunif(as.numeric(theta_season["propR_9"]), min=0, max=1, log=T)
  log_prior_propI9<- dunif(as.numeric(theta_season["propI_9"]), min=0, max=1, log=T)
  log_prior_propR10<- dunif(as.numeric(theta_season["propR_10"]), min=0, max=1, log=T)
  log_prior_propI10<- dunif(as.numeric(theta_season["propI_10"]), min=0, max=1, log=T)
  log_prior_propR11<- dunif(as.numeric(theta_season["propR_11"]), min=0, max=1, log=T)
  log_prior_propI11<- dunif(as.numeric(theta_season["propI_11"]), min=0, max=1, log=T)
  log_prior_betaI <- dunif((as.numeric(theta_season["betaI"])*R0ratI), min=0.8, max=4,log=T)
  log_prior_Rdetect2 <- dunif(as.numeric(theta_season["Rdetect2"]), min=0, max=0.5, log=TRUE)
  log_prior_Rdetect5 <- dunif(as.numeric(theta_season["Rdetect5"]), min=0, max=0.5, log=TRUE)
  log_prior_Idetect <- dunif(as.numeric(theta_season["Idetect"]), min=0, max=0.5, log=TRUE)
  log_prior_sigma <- dunif(as.numeric(theta_season["sig"]), min=0, max=1,log=T)
  #log_prior_sigma <- dnorm(as.numeric(theta_season["sig"]), mean=0.8, sd=0.5,log=T)
  log_prior_rho <- dunif(as.numeric(theta_season["rho"]), min=0, max = 0.5,log=T)
  log_prior_inf_sus1 <- dunif(as.numeric(theta_season["inf_sus_1"]), min=0, max=1, log=T)
  log_prior_inf_sus2 <- dunif(as.numeric(theta_season["inf_sus_2"]), min=0, max=1, log=T)
  log_prior_inf_sus3 <- dunif(as.numeric(theta_season["inf_sus_3"]), min=0, max=1, log=T)
  log_prior_inf_sus4 <- dunif(as.numeric(theta_season["inf_sus_4"]), min=0, max=1, log=T)
  log_prior_inf_sus5 <- dunif(as.numeric(theta_season["inf_sus_5"]), min=0, max=1, log=T)
  log_prior_inf_sus6 <- dunif(as.numeric(theta_season["inf_sus_6"]), min=0, max=1, log=T)
  log_prior_inf_sus7 <- dunif(as.numeric(theta_season["inf_sus_7"]), min=0, max=1, log=T)
  log_prior_inf_sus8 <- dunif(as.numeric(theta_season["inf_sus_8"]), min=0, max=1, log=T)
  log_prior_inf_sus9 <- dunif(as.numeric(theta_season["inf_sus_9"]), min=0, max=1, log=T)
  log_prior_inf_sus10 <- dunif(as.numeric(theta_season["inf_sus_10"]), min=0, max=1, log=T)
  log_prior_inf_sus11 <- dunif(as.numeric(theta_season["inf_sus_11"]), min=0, max=1, log=T)
  
  # Sum the priors for individual parameters
  log_sum <- log_prior_betaR + log_prior_betaI +
    log_prior_Rdetect2 + log_prior_Rdetect5 + log_prior_Idetect +
    log_prior_propR1 +  log_prior_propR2 +  log_prior_propR3 +  log_prior_propR4 +  log_prior_propR5 + 
    log_prior_propR6 +  log_prior_propR7 +  log_prior_propR8 +  log_prior_propR9 +  log_prior_propR10 +
    log_prior_propR11 +
    log_prior_propI1 +  log_prior_propI2 +  log_prior_propI3 +  log_prior_propI4 +  log_prior_propI5 + 
    log_prior_propI6 +  log_prior_propI7 +  log_prior_propI8 +  log_prior_propI9 +  log_prior_propI10 +
    log_prior_propI11 +
    log_prior_sigma + log_prior_rho +
    log_prior_inf_sus1 + log_prior_inf_sus4 + log_prior_inf_sus7 + log_prior_inf_sus10 + 
    log_prior_inf_sus2 + log_prior_inf_sus5 + log_prior_inf_sus8 + log_prior_inf_sus11 + 
    log_prior_inf_sus3 + log_prior_inf_sus6 + log_prior_inf_sus9 
  
  # if(!is.finite(log_sum)){browser()}
  
  return(ifelse(log, log_sum, exp(log_sum)))
}

calculate_init_state <- function(pop_sizes, theta_season){
  init_state_season <- c()
  start_sus <- c(dexp(0,theta_season[["inf_sus"]])+1-theta_season[["inf_sus"]],
                 dexp(1,theta_season[["inf_sus"]])+1-theta_season[["inf_sus"]],
                 dexp(2,theta_season[["inf_sus"]])+1-theta_season[["inf_sus"]],
                 dexp(3,theta_season[["inf_sus"]])+1-theta_season[["inf_sus"]],
                 dexp(4,theta_season[["inf_sus"]])+1-theta_season[["inf_sus"]])
  
  
  propR <- theta_season[["propR"]]
  propI <- theta_season[["propI"]]
  
  # for each age group
  for(a in c(1:5)) {
    temp = c((pop_sizes[a] * start_sus[a])-(pop_sizes[a]*propR)- (pop_sizes[a]*propI),
             #SS
             (pop_sizes[a]*propR),0, 0,(pop_sizes[a]*propI),0,0, 0,0,(pop_sizes[a] * (1 - start_sus[a])),
             #IS,PS,RS,SI,II,PI,SP,IP,SR
             0, 0, 0, 0)
    #,RR,R_cases, I_cases, Dual cases
    init_state_season <- c(init_state_season, temp)
    
  }
  return(init_state_season)
}



## Simulate the SIPR model, with detections
simulate_SIPR <- function (theta_season, init_state, times)
{
  outall <- data.table(ode(y = init_state, 
                           t = times, 
                           initfunc = "initmod",
                           dllname = "NT_model_dual",
                           func = "derivatives", 
                           parms = theta_season, 
                           method = "ode23"))
  
  #Calculate the incidence
  trajectory <- summary_stats(outall)
  
  # ][, Hos_INF_0 := INF_incidence_1*theta_season[["Idetect"]]
  # ][, Hos_RSV_1 := RSV_incidence_2*theta_season[["Rdetect5"]]
  # ][, Hos_INF_1 := INF_incidence_2*theta_season[["Idetect"]]*
  #     theta_season[["detect_multiplier"]]]
  
  # trajectory[, Hos_RSV_0 := RSV_incidence_1 + Dual_incidence_1
  # ][, Hos_INF_0 := INF_incidence_1 + Dual_incidence_1
  # ][, Hos_RSV_1 := RSV_incidence_2 + Dual_incidence_2
  # ][, Hos_INF_1 := INF_incidence_2 + Dual_incidence_2]
  
  # 
  return(trajectory)
}


get_attack_rates <- function(outall){
  attack_rates <- matrix(nrow = 5, ncol = 2)
  attack_rates[1,1] <- unlist(outall[nrow(outall),"Rcases1"]/population_numbers[1])
  attack_rates[1,2] <- unlist(outall[nrow(outall),"Icases1"]/population_numbers[1])
  
  attack_rates[2,1] <- unlist(outall[nrow(outall),"Rcases2"]/population_numbers[2])
  attack_rates[2,2] <- unlist(outall[nrow(outall),"Icases2"]/population_numbers[2])
  
  attack_rates[3,1] <- unlist(outall[nrow(outall),"Rcases3"]/population_numbers[3])
  attack_rates[3,2] <- unlist(outall[nrow(outall),"Icases3"]/population_numbers[3])
  
  attack_rates[4,1] <- unlist(outall[nrow(outall),"Rcases4"]/population_numbers[4])
  attack_rates[4,2] <- unlist(outall[nrow(outall),"Icases4"]/population_numbers[4])
  
  attack_rates[5,1] <- unlist(outall[nrow(outall),"Rcases5"]/population_numbers[5])
  attack_rates[5,2] <- unlist(outall[nrow(outall),"Icases5"]/population_numbers[5])
  
  colnames(attack_rates) <- c("RSV", "FLU")
  rownames(attack_rates) <- c("0-1", "2-4", "5-15", "16-64", "65+")
  
  print(attack_rates)
}

## Calculate summary statistics from model output.
summary_stats <- function(outall) {
  # Specify state names
  Names <- c("time")
  for(i in 1:5){
    x<-c(paste0("SS", i), paste0("IS", i),paste0("PS", i),paste0("RS", i),
         paste0("SI", i),paste0("II", i),paste0("PI", i),paste0("SP", i),
         paste0("IP", i), paste0("SR", i),paste0("RR", i),paste0("Rcases", i),
         paste0("Icases", i), paste0("Dual_cases",i))
    Names <- c(Names, x)
  }
  #Calculate the incidence
  colnames(outall) <- Names
  
  outall[, RSV_incidence_1 := (Rcases1 - shift(Rcases1, 1L, type = "lag"))
  ][, INF_incidence_1 := (Icases1 - shift(Icases1, 1L, type = "lag"))
  ][, RSV_incidence_2 := (Rcases2 - shift(Rcases2 ,1L, type = "lag"))
  ][, INF_incidence_2 := (Icases2 - shift(Icases2, 1L, type = "lag"))
  ][, Dual_incidence_1 := (Dual_cases1 - shift(Dual_cases1, 1L, type = "lag"))
  ][, Dual_incidence_2 := (Dual_cases2 - shift(Dual_cases2, 1L, type = "lag"))]
  
  outall[1, c("RSV_incidence_1", "INF_incidence_1",
              "RSV_incidence_2", "INF_incidence_2", 
              "Dual_incidence_1", "Dual_incidence_2")] = 0
  
  
  #print(get_attack_rates(outall))
  
  #remove all but the incidence columsn
  outall[, setdiff(names(outall), c("RSV_incidence_1","INF_incidence_1",
                                    "RSV_incidence_2", "INF_incidence_2",
                                    "Dual_incidence_1", "Dual_incidence_2")) := NULL][] 
  #Convert to weekly
  outall$week <- c(sapply(1:(nrow(outall)/7), function(x){return(rep(x,7))}))
  outall_week <- outall[, lapply(.SD, sum), by=list(week)]
  
  
  return(as.data.table(outall_week))
}
# temp_func <- function(Temperature, i){
#   if(i == 1){beta_i <- 1} else{
#     beta_i <- 1 / ( (1 + Temperature) * (i - 1) )
#   }
#   return(beta_i)
# }

temp_func <- function(Temperatures, i){
  beta_i <- 1/Temperatures[i]
  return(beta_i)
}



step_fn <- function(parameters, covmat,lower_bounds, upper_bounds, theta_names){
  #print(theta_names)
  parameters_new <- c(rtmvnorm(1,
                               mean = unlist(parameters[theta_names], use.names=F),
                               sigma = covmat[theta_names,
                                              theta_names],
                               lower = lower_bounds,
                               upper = upper_bounds
  ))
  
  names(parameters_new) <- theta_names
  
  for(param_name in names(parameters_new)){
    parameters[param_name] <- parameters_new[param_name]}
  
  return(parameters)
}



adjust_unsymmetrical <- function(log_acceptance, parameters_old,
                                 parameters_new, covmat, lower_bounds, upper_bounds,
                                 theta_names){
  
  
  log_acceptance <- log_acceptance + dtmvnorm(x = unlist(parameters_old[theta_names]),
                                              mean = unlist(parameters_new[theta_names]),
                                              sigma = covmat[theta_names,
                                                             theta_names],
                                              lower = lower_bounds[theta_names],
                                              upper = upper_bounds[theta_names],
                                              log = TRUE)
  
  log_acceptance <- log_acceptance - dtmvnorm(x = unlist(parameters_new[theta_names]),
                                              mean = unlist(parameters_old[theta_names]),
                                              sigma = covmat[theta_names,
                                                             theta_names],
                                              lower = lower_bounds[theta_names],
                                              upper = upper_bounds[theta_names],
                                              log = TRUE)
  
  return(log_acceptance)
}


swapping_rates <- function(trace_in, n_pars, prev_swap_proposed, prev_swap_accepted){
  
  swap_tot <- sum(trace_in[,n_pars+4], na.rm = T) # total swaps proposed
  swap_tot <- swap_tot + prev_swap_proposed
  swap_accep <- length(which(trace_in[,n_pars+3]!= 0)) # total swaps accepted
  swap_accep <- swap_accep + prev_swap_accepted
  
  if(swap_tot ==0){
    swap_rate <- 0} else {
      swap_rate <- swap_accep/swap_tot
      
    }
  
  return(swap_rate)
}

swapping_props <- function(trace_in, n_pars, prev_swap_proposed){
  
  swap_tot <- sum(trace_in[,n_pars+4], na.rm = T) # total swaps proposed
  swap_tot <- swap_tot + prev_swap_proposed
  
  return(swap_tot)
}

swapping_acceps <- function(trace_in, n_pars, prev_swap_accepted){
  
  swap_accep <- length(which(trace_in[,n_pars+3]!= 0)) # total swaps accepted
  swap_accep <- swap_accep + prev_swap_accepted
  
  return(swap_accep)
}


add_on <- function(tt, t, pt, each_run){
  tt[(((pt-1)*each_run)+1):(((pt-1)*each_run)+each_run),] <- t
  return(tt)
  
}

trace_wrapper <- function(length_run ,
                          each_run ,
                          n_chains,
                          init_theta,
                          data_to_fit,
                          population_numbers,
                          indep,
                          covmat,
                          adapt_rate,
                          adaption_starter,
                          proposal_sd,
                          Temperatures,
                          name,
                          lower_bounds = lower_bounds,
                          upper_bounds = upper_bounds
){
  
  
  #calculate in which pt to start the adaption
  # which one should it start in and the remainder
  pt_start <- ceiling(adaption_starter/each_run)
  pt_rem <- adaption_starter %% each_run
  
  total_trace <- lapply(1:(n_chains), function(i)  matrix(NA, nrow  = length_run,
                                                          ncol = (4 + length(init_theta))) )
  temperature_store <- matrix(nrow = length_run/each_run, 
                              ncol = length(Temperatures))
  # wrapper around traces
  for (pt in 1:(length_run/each_run)) {
    
    # specify the adapt start for each run.
    if(pt == pt_start){
      if( pt_rem == 0 ) { adapt_start <- each_run } else{
        adapt_start <- pt_rem
      }
    } else if (pt < pt_start) {adapt_start <- Inf
    } else if (pt > pt_start) {adapt_start <- 0}
    # specify the init_theta as the main one if the first one
    if( pt == 1){
      init_theta <- lapply(1:n_chains, function(x) as.list(init_theta))
      prev_swap_proposed <- rep(0, n_chains)
      prev_swap_accepted <- rep(0, n_chains)
    }
    # Run the mcmc for the pt section
    
    trace <- mcmcmc_fn(init_theta = init_theta,
                       data_to_fit = data_to_fit,
                       pop_size = population_numbers,
                       MaxTime = each_run,
                       indep = indep,
                       covmat = covmat,
                       n_chains = n_chains,
                       adapt_rate = adapt_rate,
                       adapt_start = adapt_start,
                       prev_swap_proposed = prev_swap_proposed,
                       prev_swap_accepted = prev_swap_accepted,
                       pt = pt,
                       Temperatures = Temperatures,
                       lower_bounds = lower_bounds,
                       upper_bounds = upper_bounds)
    # save the trace into total_trace
    total_trace <- lapply(1:(length(total_trace)), function(x)
      add_on(total_trace[[x]], trace[[x]], pt, each_run))
    # update the init_theta to the new one
    end_trace <- lapply(1:n_chains, function (x)
      as.list(tail(na.omit(total_trace[[x]]),1)[,2:(length(proposal_sd)+3)]))
    for ( section in 1:n_chains){
      names(end_trace[[section]]) <- names(proposal_sd)}
    init_theta <- end_trace
    # update the temperatures to the correct ones
    Temperatures <- c(tail(trace[["temperatures"]],1))
    prev_swap_proposed <- trace[["swaps_proposed"]]
    prev_swap_accepted <- trace[["swaps_accepted"]]
    
    temperature_store[pt,] <- Temperatures
    
    print(paste0("Overall: ", pt/(length_run/each_run)*100, "%"))
    save(total_trace, file=paste0(name, "_", array_num, "_", Sys.Date(),
                                  ".Rdata"))
    #print(total_trace)
  }
  total_trace$temperatures <- temperature_store
  return(total_trace)
}


subset_parameters <- function(parameters, season){
  parameters_season <- c(
    betaI = unname(parameters[season_params[season,"betaI"]]),
    Idetect = unname(parameters[season_params[season,"Idetect"]]),
    propR = unname(parameters[ season_params[season, "propR"]]),
    inf_sus = unname(parameters[season_params[season, "inf_sus"]]),
    Rdetect2 = unname(parameters[season_params[season, "Rdetect2"]]),
    Rdetect5 = unname(parameters[season_params[season, "Rdetect5"]]),
    sig = unname(parameters[season_params[season, "sig"]]),
    rho = unname(parameters[season_params[season, "rho"]]),
    betaR = unname(parameters[season_params[season, "betaR"]]),
    propI = unname(parameters[season_params[season, "propI"]]),
    detect_multiplier = unname(parameters[season_params[season, "multiplier"]]),
    RSV_mult = unname(parameters[season_params[season, "RSV_detect_multiplier"]]), 
    overdispersion = unname(parameters["overdispersion"]), 
    Dual_mult = unname(parameters["Dual_mult"])
  )
  
  return(parameters_season)
}



trace_wrapper_cont <- function(trace_previous,
                               additional_run = 30000,
                               old_length_run = 30000,
                               each_run = 5000,
                               n_chains = n_chains,
                               data_to_fit = multi_data,
                               population_numbers = population_numbers,
                               indep = 1,
                               covmat = covmat,
                               adapt_rate = 0.5,
                               adaption_starter = 10000,
                               proposal_sd = proposal_sd,
                               Temperatures = Temperatures,
                               name = name,
                               save = T
){
  
  total_trace <- lapply(1:(n_chains), function(i)  matrix(NA, nrow  = additional_run,
                                                          ncol = (4 + length(init_theta))) )
  temperature_store <- matrix(nrow = additional_run/each_run, 
                              ncol = length(Temperatures))
  # check adatption was on
  if (adaption_starter < old_length_run){
    adapt_start = 0
  } else { adapt_start = Inf}
  
  # update the init_theta from  the old trace
  end_trace <- lapply(1:n_chains, function (x)
    as.list(tail(trace_previous[[x]],1)[,2:(length(lower_bounds)+1)]))
  for (section in 1:n_chains){
    names(end_trace[[section]]) <- names(lower_bounds)}
  init_theta <- end_trace
  # update the temperatures from the old trace
  Temperatures <- c(tail(trace_previous[[length(trace_previous)]],1))
  prev_swap_proposed <- sapply(1:n_chains, function(x)
    sum(trace_previous[[x]][,length(lower_bounds)+3]))
  prev_swap_accepted <-  sapply(1:n_chains, function(x)
    sum(trace_previous[[x]][,length(lower_bounds)+2]))
  
  # wrapper around traces
  for (pt in 1:(additional_run/each_run)) {
    
    trace <- mcmcmc_fn(init_theta = init_theta,
                       data_to_fit = data_to_fit,
                       pop_size = population_numbers,
                       MaxTime = each_run,
                       indep = indep,
                       covmat = covmat,
                       n_chains = n_chains,
                       adapt_rate = adapt_rate,
                       adapt_start = adapt_start,
                       prev_swap_proposed = prev_swap_proposed,
                       prev_swap_accepted = prev_swap_accepted,
                       pt = pt,
                       Temperatures = Temperatures,
                       lower_bounds = lower_bounds,
                       upper_bounds = upper_bounds,
                       previous_trace = old_length_run)
    # save the trace into total_trace
    total_trace <- lapply(1:(length(total_trace)), function(x)
      add_on(total_trace[[x]], trace[[x]], pt, each_run))
    
    # update the init_theta to the new one
    end_trace <- lapply(1:n_chains, function(x)
      as.list(tail(na.omit(total_trace[[x]]),1)[,2:(length(lower_bounds)+1)]))
    for ( section in 1:n_chains){
      names(end_trace[[section]]) <- names(lower_bounds)}
    init_theta <- end_trace
    
    # update the temperatures from the old trace
    Temperatures <- c(tail(trace[["temperatures"]],1))
    prev_swap_proposed <- trace[["swaps_proposed"]]
    prev_swap_accepted <- trace[["swaps_accepted"]]
    
    temperature_store[pt,] <- Temperatures
    
    print(paste0("Overall: ", pt/(additional_run/each_run)*100, "%"))
    save(total_trace, file=paste0(name, "_", array_num, "_", Sys.Date(),
                                  ".Rdata"))
    #  browser()
  }
  total_trace$temperatures <- temperature_store
  return(total_trace)
}

