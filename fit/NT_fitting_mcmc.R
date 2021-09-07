# Metropolis Coupled Markov Chain Monte Carlo

 # To use for transformed parameters (i.e. adjusting for unsymnetical distributions)

# Args:
#   init_theta: either a vector of intial paramter values
#               or a list of initial parameter values for each chaing
#   data_to_fit: the data!
#   pop_size: population sizes by age group (vector)
#   MaxTime: length of time to run the chain
#   indep: period of time for which chains should wander independently
#   covmat: covariance matrix for the parameter proposals
#   n_chains: number of chains to run. Must equal length of temperature vector
#   adapt_rate: at what rate to adapt the temperatures. Default 1.
#   adapt_start: at what point to start adapting. Default infinite


# Returns:
#   chains: list containing matrix for each chain, first col is loglik + log prior prob,
#           remaining columns are fn parameters in order given in the pars[[i]]
#   Temperatures: list of temperature values every indep step

mcmcmc_fn <- function(init_theta, data_to_fit, pop_size, MaxTime = 1e3,
                      indep = 100, covmat, n_chains,
                      adapt_rate=1, adapt_start = Inf,
                      prev_swap_proposed, prev_swap_accepted, pt, Temperatures,
                      previous_trace = 0,
                      ...){

  theta_names <- names(lower_bounds)

  Temperatures_store <- matrix(nrow=MaxTime/indep, ncol=n_chains)
  adapt_switch  <- T
  #  covmat <- covmat + diag(ncol(covmat))*0.01
  pars <- init_theta
  n_pars <- length(theta_names)
  # list to store likelihood of each chain
  if(n_chains == 1){
    likelihood <-data.frame()
  } else{
    likelihood <- list(length = n_chains) }
  # create the empty traces for each chain
  chains <- lapply(1:n_chains, function(i)  matrix(NA, nrow  = MaxTime,
                                                   ncol = (4 + n_pars)) )
  # set the first row of each trace to 0
  for (i in n_chains){ chains[[i]][,n_pars + 3] <- 0 }

  # The independent intervals, run model in parallel for each chain
  Interval <- matrix(1:MaxTime, nrow = indep)

  print(Sys.time())

  # the outer time loop over intervals
  for(s in 1:(MaxTime / indep)){
    # print(paste0("s is ", s))
    # Evolve chains independently for "indep" time steps

    out <- lapply( 1:n_chains,
                   function(i){
                     # matrix to store output
                     out <- matrix(NA, ncol = length(theta_names) + 4,
                                   nrow = indep)
                     # For each time step in indep
                     #    print(paste0("i is ", i))
                     for(t in 1:indep){
                       #  print(paste0("t is ", t))

                       #if first step, calculate the likelhiood of the parameters
                       if (t == 1){
                         #transform the parameters
                         pars_transformed <- pars_transformation(pars[[i]])
                         #for each season
                         season_out <-  sum(sapply(1:length(data_to_fit), function(season){
                           #subset parameters
                           pars_season <- subset_parameters(pars_transformed, season)
                           data_season <- data_to_fit[[season]]
                           #calculate initial state
                           init_state_season <- calculate_init_state(pop_size, pars_season)
                           #print(init_state_season
                       #    print(pars_transformed)
                           
                           #calculate likelihood of paraeters
                           liklihood_init <- calc_likelihood_overall(pars_season,
                                                                     init_state_season,
                                                                     data_season,
                                                                     log = T)
                           #   if(!is.finite(liklihood_init)){browser()}
                           # if(is.na(liklihood_init)){browser()}
                           return(liklihood_init)

                         }))

                         #add together likelihood and priors
                         Pi <- season_out + calc_priors(pars_transformed,
                                                        R0ratI,
                                                        R0ratR,
                                                        log = T)
                         #else for all other steps extract previous likelihood from out.
                       } else {
                         Pi <- out[t-1,1]
                       }

                       #propose new parameters
                     #  covmat <- as.matrix(nearPD(covmat)[[1]])
                       proposed <- step_fn(pars[[i]], covmat,
                                           lower_bounds, upper_bounds,
                                           theta_names)

                       # prior of proposed parameters
                       proposed_transformed <- pars_transformation(proposed)
                       #   print(unlist(proposed_transformed))
                       Prior_proposed <- calc_priors(proposed_transformed,
                                                     R0ratI,
                                                     R0ratR,
                                                     log = T)
                       # if(!is.finite(Prior_proposed)){browser()}

                       #if priors are infinite, rejct step
                       if (is.finite(Prior_proposed)){

                         # calculate the temperature weighting for the chain
                         beta_to_use <- temp_func(Temperatures, i)

                         season_out_proposed <-  sum(sapply(1:length(data_to_fit), function(season){
                           #subset parameters
                           pars_season <- subset_parameters(proposed_transformed, season)
                           data_season <- data_to_fit[[season]]

                           init_state_proposed <- calculate_init_state(pop_size,
                                                                       pars_season)
                           likelihood_proposed <- calc_likelihood_overall(pars_season,
                                                                          init_state_proposed,
                                                                          data_season,
                                                                          log = T)
                           #  if(!is.finite(likelihood_proposed)){browser()}
                           return(likelihood_proposed)
                         }))
                         # Calculate acceptance ratio
                         log_acceptance <- season_out_proposed + Prior_proposed - Pi
                         #    if(!is.finite(log_acceptance)){browser()}
                         #adjst the acceptance for the unsymmetrical proposal distribution

                         ## DO not need to adjust if using symmetric proposals.

                         log_acceptance <- adjust_unsymmetrical(log_acceptance = log_acceptance,
                                                                  parameters_old = pars[[i]],
                                                                  parameters_new = proposed,
                                                                  covmat = covmat,
                                                                  lower_bounds =lower_bounds,
                                                                upper_bounds = upper_bounds,
                                                                  theta_names)

                         #alpha is the weighted liklihood of the proposal by the temperature
                         alpha <- exp(beta_to_use * log_acceptance )

                         #check alpha is a number, if not make -inf
                         if(is.nan(alpha)){alpha <- -Inf}
                         # if new parameters accepted
                         if ( alpha  > runif(1) ){
                           #update current parameters and likelihood
                           pars[[i]] <- proposed
                           Liklihood_stored <- (season_out_proposed + Prior_proposed)
                           #   if(!is.finite(season_out_proposed)){browser()}
                           #  if(!is.finite(Prior_proposed)){browser()}
                           Accepted <- 1
                           #else store the likelihoood of the previous parameters
                         }else {

                           Liklihood_stored <- Pi
                           Accepted <- 0}} else{
                             #if the prior was infinite
                             Liklihood_stored <- Pi
                             Accepted <- 0
                           #  print("rejected at prior")
                             }
                       #update the trace
                       out[t,] <- c(Liklihood_stored, unlist(pars[[i]][theta_names]),
                                    Accepted, 0,0)
                     }
                     # print(paste0("t is ", t))
                     out
                   })
    #print(paste0("s is ", s))
    # save the output from each indep into the main chains
    for(i in 1:n_chains){

      #copy the trace for the interval into the overall trace (chain)
      chains[[i]][Interval[,s],] <- out[[i]]

      # copy parameters from last step of indep trace
      pars[[i]] <- as.list(out[[i]][indep,][-c(1,n_pars+2, n_pars+3, n_pars+4)])
      names(pars[[i]]) <- theta_names
      # store the likelihood value for each so don't have to recaluclate later

      if(n_chains ==1){
        likelihood <- out[[i]][indep,][1]
      } else{
        likelihood[[i]]<- out[[i]][indep,][1]
      }
    }
    # Propose a chain swaps every "indep" time steps. num swaps - num chains
    if(n_chains > 1){
      for (swap_test in 1:n_chains){
        # Choose two out of total number of chains and store their numbers in 'pick'
        pick <- sample(1:(n_chains-1), 1)
        i <- pick; j <- pick+1 # for convience
        # calculate the swap potential
        #TODO is this the correct acceptance?
        if(is.finite(likelihood[[i]]) & is.finite(likelihood[[j]])){
          R <- exp((likelihood[[i]] - likelihood[[j]]) *
                     ( temp_func(Temperatures , j) - temp_func(Temperatures , i)))
        } else { R = -Inf}

        #store swap attempted
        chains[[i]][Interval[dim(Interval)[1],s], n_pars + 4] <- 1

        # accept or reject swap
        if(R > runif(1)){
          #if accepted, change parameters i to parameters j.
          pars_temp <- pars[[i]]
          pars[[i]] <- pars[[j]]
          pars[[j]] <- pars_temp
          #store the number it swapped to
          chains[[i]][Interval[dim(Interval)[1],s],n_pars+3]<-j

        } } }

    if(Interval[1,s] > adapt_start){
      if(adapt_switch == T){
        print("Starting to adapt at this point.......")
        adapt_switch <- F
      }
      # calculate swap rate for each temperatre
      # print(prev_swap_proposed)
      # print(prev_swap_accepted)
      swap_rates <- sapply(1:n_chains, function(i)
        swapping_rates(chains[[i]], n_pars = n_pars,
                       prev_swap_proposed = prev_swap_proposed[i],
                       prev_swap_accepted = prev_swap_accepted[i]))
      # swap_rates[1] is between 1 and 2, swap_rates[2] is between 2 and 3 etc.
      # calculate kappa
      if(previous_trace==0){
        time_adapting<- Interval[dim(Interval)[1],s] + (pt-1)*MaxTime - adapt_start
      } else{
        time_adapting <- (Interval[dim(Interval)[1],s] +
                            (pt-1)*MaxTime - adapt_start) +
          (previous_trace - adapt_start)
      }

      kappa <- 1/(1+time_adapting)^(adapt_rate)
      # store current temperatures
      temps_store <- Temperatures
      for(tem in 2:(n_chains-1)){

        # calculate s = log(Ti - Ti+1) for first chain
        S2 <- log(temps_store[tem]- temps_store[tem-1])
        # calculate change in Schange = k[Ai - Ai+1] + S for first chain
        s2new <- kappa*(swap_rates[tem-1] - swap_rates[tem])+S2
        #calculate new temperature, based on new temperature for previous one
        #TODO check this is the right way round?
        Tnew <- exp(s2new)+Temperatures[tem-1]
        Temperatures[tem] <- Tnew
      }}
    #print at right point in time
    # percent_points <- seq(from = 0, to = MaxTime / indep, by = (MaxTime / indep) /10 )
    # if(s %in% percent_points){
    #   print(paste0("pt complete: ", s/(MaxTime/indep)*100, "% - ",
    #                " time is ", Sys.time()))}
    Temperatures_store[s,] <- Temperatures

  }
  # Returns the full history of all chains
  swaps_proposed <- sapply(1:n_chains, function(i) swapping_props(chains[[i]], n_pars=n_pars,
                                                                  prev_swap_proposed = prev_swap_proposed[i]))
  swaps_accepted <- sapply(1:n_chains, function(i) swapping_acceps(chains[[i]], n_pars=n_pars,
                                                                   prev_swap_accepted = prev_swap_accepted[i]))
  #TODO return the swap attempts and swaps accepted for each chain
  chains[["temperatures"]] <- Temperatures_store
  chains[["swaps_proposed"]] <- swaps_proposed
  chains[["swaps_accepted"]] <- swaps_accepted
  return(chains)
}
