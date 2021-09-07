# FIgure generation functions
calc_start_sus <- function(value){
start_sus <- c(dexp(0,value)+1-value,
               dexp(1,value)+1-value,
               dexp(2,value)+1-value,
               dexp(3,value)+1-value,
               dexp(4,value)+1-value)
}


get_attack_rates <- function(outall){
  attack_rates <- matrix(nrow = 5, ncol = 2)
  attack_rates[1,1] <- unlist((outall[nrow(outall),"Rcases1"]+
                                outall[nrow(outall),"Dual_cases1"]) /population_numbers[1])
  attack_rates[1,2] <- unlist((outall[nrow(outall),"Icases1"] +
                                outall[nrow(outall),"Dual_cases1"])/population_numbers[1])
  
  attack_rates[2,1] <- unlist((outall[nrow(outall),"Rcases2"] +
                                outall[nrow(outall),"Dual_cases2"])/population_numbers[2])
  attack_rates[2,2] <- unlist((outall[nrow(outall),"Icases2"]+
                                outall[nrow(outall),"Dual_cases2"])/population_numbers[2])
  
  attack_rates[3,1] <- unlist((outall[nrow(outall),"Rcases3"] +
                                outall[nrow(outall),"Dual_cases3"])/population_numbers[3])
  attack_rates[3,2] <- unlist((outall[nrow(outall),"Icases3"]+
                                outall[nrow(outall),"Dual_cases3"])/population_numbers[3])
  
  attack_rates[4,1] <- unlist((outall[nrow(outall),"Rcases4"]+
                                outall[nrow(outall),"Dual_cases4"])/population_numbers[4])
  attack_rates[4,2] <- unlist((outall[nrow(outall),"Icases4"]+
                                outall[nrow(outall),"Dual_cases4"])/population_numbers[4])
  
  attack_rates[5,1] <- unlist((outall[nrow(outall),"Rcases5"]+
                                outall[nrow(outall),"Dual_cases5"])/population_numbers[5])
  attack_rates[5,2] <- unlist((outall[nrow(outall),"Icases5"]+
                                outall[nrow(outall),"Dual_cases5"])/population_numbers[5])
  
  colnames(attack_rates) <- c("RSV", "FLU")
  rownames(attack_rates) <- c("0-1", "2-4", "5-15", "16-64", "65+")
  
  return(attack_rates)
}

run_model_attack <- function(theta_season, init_state, season_data){
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

# theta_season <- pars_transformation_trans_int(theta_season)
times <- seq(from = 1, to = nrow(season_data)*7, by = 1)
# run the mod
outall <- data.table(ode(y = init_state, 
                         t = times, 
                         initfunc = "initmod",
                         dllname = "NT_model_dual",
                         func = "derivatives", 
                         parms = theta_season, 
                         method = "ode23"))

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

attack_rates <- get_attack_rates(outall)

case_totals <- get_cases_season(outall)
case_totals[1,1] <- case_totals[1,1] * theta_season[["Rdetect2"]]*theta_season[["RSV_mult"]]
case_totals[2,1] <- case_totals[2,1] * theta_season[["Rdetect5"]]*theta_season[["RSV_mult"]]

case_totals[1,2] <- case_totals[1,2] * theta_season[["Idetect"]]
case_totals[2,2] <- case_totals[2,2] * theta_season[["Idetect"]]*theta_season[["detect_multiplier"]]

case_totals[1,3] <- case_totals[1,3] * theta_season[["Rdetect2"]]*theta_season[["RSV_mult"]]
case_totals[2,3] <- case_totals[2,3] * theta_season[["Rdetect5"]]*theta_season[["RSV_mult"]]
return(list(Attack = attack_rates, 
               Cases = case_totals))
}


#### difference in attack rates with vaccines #####

attack_rates_vaccination <- function(total_trace, fixed_sim, sample_size,
                                together = F,
                                rolling_mean = 1, from, to){
  for(season in 1:11){
    # Need to change numbers if number of parameters changes.
    if (together == T){
      trace_burned <- as.data.frame(total_trace[c(from:to),])
    }else {
      trace_burned <- as.data.frame(total_trace[ c(from:to),2:11])
    }
    # select a sample of 1000 parameter combinations
    sample_trace <- Momocs::sample_n(tbl= trace_burned, size = sample_size, replace = T)
    # for each combination
    Hos_RSV_0_store <- data.frame(NULL)
    Hos_RSV_1_store <- data.frame(NULL)
    noflu_RSV_0_store <- data.frame(NULL)
    noflu_RSV_1_store <- data.frame(NULL)
    Hos_INF_0_store <- data.frame(NULL)
    Hos_INF_1_store <- data.frame(NULL)
    norsv_INF_0_store <- data.frame(NULL)
    norsv_INF_1_store <- data.frame(NULL)
    
    
    for (r in 1:nrow(sample_trace)){
      
      # run the model
      # not sure that this will work as the input probably a different format
      
      sample_trace$RSV_default <- 1
      pars_out <-unlist(subset_parameters(sample_trace[r,], season))
      pars_df <- matrix(pars_out,nrow=1)
      colnames(pars_df) <- names(pars_out)
      pars_df[,"sig"] = 1-pars_df[,"sig"] 
      attack_rates <- s(pars_df,season)
      
      pars_out1 <- pars_out
      pars_out1["propI"] <- 0
      pars_df <- matrix(pars_out1,nrow=1)
      colnames(pars_df) <- names(pars_out)
      sim_out_noflu <- run_simulation_quantiles(pars_df,season)
      pars_out2 <- pars_out
      pars_out2["propR"] <- 0
      pars_df <- matrix(pars_out2,nrow=1)
      colnames(pars_df) <- names(pars_out)
      sim_out_norsv <- run_simulation_quantiles(pars_df,season)
      # save each one into dataframe with other samples
      #each sample saved as row, so columns are timestep
      Hos_RSV_0_store <- rbind(Hos_RSV_0_store,
                               unname(unlist(sim_out[1:63,"Hos_RSV_0"])))
      Hos_RSV_1_store <- rbind(Hos_RSV_1_store,
                               unname(unlist(sim_out[1:63,"Hos_RSV_1"])))
      noflu_RSV_0_store <- rbind(noflu_RSV_0_store,
                                 unname(unlist(sim_out_noflu[1:63,"Hos_RSV_0"])))
      noflu_RSV_1_store <- rbind(noflu_RSV_1_store,
                                 unname(unlist(sim_out_noflu[1:63,"Hos_RSV_1"])))
      Hos_INF_0_store <- rbind(Hos_INF_0_store,
                               unname(unlist(sim_out[1:63,"Hos_INF_0"])))
      Hos_INF_1_store <- rbind(Hos_INF_1_store,
                               unname(unlist(sim_out[1:63,"Hos_INF_1"])))
      norsv_INF_0_store <- rbind(norsv_INF_0_store,
                                 unname(unlist(sim_out_norsv[1:63,"Hos_INF_0"])))
      norsv_INF_1_store <- rbind(norsv_INF_1_store,
                                 unname(unlist(sim_out_norsv[1:63,"Hos_INF_1"])))
      
      colnames(Hos_RSV_0_store) <- seq(1:63)
      colnames(Hos_RSV_1_store) <- seq(1:63)
      colnames(noflu_RSV_0_store) <- seq(1:63)
      colnames(noflu_RSV_1_store) <- seq(1:63)
      colnames(Hos_INF_0_store) <- seq(1:63)
      colnames(Hos_INF_1_store) <- seq(1:63)
      colnames(norsv_INF_0_store) <- seq(1:63)
      colnames(norsv_INF_1_store) <- seq(1:63)
    }
    # take the 95% quantiles for each group and format back into one table
    temp_quantiles_RSV_0 <- as.data.table(apply(Hos_RSV_0_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    temp_quantiles_RSV_0[,Quantile := c("p0.025","p0.5","p0.975")]
    RSV_0 <- data.table(melt(temp_quantiles_RSV_0, id="Quantile"))
    RSV_0 <- dcast.data.table(RSV_0, variable~Quantile, value.var = "value")
    RSV_0[,people := "A) RSV < 2"]
    RSV_0[,type := "Current"]
    
    temp_quantiles_RSV_1 <- as.data.table(apply(Hos_RSV_1_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    temp_quantiles_RSV_1[,Quantile := c("p0.025","p0.5","p0.975")]
    RSV_1 <- data.table(melt(temp_quantiles_RSV_1, id="Quantile"))
    RSV_1 <- dcast.data.table(RSV_1, variable~Quantile, value.var = "value")
    RSV_1[,people := "B) RSV 2-5"]
    RSV_1[,type := "Current"]
    
    noflu_quantiles_RSV_0 <- as.data.table(apply(noflu_RSV_0_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    noflu_quantiles_RSV_0[,Quantile := c("p0.025","p0.5","p0.975")]
    RSV_noflu_0 <- data.table(melt(noflu_quantiles_RSV_0, id="Quantile"))
    RSV_noflu_0 <- dcast.data.table(RSV_noflu_0, variable~Quantile, value.var = "value")
    RSV_noflu_0[,people := "A) RSV < 2"]
    RSV_noflu_0[,type := "Vaccination"]
    
    noflu_quantiles_RSV_1 <- as.data.table(apply(noflu_RSV_1_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    noflu_quantiles_RSV_1[,Quantile := c("p0.025","p0.5","p0.975")]
    RSV_noflu_1 <- data.table(melt(noflu_quantiles_RSV_1, id="Quantile"))
    RSV_noflu_1 <- dcast.data.table(RSV_noflu_1, variable~Quantile, value.var = "value")
    RSV_noflu_1[,people := "B) RSV 2-5"]
    RSV_noflu_1[,type := "Vaccination"]
    
    temp_quantiles_INF_0 <- as.data.table(apply(Hos_INF_0_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    temp_quantiles_INF_0[,Quantile := c("p0.025","p0.5","p0.975")]
    INF_0 <- data.table(melt(temp_quantiles_INF_0, id="Quantile"))
    INF_0 <- dcast.data.table(INF_0, variable~Quantile, value.var = "value")
    INF_0[,people := "C) INF < 2"]
    INF_0[,type := "Current"]
    
    temp_quantiles_INF_1 <- as.data.table(apply(Hos_INF_1_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    temp_quantiles_INF_1[,Quantile := c("p0.025","p0.5","p0.975")]
    INF_1 <- data.table(melt(temp_quantiles_INF_1, id="Quantile"))
    INF_1 <- dcast.data.table(INF_1, variable~Quantile, value.var = "value")
    INF_1[,people := "D) INF 2-5"]
    INF_1[,type := "Current"]
    
    norsv_quantiles_INF_0 <- as.data.table(apply(norsv_INF_0_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    norsv_quantiles_INF_0[,Quantile := c("p0.025","p0.5","p0.975")]
    INF_norsv_0 <- data.table(melt(norsv_quantiles_INF_0, id="Quantile"))
    INF_norsv_0 <- dcast.data.table(INF_norsv_0, variable~Quantile, value.var = "value")
    INF_norsv_0[,people := "C) INF < 2"]
    INF_norsv_0[,type := "Vaccination"]
    
    norsv_quantiles_INF_1 <- as.data.table(apply(norsv_INF_1_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    norsv_quantiles_INF_1[,Quantile := c("p0.025","p0.5","p0.975")]
    INF_norsv_1 <- data.table(melt(norsv_quantiles_INF_1, id="Quantile"))
    INF_norsv_1 <- dcast.data.table(INF_norsv_1, variable~Quantile, value.var = "value")
    INF_norsv_1[,people := "D) INF 2-5"]
    INF_norsv_1[,type := "Vaccination"]
    
    quantile_table <- rbind(RSV_0,RSV_1, RSV_noflu_0, RSV_noflu_1, 
                            INF_0,INF_1, INF_norsv_0, INF_norsv_1)
    quantile_table$season <- season
    quantile_table <- na.omit(quantile_table)
    quantile_table$week_begin <- rep(fixed_sim[[season]]$start_week, 8)
    
    if(season ==1){
      all_seasons <- quantile_table
    } else {
      all_seasons <- rbind(all_seasons, quantile_table)
    }
  }
  all_seasons$variable <- as.numeric(as.character(all_seasons$variable))
  all_seasons$season <- as.numeric(as.character(all_seasons$season))
  all_seasons$week_begin <- as.Date(all_seasons$week_begin)
  # plot graph]
  SIM_PLOT <- ggplot(all_seasons) +
    geom_line(aes(x=week_begin, y=p0.5, colour = type))+
    geom_ribbon(aes(x=week_begin,ymin=p0.025, ymax = p0.975, 
                    fill = type), alpha = 0.5) +
    facet_wrap(people~., ncol = 1, scales = "free")+
    labs(y="Incidence",
         colour = "Group", fill="Group") +
    theme_bw() +
    # guides(shape = guide_legend(override.aes = list(size = 0.5)),
    #        color = guide_legend(override.aes = list(size = 0.5)))+
    # legend_key_width= unit(0.05,"cm"))+
    theme(legend.title = element_text(size = 10),
          legend.text  = element_text(size = 10),
          axis.text = element_text(size=10),
          axis.title = element_text(size=10, margin = margin(l=0,r=100,b=0,t=0)),
          legend.key.width = unit(0.25,"cm"),
          axis.title.x =  element_blank(),
          # legend.margin = margin(c(1,1,1,1))
          legend.key.height = unit(0.3,"cm"),
          title = element_text(size=10),
          strip.text.y = element_text(angle = 0), 
          strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(colour = 'black',hjust = 0))+
    geom_vline(xintercept = dates_to_cut, alpha= 0.5)
  return(SIM_PLOT)
  
}


get_cases_season <- function(outall){
  tot_cases <- matrix(nrow = 2, ncol = 3)
  tot_cases[1,1] <- unlist(outall[nrow(outall),"Rcases1"])
  tot_cases[1,2] <- unlist(outall[nrow(outall),"Icases1"])
  tot_cases[1,3] <- unlist(outall[nrow(outall),"Dual_cases1"])
  
  tot_cases[2,1] <- unlist(outall[nrow(outall),"Rcases2"])
  tot_cases[2,2] <- unlist(outall[nrow(outall),"Icases2"])
  tot_cases[2,3] <- unlist(outall[nrow(outall),"Dual_cases2"])

  colnames(tot_cases) <- c("RSV", "FLU", "Dual")
  rownames(tot_cases) <- c("0-1", "2-4")
  
  return(tot_cases)
}
