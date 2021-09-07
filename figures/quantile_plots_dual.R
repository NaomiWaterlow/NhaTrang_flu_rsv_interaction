# quantiles plots dual 

# plot the quantiles of the fit: functions needed


# take samples from the trace and generate quantiles + plot
fit_quantiles_together <- function(total_trace, fixed_sim, sample_size,
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
    sample_trace[,"RSV_default"] <- 1
    # for each combination
    Hos_RSV_0_store <- data.frame(NULL)
    Hos_RSV_1_store <- data.frame(NULL)
    Hos_INF_0_store <- data.frame(NULL)
    Hos_INF_1_store <- data.frame(NULL)
    Dual_0_store <- data.frame(NULL)
    Dual_1_store <- data.frame(NULL)
    
    for (r in 1:nrow(sample_trace)){
      # run the model
      # not sure that this will work as the input probably a different format
      pars_out <-unlist(subset_parameters(sample_trace[r,], season))
      pars_df <- matrix(pars_out,nrow=1)
      colnames(pars_df) <- names(pars_out)
      pars_df[,"sig"] = 1-pars_df[,"sig"] 

      sim_out <- run_simulation_quantiles(pars_df,season)
      # save each one into dataframe with other samples
      #each sample saved as row, so columns are timestep
      Hos_RSV_0_store <- rbind(Hos_RSV_0_store,
                               unname(unlist(sim_out[1:66,"Hos_RSV_0"])))
      Hos_RSV_1_store <- rbind(Hos_RSV_1_store,
                               unname(unlist(sim_out[1:66,"Hos_RSV_1"])))
      Hos_INF_0_store <- rbind(Hos_INF_0_store,
                               unname(unlist(sim_out[1:66,"Hos_INF_0"])))
      Hos_INF_1_store <- rbind(Hos_INF_1_store,
                               unname(unlist(sim_out[1:66,"Hos_INF_1"])))
      Dual_0_store <- rbind(Dual_0_store,
                               unname(unlist(sim_out[1:66,"Dual_0"])))
      Dual_1_store <- rbind(Dual_1_store,
                               unname(unlist(sim_out[1:66,"Dual_1"])))
      colnames(Hos_RSV_0_store) <- seq(1:66)
      colnames(Hos_RSV_1_store) <- seq(1:66)
      colnames(Hos_INF_0_store) <- seq(1:66)
      colnames(Hos_INF_1_store) <- seq(1:66)
      colnames(Dual_0_store) <- seq(1:66)
      colnames(Dual_1_store) <- seq(1:66)
    }
    # take the 95% quantiles for each group and format back into one table
    temp_quantiles_RSV_0 <- as.data.table(apply(Hos_RSV_0_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    temp_quantiles_RSV_0[,Quantile := c("p0.025","p0.5","p0.975")]
    RSV_0 <- data.table(melt(temp_quantiles_RSV_0, id="Quantile"))
    RSV_0 <- dcast.data.table(RSV_0, variable~Quantile, value.var = "value")
    RSV_0[,people := "RSV < 2"]
    
    temp_quantiles_RSV_1 <- as.data.table(apply(Hos_RSV_1_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    temp_quantiles_RSV_1[,Quantile := c("p0.025","p0.5","p0.975")]
    RSV_1 <- data.table(melt(temp_quantiles_RSV_1, id="Quantile"))
    RSV_1 <- dcast.data.table(RSV_1, variable~Quantile, value.var = "value")
    RSV_1[,people := "RSV 2-5"]
    
    temp_quantiles_INF_0 <- as.data.table(apply(Hos_INF_0_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    temp_quantiles_INF_0[,Quantile := c("p0.025","p0.5","p0.975")]
    INF_0 <- data.table(melt(temp_quantiles_INF_0, id="Quantile"))
    INF_0 <- dcast.data.table(INF_0, variable~Quantile, value.var = "value")
    INF_0[,people := "INF < 2"]
    
    temp_quantiles_INF_1 <- as.data.table(apply(Hos_INF_1_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    temp_quantiles_INF_1[,Quantile := c("p0.025","p0.5","p0.975")]
    INF_1 <- data.table(melt(temp_quantiles_INF_1, id="Quantile"))
    INF_1 <- dcast.data.table(INF_1, variable~Quantile, value.var = "value")
    INF_1[,people := "INF 2-5"]
    
    temp_quantiles_dual_0 <- as.data.table(apply(Dual_0_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    temp_quantiles_dual_0[,Quantile := c("p0.025","p0.5","p0.975")]
    dual_0 <- data.table(melt(temp_quantiles_dual_0, id="Quantile"))
    dual_0 <- dcast.data.table(dual_0, variable~Quantile, value.var = "value")
    dual_0[,people := "Dual < 2"]
    
    temp_quantiles_dual_1 <- as.data.table(apply(Dual_1_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    temp_quantiles_dual_1[,Quantile := c("p0.025","p0.5","p0.975")]
    dual_1 <- data.table(melt(temp_quantiles_dual_1, id="Quantile"))
    dual_1 <- dcast.data.table(dual_1, variable~Quantile, value.var = "value")
    dual_1[,people := "Dual 2-5"]
    
    quantile_table <- rbind(RSV_0,RSV_1,INF_0,INF_1, dual_0, dual_1)
    quantile_table$season <- season
    quantile_table <- na.omit(quantile_table)
    quantile_table$week_begin <- rep(fixed_sim[[season]]$start_week, 6)
    
    if(season ==1){
      all_seasons <- quantile_table
    } else {
      all_seasons <- rbind(all_seasons, quantile_table)
    }
  }
  all_seasons$variable <- as.numeric(as.character(all_seasons$variable))
  all_seasons$season <- as.numeric(as.character(all_seasons$season))
  # plot graph
  SIM_PLOT <- create_quantile_plot_together(quantile_table = all_seasons,
                                            fixed_sim = fixed_sim,
                                            rolling_mean = rolling_mean)
  return(SIM_PLOT)
  
}

# run the simulations and extract the quantiles
run_simulation_quantiles <- function(param_values, season){
  #Age group specific variables
  RSV_Sus <- c(1, 0.75, rep(0.65, (length(age_groups) - 2)))
  Pop_sizes <- population_numbers
  #Fixed parameters
  times <- seq(from = 1, to = nrow(multi_data[[season]])*7, by = 1)

    Parameters <- list(
    betaR = unname(param_values[,"betaR"]),
    betaI = unname(param_values[,"betaI"]),
    sig = unname(param_values[,"sig"]),
    gammaR =  0.1111,
    gammaI = 0.2632,
    rho = unname(param_values[,"rho"]),
    Rdetect2 = unname(param_values[,"Rdetect2"]),
    Rdetect5 = unname(param_values[,"Rdetect5"]),
    Idetect = unname(param_values[,"Idetect"]),
    RSV_Sus = RSV_Sus,
    Contact_Structure = population_contacts,
    num_grps = length(age_groups),
    inf_sus = unname(param_values[,"inf_sus"]),
    propR = unname(param_values[,"propR"]),
    propI = unname(param_values[,"propI"]),
    seedR = 0,
    seedI = 0,
    pop_sizes = population_numbers,
    trickleI = 0,#unname(param_values[,"trickleI"]),
    trickleR = 0,#unname(param_values[,"trickleR"]),
    detect_multiplier = unname(param_values[,"detect_multiplier"]), 
    # overdispersion = unname(param_values[,"overdispersion"]), 
    RSV_mult = unname(param_values[,"RSV_mult"])#, 
  #  Dual_mult = unname(param_values[,"Dual_mult"])
  )
  
  # if (detectI == T) {
  #   Parameters$detect_multiplier <- param_values[,"detect_multiplier"]
  # }
 
  initials <- calculate_init_state(population_numbers, Parameters)

  trajectory <- simulate_SIPR(theta_season = Parameters,
                              init_state = initials,
                              times = times)

  trajectory[, Hos_RSV_0 := RSV_incidence_1*Parameters[["Rdetect2"]]*Parameters[["RSV_mult"]]
  ][, Hos_INF_0 := INF_incidence_1*Parameters[["Idetect"]]
  ][, Hos_RSV_1 := RSV_incidence_2*Parameters[["Rdetect5"]]*Parameters[["RSV_mult"]]
  ][, Hos_INF_1 := INF_incidence_2*Parameters[["Idetect"]]*Parameters[["detect_multiplier"]]
  ][, Dual_0 := Dual_incidence_1*Parameters[["Rdetect2"]]*Parameters[["RSV_mult"]] #*Parameters[["Dual_mult"]]
  ][, Dual_1 := Dual_incidence_2*Parameters[["Rdetect5"]]]*Parameters[["RSV_mult"]]#*Parameters[["Dual_mult"]]
  return(trajectory)
  
}

# create the plot of the quantiles

create_quantile_plot_together <- function(quantile_table, fixed_sim,
                                          rolling_mean = 1){
  
  #combine all the seasons of real data
  for(i in 1:11){
    fixed_sim[[i]]$season_number <- i
  }
  df <- do.call(rbind.data.frame,fixed_sim)
  fixed_sim <- df[,c("start_week", "Hos_INF_1", "Hos_INF_0", "Hos_RSV_0",   
                     "Hos_RSV_1", "both_0", "both_1" )]
  
  fixed_sim <- data.table(melt(fixed_sim, id=c( "start_week")))
  fixed_sim[which(fixed_sim$variable=="Hos_RSV_0"), variable := "RSV < 2"]
  fixed_sim[which(fixed_sim$variable=="Hos_RSV_1"), variable := "RSV 2-5"]
  fixed_sim[which(fixed_sim$variable=="Hos_INF_0"), variable := "INF < 2"]
  fixed_sim[which(fixed_sim$variable=="Hos_INF_1"), variable := "INF 2-5"]
  fixed_sim[which(fixed_sim$variable=="both_0"), variable := "Dual < 2"]
  fixed_sim[which(fixed_sim$variable=="both_1"), variable := "Dual 2-5"]
  
  colnames(quantile_table) <- c("week_begin", "p0.025", "p0.5", "p0.975", "variable",
                                "season", "start_week")
  merged_table <- merge(quantile_table, fixed_sim, by = c("start_week","variable"))
  merged_table$start_week <- as.Date(merged_table$start_week)
  merged_table$variable <- fct_rev(merged_table$variable)

  merged_table
  PLOT <- ggplot(merged_table) +
    geom_line(aes(x=start_week, y=zoo::rollmean(value, k=rolling_mean, na.pad = T)),
              size=0.5)+
    geom_ribbon(aes(x=start_week,ymin=p0.025, ymax = p0.975, fill = variable), alpha = 0.8) +
    scale_fill_manual(values =c("#66c2a5","#66c2a5","#fc8d62","#fc8d62", "#8da0cb", "#8da0cb")) +
    geom_line(aes(x=start_week, y=p0.5, 
                  group = variable, colour = variable), size=0.5)+
    scale_colour_manual(values =c("#66c2a5","#66c2a5","#fc8d62","#fc8d62", "#8da0cb", "#8da0cb")) +
    facet_wrap(variable~., ncol = 1, scales="free_y")+
    labs(y="Incidence",
         colour = "Group", fill="Group") +
    theme_bw() +
    guides(shape = guide_legend(override.aes = list(size = 0.5)),
           color = guide_legend(override.aes = list(size = 0.5)))+
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
          legend.position = "none", strip.text.y = element_text(angle = 0), 
          strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(colour = 'black',hjust = 0),
          axis.text.x = element_text(vjust = 2),
    panel.spacing.x=unit(0.1, "lines"),panel.spacing.y=unit(-0, "lines"))+
    
    # scale_x_continuous(breaks = seq(1,by=52, length.out = 11), labels = 
    #                      c("Jan 2007", "Jan 2008", "Jan 2009", "Jan 2010", 
    #                        "Jan 2011", "Jan 2012", "Jan 2013",
    #                        "Jan 2014", "Jan 2015", "Jan 2016", "Jan 2017")) +
    geom_vline(xintercept = dates_to_cut, alpha= 0.5)
  
  return(PLOT)
  
}
dates_to_cut <- c(as.Date("2007-02-07"),
                  as.Date("2007-12-17"),
                  as.Date("2009-03-23"),
                  as.Date("2010-01-25"),
                  as.Date("2011-01-17"),
                  as.Date("2012-01-09"),
                  as.Date("2012-12-10"),
                  as.Date("2014-01-20"),
                  as.Date("2014-12-15"),
                  as.Date("2016-02-15"),
                  as.Date("2017-03-06"), 
                  as.Date("2018-01-01"))

# 
# quantiles_noflu_rel <- function(total_trace, fixed_sim, sample_size,
#                                    together = F,
#                                    rolling_mean = 1, from, to){
#   for(season in 1:11){
#     # Need to change numbers if number of parameters changes.
#     if (together == T){
#       trace_burned <- as.data.frame(total_trace[c(from:to),])
#     }else {
#       trace_burned <- as.data.frame(total_trace[ c(from:to),2:11])
#     }
#     # select a sample of 1000 parameter combinations
#     sample_trace <- Momocs::sample_n(tbl= trace_burned, size = sample_size, replace = T)
#     # for each combination
#     Hos_RSV_0_store <- data.frame(NULL)
#     Hos_RSV_1_store <- data.frame(NULL)
#     noflu_RSV_0_store <- data.frame(NULL)
#     noflu_RSV_1_store <- data.frame(NULL)
#     noflu_INF_0_store <- data.frame(NULL)
#     noflu_INF_1_store <- data.frame(NULL)
#     Hos_INF_0_store <- data.frame(NULL)
#     Hos_INF_1_store <- data.frame(NULL)
#     norsv_INF_0_store <- data.frame(NULL)
#     norsv_INF_1_store <- data.frame(NULL)
#     norsv_RSV_0_store <- data.frame(NULL)
#     norsv_RSV_1_store <- data.frame(NULL)
#     
#     for (r in 1:nrow(sample_trace)){
#       
#       # run the model
#       sample_trace$RSV_default <- 1
#       pars_out <-unlist(subset_parameters(sample_trace[r,], season))
#       pars_df <- matrix(pars_out,nrow=1)
#       colnames(pars_df) <- names(pars_out)
#       pars_df[,"sig"] = 1-pars_df[,"sig"] 
#       sim_out <- run_simulation_quantiles(pars_df,season)
#  
#       pars_out1 <- pars_out
#       pars_out1["propI"] <- 0
#       pars_df <- matrix(pars_out1,nrow=1)
#       colnames(pars_df) <- names(pars_out)
#       sim_out_noflu <- run_simulation_quantiles(pars_df,season)
#       pars_out2 <- pars_out
#       pars_out2["propR"] <- 0
#       pars_df <- matrix(pars_out2,nrow=1)
#       colnames(pars_df) <- names(pars_out)
#       sim_out_norsv <- run_simulation_quantiles(pars_df,season)
#       # save each one into dataframe with other samples
#       #each sample saved as row, so columns are timestep
#       Hos_RSV_0_store <- rbind(Hos_RSV_0_store,
#                                unname(unlist(sim_out[1:66,"Hos_RSV_0"])))
#       Hos_RSV_1_store <- rbind(Hos_RSV_1_store,
#                                unname(unlist(sim_out[1:66,"Hos_RSV_1"])))
#       noflu_RSV_0_store <- rbind(noflu_RSV_0_store,
#                                unname(unlist(sim_out_noflu[1:66,"Hos_RSV_0"])))
#       noflu_RSV_1_store <- rbind(noflu_RSV_1_store,
#                                unname(unlist(sim_out_noflu[1:66,"Hos_RSV_1"])))
#       Hos_INF_0_store <- rbind(Hos_INF_0_store,
#                                unname(unlist(sim_out[1:66,"Hos_INF_0"])))
#       Hos_INF_1_store <- rbind(Hos_INF_1_store,
#                                unname(unlist(sim_out[1:66,"Hos_INF_1"])))
#       norsv_INF_0_store <- rbind(norsv_INF_0_store,
#                                  unname(unlist(sim_out_norsv[1:66,"Hos_INF_0"])))
#       norsv_INF_1_store <- rbind(norsv_INF_1_store,
#                                  unname(unlist(sim_out_norsv[1:66,"Hos_INF_1"])))
#  
#       colnames(Hos_RSV_0_store) <- seq(1:66)
#       colnames(Hos_RSV_1_store) <- seq(1:66)
#       colnames(noflu_RSV_0_store) <- seq(1:66)
#       colnames(noflu_RSV_1_store) <- seq(1:66)
#       colnames(Hos_INF_0_store) <- seq(1:66)
#       colnames(Hos_INF_1_store) <- seq(1:66)
#       colnames(norsv_INF_0_store) <- seq(1:66)
#       colnames(norsv_INF_1_store) <- seq(1:66)
#     }
#     # take the 95% quantiles for each group and format back into one table
#     temp_quantiles_RSV_0 <- as.data.table(apply(Hos_RSV_0_store, 2, function(x)
#       quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
#     temp_quantiles_RSV_0[,Quantile := c("p0.025","p0.5","p0.975")]
#     RSV_0 <- data.table(melt(temp_quantiles_RSV_0, id="Quantile"))
#     RSV_0 <- dcast.data.table(RSV_0, variable~Quantile, value.var = "value")
#     RSV_0[,people := "A) RSV < 2"]
#     RSV_0[,type := "Current"]
#     
#     temp_quantiles_RSV_1 <- as.data.table(apply(Hos_RSV_1_store, 2, function(x)
#       quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
#     temp_quantiles_RSV_1[,Quantile := c("p0.025","p0.5","p0.975")]
#     RSV_1 <- data.table(melt(temp_quantiles_RSV_1, id="Quantile"))
#     RSV_1 <- dcast.data.table(RSV_1, variable~Quantile, value.var = "value")
#     RSV_1[,people := "B) RSV 2-5"]
#     RSV_1[,type := "Current"]
#     
#     noflu_quantiles_RSV_0 <- as.data.table(apply(noflu_RSV_0_store, 2, function(x)
#       quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
#     noflu_quantiles_RSV_0[,Quantile := c("p0.025","p0.5","p0.975")]
#     RSV_noflu_0 <- data.table(melt(noflu_quantiles_RSV_0, id="Quantile"))
#     RSV_noflu_0 <- dcast.data.table(RSV_noflu_0, variable~Quantile, value.var = "value")
#     RSV_noflu_0[,people := "A) RSV < 2"]
#     RSV_noflu_0[,type := "Vaccination"]
#     
#     noflu_quantiles_RSV_1 <- as.data.table(apply(noflu_RSV_1_store, 2, function(x)
#       quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
#     noflu_quantiles_RSV_1[,Quantile := c("p0.025","p0.5","p0.975")]
#     RSV_noflu_1 <- data.table(melt(noflu_quantiles_RSV_1, id="Quantile"))
#     RSV_noflu_1 <- dcast.data.table(RSV_noflu_1, variable~Quantile, value.var = "value")
#     RSV_noflu_1[,people := "B) RSV 2-5"]
#     RSV_noflu_1[,type := "Vaccination"]
# 
#     temp_quantiles_INF_0 <- as.data.table(apply(Hos_INF_0_store, 2, function(x)
#       quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
#     temp_quantiles_INF_0[,Quantile := c("p0.025","p0.5","p0.975")]
#     INF_0 <- data.table(melt(temp_quantiles_INF_0, id="Quantile"))
#     INF_0 <- dcast.data.table(INF_0, variable~Quantile, value.var = "value")
#     INF_0[,people := "C) INF < 2"]
#     INF_0[,type := "Current"]
#     
#     temp_quantiles_INF_1 <- as.data.table(apply(Hos_INF_1_store, 2, function(x)
#       quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
#     temp_quantiles_INF_1[,Quantile := c("p0.025","p0.5","p0.975")]
#     INF_1 <- data.table(melt(temp_quantiles_INF_1, id="Quantile"))
#     INF_1 <- dcast.data.table(INF_1, variable~Quantile, value.var = "value")
#     INF_1[,people := "D) INF 2-5"]
#     INF_1[,type := "Current"]
#     
#     norsv_quantiles_INF_0 <- as.data.table(apply(norsv_INF_0_store, 2, function(x)
#       quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
#     norsv_quantiles_INF_0[,Quantile := c("p0.025","p0.5","p0.975")]
#     INF_norsv_0 <- data.table(melt(norsv_quantiles_INF_0, id="Quantile"))
#     INF_norsv_0 <- dcast.data.table(INF_norsv_0, variable~Quantile, value.var = "value")
#     INF_norsv_0[,people := "C) INF < 2"]
#     INF_norsv_0[,type := "Vaccination"]
#     
#     norsv_quantiles_INF_1 <- as.data.table(apply(norsv_INF_1_store, 2, function(x)
#       quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
#     norsv_quantiles_INF_1[,Quantile := c("p0.025","p0.5","p0.975")]
#     INF_norsv_1 <- data.table(melt(norsv_quantiles_INF_1, id="Quantile"))
#     INF_norsv_1 <- dcast.data.table(INF_norsv_1, variable~Quantile, value.var = "value")
#     INF_norsv_1[,people := "D) INF 2-5"]
#     INF_norsv_1[,type := "Vaccination"]
# 
#     quantile_table <- rbind(RSV_0,RSV_1, RSV_noflu_0, RSV_noflu_1, 
#                             INF_0,INF_1, INF_norsv_0, INF_norsv_1)
#     quantile_table$season <- season
#     quantile_table <- na.omit(quantile_table)
#     quantile_table$week_begin <- rep(fixed_sim[[season]]$start_week, 8)
#     
#     if(season ==1){
#       all_seasons <- quantile_table
#     } else {
#       all_seasons <- rbind(all_seasons, quantile_table)
#     }
#   }
#   all_seasons$variable <- as.numeric(as.character(all_seasons$variable))
#   all_seasons$season <- as.numeric(as.character(all_seasons$season))
#   all_seasons$week_begin <- as.Date(all_seasons$week_begin)
#   browser()
#   
#   all_seasons[, sum(p0.5), by = c("type", "people")]
#   # plot graph]
#   SIM_PLOT <- ggplot(all_seasons) +
#     geom_line(aes(x=week_begin, y=p0.5, colour = type))+
#     geom_ribbon(aes(x=week_begin,ymin=p0.025, ymax = p0.975, 
#                     fill = type), alpha = 0.5) +
#     facet_wrap(people~., ncol = 1, scales = "free")+
#     labs(y="Incidence",
#          colour = "Group", fill="Group") +
#     theme_bw() +
#     # guides(shape = guide_legend(override.aes = list(size = 0.5)),
#     #        color = guide_legend(override.aes = list(size = 0.5)))+
#     # legend_key_width= unit(0.05,"cm"))+
#     theme(legend.title = element_text(size = 10),
#           legend.text  = element_text(size = 10),
#           axis.text = element_text(size=10),
#           axis.title = element_text(size=10, margin = margin(l=0,r=100,b=0,t=0)),
#           legend.key.width = unit(0.25,"cm"),
#           axis.title.x =  element_blank(),
#           # legend.margin = margin(c(1,1,1,1))
#           legend.key.height = unit(0.3,"cm"),
#           title = element_text(size=10),
#            strip.text.y = element_text(angle = 0), 
#           strip.background = element_rect(colour="white", fill="white"),
#           strip.text = element_text(colour = 'black',hjust = 0))+
#     geom_vline(xintercept = dates_to_cut, alpha= 0.5)
#   return(SIM_PLOT)
#   
# }
# 


# quantiles_noflu <- function(total_trace, fixed_sim, sample_size,
#                             together = F,
#                             rolling_mean = 1, from, to){
#   for(season in 1:11){
#     # Need to change numbers if number of parameters changes.
#     if (together == T){
#       trace_burned <- as.data.frame(total_trace[c(from:to),])
#     }else {
#       trace_burned <- as.data.frame(total_trace[ c(from:to),2:11])
#     }
#     # select a sample of 1000 parameter combinations
#     sample_trace <- Momocs::sample_n(tbl= trace_burned, size = sample_size, replace = T)
#     # for each combination
#     Hos_RSV_0_store <- data.frame(NULL)
#     Hos_RSV_1_store <- data.frame(NULL)
#     noflu_RSV_0_store <- data.frame(NULL)
#     noflu_RSV_1_store <- data.frame(NULL)
#     
#     for (r in 1:nrow(sample_trace)){
#       
#       # run the model
#       # not sure that this will work as the input probably a different format
#       pars_out <-unlist(subset_parameters(sample_trace[r,], season))
#       pars_df <- matrix(pars_out,nrow=1)
#       colnames(pars_df) <- names(pars_out)
#       sim_out <- run_simulation_quantiles(pars_df,season)
#       pars_out["propI"] <- 0
#       pars_out["trickleI"] <- 0
#       pars_df <- matrix(pars_out,nrow=1)
#       colnames(pars_df) <- names(pars_out)
#       sim_out_noflu <- run_simulation_quantiles(pars_df,season)
#       # save each one into dataframe with other samples
#       #each sample saved as row, so columns are timestep
#       Hos_RSV_0_store <- rbind(Hos_RSV_0_store,
#                                unname(unlist(sim_out[1:63,"Hos_RSV_0"])))
#       Hos_RSV_1_store <- rbind(Hos_RSV_1_store,
#                                unname(unlist(sim_out[1:63,"Hos_RSV_1"])))
#       noflu_RSV_0_store <- rbind(noflu_RSV_0_store,
#                                  unname(unlist(sim_out_noflu[1:63,"Hos_RSV_0"])))
#       noflu_RSV_1_store <- rbind(noflu_RSV_1_store,
#                                  unname(unlist(sim_out_noflu[1:63,"Hos_RSV_1"])))
#       
#       colnames(Hos_RSV_0_store) <- seq(1:63)
#       colnames(Hos_RSV_1_store) <- seq(1:63)
#       colnames(noflu_RSV_0_store) <- seq(1:63)
#       colnames(noflu_RSV_1_store) <- seq(1:63)
#     }
#     # take the 95% quantiles for each group and format back into one table
#     temp_quantiles_RSV_0 <- as.data.table(apply(Hos_RSV_0_store, 2, function(x)
#       quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
#     temp_quantiles_RSV_0[,Quantile := c("p0.025","p0.5","p0.975")]
#     RSV_0 <- data.table(melt(temp_quantiles_RSV_0, id="Quantile"))
#     RSV_0 <- dcast.data.table(RSV_0, variable~Quantile, value.var = "value")
#     RSV_0[,people := "RSV < 2"]
#     RSV_0[,type := "Current"]
#     
#     temp_quantiles_RSV_1 <- as.data.table(apply(Hos_RSV_1_store, 2, function(x)
#       quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
#     temp_quantiles_RSV_1[,Quantile := c("p0.025","p0.5","p0.975")]
#     RSV_1 <- data.table(melt(temp_quantiles_RSV_1, id="Quantile"))
#     RSV_1 <- dcast.data.table(RSV_1, variable~Quantile, value.var = "value")
#     RSV_1[,people := "RSV 2-5"]
#     RSV_1[,type := "Current"]
#     
#     noflu_quantiles_RSV_0 <- as.data.table(apply(noflu_RSV_0_store, 2, function(x)
#       quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
#     noflu_quantiles_RSV_0[,Quantile := c("p0.025","p0.5","p0.975")]
#     RSV_noflu_0 <- data.table(melt(noflu_quantiles_RSV_0, id="Quantile"))
#     RSV_noflu_0 <- dcast.data.table(RSV_noflu_0, variable~Quantile, value.var = "value")
#     RSV_noflu_0[,people := "RSV < 2"]
#     RSV_noflu_0[,type := "Vaccination"]
#     
#     noflu_quantiles_RSV_1 <- as.data.table(apply(noflu_RSV_1_store, 2, function(x)
#       quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
#     noflu_quantiles_RSV_1[,Quantile := c("p0.025","p0.5","p0.975")]
#     RSV_noflu_1 <- data.table(melt(noflu_quantiles_RSV_1, id="Quantile"))
#     RSV_noflu_1 <- dcast.data.table(RSV_noflu_1, variable~Quantile, value.var = "value")
#     RSV_noflu_1[,people := "RSV 2-5"]
#     RSV_noflu_1[,type := "Vaccination"]
#     
#     
#     quantile_table <- rbind(RSV_0,RSV_1, RSV_noflu_0, RSV_noflu_1)
#     quantile_table$season <- season
#     quantile_table <- na.omit(quantile_table)
#     quantile_table$week_begin <- rep(fixed_sim[[season]]$start_week, 4)
#     
#     if(season ==1){
#       all_seasons <- quantile_table
#     } else {
#       all_seasons <- rbind(all_seasons, quantile_table)
#     }
#   }
#   all_seasons$variable <- as.numeric(as.character(all_seasons$variable))
#   all_seasons$season <- as.numeric(as.character(all_seasons$season))
#   all_seasons$week_begin <- as.Date(all_seasons$week_begin)
#   # plot graph]
#   browser()
#   SIM_PLOT <- ggplot(all_seasons) +
#     geom_line(aes(x=week_begin, y=p0.5, colour = type))+
#     geom_ribbon(aes(x=week_begin,ymin=p0.025, ymax = p0.975, 
#                     fill = type), alpha = 0.5) +
#     facet_wrap(people~., ncol = 1, scales = "free")+
#     labs(y="Incidence",
#          colour = "Group", fill="Group") +
#     theme_bw() +
#     # guides(shape = guide_legend(override.aes = list(size = 0.5)),
#     #        color = guide_legend(override.aes = list(size = 0.5)))+
#     # legend_key_width= unit(0.05,"cm"))+
#     theme(legend.title = element_text(size = 10),
#           legend.text  = element_text(size = 10),
#           axis.text = element_text(size=10),
#           axis.title = element_text(size=10, margin = margin(l=0,r=100,b=0,t=0)),
#           legend.key.width = unit(0.25,"cm"),
#           axis.title.x =  element_blank(),
#           # legend.margin = margin(c(1,1,1,1))
#           legend.key.height = unit(0.3,"cm"),
#           title = element_text(size=10),
#           strip.text.y = element_text(angle = 0), 
#           strip.background = element_rect(colour="white", fill="white"),
#           strip.text = element_text(colour = 'black',hjust = 0))+
#     geom_vline(xintercept = dates_to_cut, alpha= 0.5)
#   return(SIM_PLOT)
#   
# }
# 
# plot_noflu <- function(quantile_table, fixed_sim,
#                        rolling_mean = 1){
#   
#   #combine all the seasons of real data
#   for(i in 1:11){
#     fixed_sim[[i]]$season_number <- i
#   }
#   df <- do.call(rbind.data.frame,fixed_sim)
#   fixed_sim <- df[,c("start_week", "Hos_INF_1", "Hos_INF_0", "Hos_RSV_0",   
#                      "Hos_RSV_1", "both_0", "both_1" )]
#   
#   fixed_sim <- data.table(melt(fixed_sim, id=c( "start_week")))
#   fixed_sim[which(fixed_sim$variable=="Hos_RSV_0"), variable := "RSV < 2"]
#   fixed_sim[which(fixed_sim$variable=="Hos_RSV_1"), variable := "RSV 2-5"]
#   fixed_sim[which(fixed_sim$variable=="Hos_INF_0"), variable := "INF < 2"]
#   fixed_sim[which(fixed_sim$variable=="Hos_INF_1"), variable := "INF 2-5"]
#   fixed_sim[which(fixed_sim$variable=="both_0"), variable := "Dual < 2"]
#   fixed_sim[which(fixed_sim$variable=="both_1"), variable := "Dual 2-5"]
#   
#   colnames(quantile_table) <- c("week_begin", "p0.025", "p0.5", "p0.975", "variable",
#                                 "season", "start_week")
#   merged_table <- merge(quantile_table, fixed_sim, by = c("start_week","variable"))
#   merged_table$start_week <- as.Date(merged_table$start_week)
#   merged_table$variable <- fct_rev(merged_table$variable)
#   
#   PLOT <- ggplot(merged_table) +
#     geom_line(aes(x=start_week, y=zoo::rollmean(value, k=rolling_mean, na.pad = T)),
#               size=0.5)+
#     geom_ribbon(aes(x=start_week,ymin=p0.025, ymax = p0.975), alpha = 0.5, 
#                 fill = "magenta") +
#     geom_line(aes(x=start_week, y=p0.5, 
#                   group = variable), colour = "magenta", size=0.3)+ 
#     facet_wrap(variable~., ncol = 1, scales = "free")+
#     labs(y="Incidence",
#          colour = "Group", fill="Group") +
#     theme_bw() +
#     guides(shape = guide_legend(override.aes = list(size = 0.5)),
#            color = guide_legend(override.aes = list(size = 0.5)))+
#     # legend_key_width= unit(0.05,"cm"))+
#     theme(legend.title = element_text(size = 10),
#           legend.text  = element_text(size = 10),
#           axis.text = element_text(size=10),
#           axis.title = element_text(size=10, margin = margin(l=0,r=100,b=0,t=0)),
#           legend.key.width = unit(0.25,"cm"),
#           axis.title.x =  element_blank(),
#           # legend.margin = margin(c(1,1,1,1))
#           legend.key.height = unit(0.3,"cm"),
#           title = element_text(size=10),
#           legend.position = "none", strip.text.y = element_text(angle = 0), 
#           strip.background = element_rect(colour="white", fill="white"),
#           strip.text = element_text(colour = 'black',hjust = 0),
#           axis.text.x = element_text(vjust = 2),
#           panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(-0.5, "lines"))+
#     
#     # scale_x_continuous(breaks = seq(1,by=52, length.out = 11), labels = 
#     #                      c("Jan 2007", "Jan 2008", "Jan 2009", "Jan 2010", 
#     #                        "Jan 2011", "Jan 2012", "Jan 2013",
#     #                        "Jan 2014", "Jan 2015", "Jan 2016", "Jan 2017")) +
#     geom_vline(xintercept = dates_to_cut, alpha= 0.5)
#   
#   return(PLOT)
#   
# }


quantiles_noflu <- function(total_trace, fixed_sim, sample_size,
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
  #    browser()
      # not sure that this will work as the input probably a different format
      pars_out <-unlist(subset_parameters(sample_trace[r,], season))
      if(season_params[season, "RSV_detect_multiplier"]=="RSV_default"){
        pars_season$RSV_mult <- 1
      } else if(season_params[season, "RSV_detect_multiplier"]=="RSV_detect_muliplier"){
        pars_season$RSV_mult <- input_params$RSV_detect_multiplier
      }
      pars_df <- matrix(pars_out,nrow=1)
      colnames(pars_df) <- names(pars_out)
      pars_df[,"sig"] = 1-pars_df[,"sig"] 
      sim_out <- run_simulation_quantiles(pars_df,season)
      pars_out["propI"] <- 0
      pars_out["trickleI"] <- 0
      pars_df <- matrix(pars_out,nrow=1)
      colnames(pars_df) <- names(pars_out)
      sim_out_noflu <- run_simulation_quantiles(pars_df,season)
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

      colnames(Hos_RSV_0_store) <- seq(1:63)
      colnames(Hos_RSV_1_store) <- seq(1:63)
      colnames(noflu_RSV_0_store) <- seq(1:63)
      colnames(noflu_RSV_1_store) <- seq(1:63)
    }
    # take the 95% quantiles for each group and format back into one table
    temp_quantiles_RSV_0 <- as.data.table(apply(Hos_RSV_0_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    temp_quantiles_RSV_0[,Quantile := c("p0.025","p0.5","p0.975")]
    RSV_0 <- data.table(melt(temp_quantiles_RSV_0, id="Quantile"))
    RSV_0 <- dcast.data.table(RSV_0, variable~Quantile, value.var = "value")
    RSV_0[,people := "RSV < 2"]
    RSV_0[,type := "Current"]

    temp_quantiles_RSV_1 <- as.data.table(apply(Hos_RSV_1_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    temp_quantiles_RSV_1[,Quantile := c("p0.025","p0.5","p0.975")]
    RSV_1 <- data.table(melt(temp_quantiles_RSV_1, id="Quantile"))
    RSV_1 <- dcast.data.table(RSV_1, variable~Quantile, value.var = "value")
    RSV_1[,people := "RSV 2-5"]
    RSV_1[,type := "Current"]

    noflu_quantiles_RSV_0 <- as.data.table(apply(noflu_RSV_0_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    noflu_quantiles_RSV_0[,Quantile := c("p0.025","p0.5","p0.975")]
    RSV_noflu_0 <- data.table(melt(noflu_quantiles_RSV_0, id="Quantile"))
    RSV_noflu_0 <- dcast.data.table(RSV_noflu_0, variable~Quantile, value.var = "value")
    RSV_noflu_0[,people := "RSV < 2"]
    RSV_noflu_0[,type := "Vaccination"]

    noflu_quantiles_RSV_1 <- as.data.table(apply(noflu_RSV_1_store, 2, function(x)
      quantile(x, probs = c(0.025,0.5,0.975), na.rm = T)))
    noflu_quantiles_RSV_1[,Quantile := c("p0.025","p0.5","p0.975")]
    RSV_noflu_1 <- data.table(melt(noflu_quantiles_RSV_1, id="Quantile"))
    RSV_noflu_1 <- dcast.data.table(RSV_noflu_1, variable~Quantile, value.var = "value")
    RSV_noflu_1[,people := "RSV 2-5"]
    RSV_noflu_1[,type := "Vaccination"]


    quantile_table <- rbind(RSV_0,RSV_1, RSV_noflu_0, RSV_noflu_1)
    quantile_table$season <- season
    quantile_table <- na.omit(quantile_table)
    quantile_table$week_begin <- rep(fixed_sim[[season]]$start_week, 4)

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

# plot_noflu <- function(quantile_table, fixed_sim,
#                        rolling_mean = 1){
#   
#   #combine all the seasons of real data
#   for(i in 1:11){
#     fixed_sim[[i]]$season_number <- i
#   }
#   df <- do.call(rbind.data.frame,fixed_sim)
#   fixed_sim <- df[,c("start_week", "Hos_INF_1", "Hos_INF_0", "Hos_RSV_0",   
#                      "Hos_RSV_1", "both_0", "both_1" )]
#   
#   fixed_sim <- data.table(melt(fixed_sim, id=c( "start_week")))
#   fixed_sim[which(fixed_sim$variable=="Hos_RSV_0"), variable := "RSV < 2"]
#   fixed_sim[which(fixed_sim$variable=="Hos_RSV_1"), variable := "RSV 2-5"]
#   fixed_sim[which(fixed_sim$variable=="Hos_INF_0"), variable := "INF < 2"]
#   fixed_sim[which(fixed_sim$variable=="Hos_INF_1"), variable := "INF 2-5"]
#   fixed_sim[which(fixed_sim$variable=="both_0"), variable := "Dual < 2"]
#   fixed_sim[which(fixed_sim$variable=="both_1"), variable := "Dual 2-5"]
#   
#   colnames(quantile_table) <- c("week_begin", "p0.025", "p0.5", "p0.975", "variable",
#                                 "season", "start_week")
#   merged_table <- merge(quantile_table, fixed_sim, by = c("start_week","variable"))
#   merged_table$start_week <- as.Date(merged_table$start_week)
#   merged_table$variable <- fct_rev(merged_table$variable)
#   
#   PLOT <- ggplot(merged_table) +
#     geom_line(aes(x=start_week, y=zoo::rollmean(value, k=rolling_mean, na.pad = T)),
#               size=0.5)+
#     geom_ribbon(aes(x=start_week,ymin=p0.025, ymax = p0.975), alpha = 0.5, 
#                 fill = "magenta") +
#     geom_line(aes(x=start_week, y=p0.5, 
#                   group = variable), colour = "magenta", size=0.3)+ 
#     facet_wrap(variable~., ncol = 1, scales = "free")+
#     labs(y="Incidence",
#          colour = "Group", fill="Group") +
#     theme_bw() +
#     guides(shape = guide_legend(override.aes = list(size = 0.5)),
#            color = guide_legend(override.aes = list(size = 0.5)))+
#     # legend_key_width= unit(0.05,"cm"))+
#     theme(legend.title = element_text(size = 10),
#           legend.text  = element_text(size = 10),
#           axis.text = element_text(size=10),
#           axis.title = element_text(size=10, margin = margin(l=0,r=100,b=0,t=0)),
#           legend.key.width = unit(0.25,"cm"),
#           axis.title.x =  element_blank(),
#           # legend.margin = margin(c(1,1,1,1))
#           legend.key.height = unit(0.3,"cm"),
#           title = element_text(size=10),
#           legend.position = "none", strip.text.y = element_text(angle = 0), 
#           strip.background = element_rect(colour="white", fill="white"),
#           strip.text = element_text(colour = 'black',hjust = 0),
#           axis.text.x = element_text(vjust = 2),
#           panel.spacing.x=unit(0.5, "lines"),panel.spacing.y=unit(-0.5, "lines"))+
#     
#     # scale_x_continuous(breaks = seq(1,by=52, length.out = 11), labels = 
#     #                      c("Jan 2007", "Jan 2008", "Jan 2009", "Jan 2010", 
#     #                        "Jan 2011", "Jan 2012", "Jan 2013",
#     #                        "Jan 2014", "Jan 2015", "Jan 2016", "Jan 2017")) +
#     geom_vline(xintercept = dates_to_cut, alpha= 0.5)
#   
#   return(PLOT)
#   
# }
# 
