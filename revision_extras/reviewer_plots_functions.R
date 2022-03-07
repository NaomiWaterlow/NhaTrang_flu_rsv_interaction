# quantiles plots dual 

# plot the quantiles of the fit: functions needed


# take samples from the trace and generate quantiles + plot
fit_all_ages_spaghetti <- function(total_trace, fixed_sim, sample_size,
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
    RSV_0_store <- data.frame(NULL)
    RSV_1_store <- data.frame(NULL)
    INF_0_store <- data.frame(NULL)
    INF_1_store <- data.frame(NULL)
    Dual_0_store <- data.frame(NULL)
    Dual_1_store <- data.frame(NULL)
    RSV_2_store <- data.frame(NULL)
    RSV_3_store <- data.frame(NULL)
    INF_2_store <- data.frame(NULL)
    INF_3_store <- data.frame(NULL)
    Dual_2_store <- data.frame(NULL)
    Dual_3_store <- data.frame(NULL)
    RSV_4_store <- data.frame(NULL)
    INF_4_store <- data.frame(NULL)
    Dual_4_store <- data.frame(NULL)
    
    for (r in 1:nrow(sample_trace)){
      # run the model
      # not sure that this will work as the input probably a different format
      pars_out <-unlist(subset_parameters(sample_trace[r,], season))
      pars_df <- matrix(pars_out,nrow=1)
      colnames(pars_df) <- names(pars_out)
      pars_df[,"sig"] = 1-pars_df[,"sig"] 
      
      sim_out <- run_simulation_quantiles(pars_df,season)
      
      
      #each sample saved as row, so columns are timestep
      RSV_0_store <- rbind(RSV_0_store,
                           unname(unlist(sim_out[1:66,"RSV_incidence_1"])))
      RSV_1_store <- rbind(RSV_1_store,
                           unname(unlist(sim_out[1:66,"RSV_incidence_2"])))
      INF_0_store <- rbind(INF_0_store,
                           unname(unlist(sim_out[1:66,"INF_incidence_1"])))
      INF_1_store <- rbind(INF_1_store,
                               unname(unlist(sim_out[1:66,"INF_incidence_2"])))
      Dual_0_store <- rbind(Dual_0_store,
                            unname(unlist(sim_out[1:66,"Dual_incidence_1"])))
      Dual_1_store <- rbind(Dual_1_store,
                            unname(unlist(sim_out[1:66,"Dual_incidence_2"])))
      RSV_2_store <- rbind(RSV_2_store,
                           unname(unlist(sim_out[1:66,"RSV_incidence_3"])))
      RSV_3_store <- rbind(RSV_3_store,
                           unname(unlist(sim_out[1:66,"RSV_incidence_4"])))
      INF_2_store <- rbind(INF_2_store,
                           unname(unlist(sim_out[1:66,"INF_incidence_3"])))
      INF_3_store <- rbind(INF_3_store,
                           unname(unlist(sim_out[1:66,"INF_incidence_4"])))
      Dual_2_store <- rbind(Dual_2_store,
                            unname(unlist(sim_out[1:66,"Dual_incidence_3"])))
      Dual_3_store <- rbind(Dual_3_store,
                            unname(unlist(sim_out[1:66,"Dual_incidence_4"])))
      RSV_4_store <- rbind(RSV_4_store,
                          unname(unlist(sim_out[1:66,"RSV_incidence_5"])))
      INF_4_store <- rbind(INF_4_store,
                           unname(unlist(sim_out[1:66,"INF_incidence_5"])))
      Dual_4_store <- rbind(Dual_4_store,
                            unname(unlist(sim_out[1:66,"Dual_incidence_5"])))

      colnames(RSV_0_store) <- seq(1:66)
      colnames(RSV_1_store) <- seq(1:66)
      colnames(INF_0_store) <- seq(1:66)
      colnames(INF_1_store) <- seq(1:66)
      colnames(Dual_0_store) <- seq(1:66)
      colnames(Dual_1_store) <- seq(1:66)
      colnames(RSV_2_store) <- seq(1:66)
      colnames(RSV_3_store) <- seq(1:66)
      colnames(INF_2_store) <- seq(1:66)
      colnames(INF_3_store) <- seq(1:66)
      colnames(Dual_2_store) <- seq(1:66)
      colnames(Dual_3_store) <- seq(1:66)
      colnames(RSV_4_store) <- seq(1:66)
      colnames(INF_4_store) <- seq(1:66)
      colnames(Dual_4_store) <- seq(1:66)

    }
    
    confidence_interval <- function(vector, interval) {
      # Standard deviation of sample
      vec_sd <- sd(vector)
      # Sample size
      n <- length(vector)
      # Mean of sample
      vec_mean <- mean(vector)
      # Error according to t distribution
      error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
      # Confidence interval as a vector
      result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
      return(result)
    }
    
    

    # take the 95% quantiles for each group and format back into one table
    temp_quantiles_RSV_0 <- as.data.table(apply(RSV_0_store, 2, function(x)
      c(quantile(x, probs = c(0.025), na.rm = T), mean(x), quantile(x, probs = c(0.975), na.rm = T))))
    temp_quantiles_RSV_0[,Quantile := c("p0.025","p0.5","p0.975")]
    RSV_0 <- data.table(melt(temp_quantiles_RSV_0, id="Quantile"))
    RSV_0 <- dcast.data.table(RSV_0, variable~Quantile, value.var = "value")
    RSV_0[,people := "< 2"]
    RSV_0[,infection := "RSV"]
    
    temp_quantiles_RSV_1 <- as.data.table(apply(RSV_1_store, 2, function(x)
      c(quantile(x, probs = c(0.025), na.rm = T), mean(x), quantile(x, probs = c(0.975), na.rm = T))))
    temp_quantiles_RSV_1[,Quantile := c("p0.025","p0.5","p0.975")]
    RSV_1 <- data.table(melt(temp_quantiles_RSV_1, id="Quantile"))
    RSV_1 <- dcast.data.table(RSV_1, variable~Quantile, value.var = "value")
    RSV_1[,people := "2-4"]
    RSV_1[,infection := "RSV"]
    
    temp_quantiles_RSV_2 <- as.data.table(apply(RSV_2_store, 2, function(x)
      c(quantile(x, probs = c(0.025), na.rm = T), mean(x), quantile(x, probs = c(0.975), na.rm = T))))
    temp_quantiles_RSV_2[,Quantile := c("p0.025","p0.5","p0.975")]
    RSV_2 <- data.table(melt(temp_quantiles_RSV_2, id="Quantile"))
    RSV_2 <- dcast.data.table(RSV_2, variable~Quantile, value.var = "value")
    RSV_2[,people := "5-15"]
    RSV_2[,infection := "RSV"]
    
    temp_quantiles_RSV_3 <- as.data.table(apply(RSV_3_store, 2, function(x)
      c(quantile(x, probs = c(0.025), na.rm = T), mean(x), quantile(x, probs = c(0.975), na.rm = T))))
    temp_quantiles_RSV_3[,Quantile := c("p0.025","p0.5","p0.975")]
    RSV_3 <- data.table(melt(temp_quantiles_RSV_3, id="Quantile"))
    RSV_3 <- dcast.data.table(RSV_3, variable~Quantile, value.var = "value")
    RSV_3[,people := "16-64"]
    RSV_3[,infection := "RSV"]
    
    temp_quantiles_RSV_4 <- as.data.table(apply(RSV_4_store, 2, function(x)
      c(quantile(x, probs = c(0.025), na.rm = T), mean(x), quantile(x, probs = c(0.975), na.rm = T))))
    temp_quantiles_RSV_4[,Quantile := c("p0.025","p0.5","p0.975")]
    RSV_4 <- data.table(melt(temp_quantiles_RSV_4, id="Quantile"))
    RSV_4 <- dcast.data.table(RSV_4, variable~Quantile, value.var = "value")
    RSV_4[,people := "65+"]
    RSV_4[,infection := "RSV"]
    
    temp_quantiles_INF_0 <- as.data.table(apply(INF_0_store, 2, function(x)
      c(quantile(x, probs = c(0.025), na.rm = T), mean(x), quantile(x, probs = c(0.975), na.rm = T))))
    temp_quantiles_INF_0[,Quantile := c("p0.025","p0.5","p0.975")]
    INF_0 <- data.table(melt(temp_quantiles_INF_0, id="Quantile"))
    INF_0 <- dcast.data.table(INF_0, variable~Quantile, value.var = "value")
    INF_0[,people := "< 2"]
    INF_0[,infection := "FLU"]
    
    temp_quantiles_INF_1 <- as.data.table(apply(INF_1_store, 2, function(x)
      c(quantile(x, probs = c(0.025), na.rm = T), mean(x), quantile(x, probs = c(0.975), na.rm = T))))
    temp_quantiles_INF_1[,Quantile := c("p0.025","p0.5","p0.975")]
    INF_1 <- data.table(melt(temp_quantiles_INF_1, id="Quantile"))
    INF_1 <- dcast.data.table(INF_1, variable~Quantile, value.var = "value")
    INF_1[,people := "2-4"]
    INF_1[,infection := "FLU"]
    
    temp_quantiles_INF_2 <- as.data.table(apply(INF_2_store, 2, function(x)
      c(quantile(x, probs = c(0.025), na.rm = T), mean(x), quantile(x, probs = c(0.975), na.rm = T))))
    temp_quantiles_INF_2[,Quantile := c("p0.025","p0.5","p0.975")]
    INF_2 <- data.table(melt(temp_quantiles_INF_2, id="Quantile"))
    INF_2 <- dcast.data.table(INF_2, variable~Quantile, value.var = "value")
    INF_2[,people := "5-15"]
    INF_2[,infection := "FLU"]
    
    temp_quantiles_INF_3 <- as.data.table(apply(INF_3_store, 2, function(x)
      c(quantile(x, probs = c(0.025), na.rm = T), mean(x), quantile(x, probs = c(0.975), na.rm = T))))
    temp_quantiles_INF_3[,Quantile := c("p0.025","p0.5","p0.975")]
    INF_3 <- data.table(melt(temp_quantiles_INF_3, id="Quantile"))
    INF_3 <- dcast.data.table(INF_3, variable~Quantile, value.var = "value")
    INF_3[,people := "16-64"]
    INF_3[,infection := "FLU"]
    
    temp_quantiles_INF_4 <- as.data.table(apply(INF_4_store, 2, function(x)
      c(quantile(x, probs = c(0.025), na.rm = T), mean(x), quantile(x, probs = c(0.975), na.rm = T))))
    temp_quantiles_INF_4[,Quantile := c("p0.025","p0.5","p0.975")]
    INF_4 <- data.table(melt(temp_quantiles_INF_4, id="Quantile"))
    INF_4 <- dcast.data.table(INF_4, variable~Quantile, value.var = "value")
    INF_4[,people := "65+"]
    INF_4[,infection := "FLU"]
    
    temp_quantiles_dual_0 <- as.data.table(apply(Dual_0_store, 2, function(x)
      c(quantile(x, probs = c(0.025), na.rm = T), mean(x), quantile(x, probs = c(0.975), na.rm = T))))
    temp_quantiles_dual_0[,Quantile := c("p0.025","p0.5","p0.975")]
    dual_0 <- data.table(melt(temp_quantiles_dual_0, id="Quantile"))
    dual_0 <- dcast.data.table(dual_0, variable~Quantile, value.var = "value")
    dual_0[,people := "< 2"]
    dual_0[,infection := "Dual"]
    
    temp_quantiles_dual_1 <- as.data.table(apply(Dual_1_store, 2, function(x)
      c(quantile(x, probs = c(0.025), na.rm = T), mean(x), quantile(x, probs = c(0.975), na.rm = T))))
    temp_quantiles_dual_1[,Quantile := c("p0.025","p0.5","p0.975")]
    dual_1 <- data.table(melt(temp_quantiles_dual_1, id="Quantile"))
    dual_1 <- dcast.data.table(dual_1, variable~Quantile, value.var = "value")
    dual_1[,people := "2-4"]
    dual_1[,infection := "Dual"]
    
    temp_quantiles_dual_2 <- as.data.table(apply(Dual_2_store, 2, function(x)
      c(quantile(x, probs = c(0.025), na.rm = T), mean(x), quantile(x, probs = c(0.975), na.rm = T))))
    temp_quantiles_dual_2[,Quantile := c("p0.025","p0.5","p0.975")]
    dual_2 <- data.table(melt(temp_quantiles_dual_2, id="Quantile"))
    dual_2 <- dcast.data.table(dual_2, variable~Quantile, value.var = "value")
    dual_2[,people := "5-15"]
    dual_2[,infection := "Dual"]
    
    temp_quantiles_dual_3 <- as.data.table(apply(Dual_3_store, 2, function(x)
      c(quantile(x, probs = c(0.025), na.rm = T), mean(x), quantile(x, probs = c(0.975), na.rm = T))))
    temp_quantiles_dual_3[,Quantile := c("p0.025","p0.5","p0.975")]
    dual_3 <- data.table(melt(temp_quantiles_dual_3, id="Quantile"))
    dual_3 <- dcast.data.table(dual_3, variable~Quantile, value.var = "value")
    dual_3[,people := "16-64"]
    dual_3[,infection := "Dual"]
    
    temp_quantiles_dual_4 <- as.data.table(apply(Dual_4_store, 2, function(x)
      c(quantile(x, probs = c(0.025), na.rm = T), mean(x), quantile(x, probs = c(0.975), na.rm = T))))
    temp_quantiles_dual_4[,Quantile := c("p0.025","p0.5","p0.975")]
    dual_4 <- data.table(melt(temp_quantiles_dual_4, id="Quantile"))
    dual_4 <- dcast.data.table(dual_4, variable~Quantile, value.var = "value")
    dual_4[,people := "65+"]
    dual_4[,infection := "Dual"]
   
    quantile_table <- rbind(RSV_0,RSV_1,INF_0,INF_1, dual_0, dual_1,
                            RSV_2,RSV_4,INF_2,INF_3, dual_2, dual_3,
                            RSV_3,INF_4,dual_4)
    quantile_table$season <- season
    quantile_table <- na.omit(quantile_table)
    quantile_table$week_begin <- rep(fixed_sim[[season]]$start_week, 15)
    
    if(season ==1){
      all_seasons <- quantile_table
    } else {
      all_seasons <- rbind(all_seasons, quantile_table)
    }
  }
  all_seasons$variable <- as.numeric(as.character(all_seasons$variable))
  all_seasons$season <- as.numeric(as.character(all_seasons$season))
  all_seasons$week_begin <- as.Date(all_seasons$week_begin)
  all_seasons$people <- factor(all_seasons$people, levels = 
                                 c("< 2", "2-4", "5-15", "16-64", "65+"))
  
  all_seasons[people == "< 2", p0.025 := p0.025/(population_numbers[1]/100000)]
  all_seasons[people == "2-4", p0.025 := p0.025/(population_numbers[2]/100000)]
  all_seasons[people == "5-15", p0.025 := p0.025/(population_numbers[3]/100000)]
  all_seasons[people == "16-64", p0.025 := p0.025/(population_numbers[4]/100000)]
  all_seasons[people == "65+", p0.025 := p0.025/(population_numbers[5]/100000)]
  
  all_seasons[people == "< 2", p0.5 := p0.5/(population_numbers[1]/100000)]
  all_seasons[people == "2-4", p0.5 := p0.5/(population_numbers[2]/100000)]
  all_seasons[people == "5-15", p0.5 := p0.5/(population_numbers[3]/100000)]
  all_seasons[people == "16-64", p0.5 := p0.5/(population_numbers[4]/100000)]
  all_seasons[people == "65+", p0.5 := p0.5/(population_numbers[5]/100000)]
  
  all_seasons[people == "< 2", p0.975 := p0.975/(population_numbers[1]/100000)]
  all_seasons[people == "2-4", p0.975 := p0.975/(population_numbers[2]/100000)]
  all_seasons[people == "5-15", p0.975 := p0.975/(population_numbers[3]/100000)]
  all_seasons[people == "16-64", p0.975 := p0.975/(population_numbers[4]/100000)]
  all_seasons[people == "65+", p0.975 := p0.975/(population_numbers[5]/100000)]
  

  SIM_PLOT <- ggplot(all_seasons, aes(x =week_begin, y = p0.5,  fill = people, colour = people)) + 
    facet_grid(infection~., scales = "free_y") + 
    geom_ribbon(aes(ymin=p0.025, ymax = p0.975), alpha = 0.3, colour = NA) + 
    geom_line(size = .75) + 

    theme_linedraw() + 
    labs(y = "Incidence per 100'000 population", x = "year",
    #     title=("Proportion infected over time by age group and infection type")
    )+
    geom_vline(xintercept = dates_to_cut, alpha= 0.5)
  
  # plot graph

  return(SIM_PLOT)
  
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
  browser()
  merged_table
  PLOT <- ggplot(merged_table) +
    geom_line(aes(x=start_week, y=zoo::rollmean(value, k=rolling_mean, na.pad = T)),
              size=0.5)+
    geom_ribbon(aes(x=start_week,ymin=p0.025, ymax = p0.975, fill = variable), alpha = 0.5) +
    scale_fill_manual(values =c("#66c2a5","#66c2a5","#fc8d62","#fc8d62", "#8da0cb", "#8da0cb")) +
    geom_line(aes(x=start_week, y=p0.5, 
                  group = variable, colour = variable), size=0.6)+
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
  ][, Dual_incidence_2 := (Dual_cases2 - shift(Dual_cases2, 1L, type = "lag"))
    ][, RSV_incidence_3 := (Rcases3 - shift(Rcases3, 1L, type = "lag"))
    ][, INF_incidence_3 := (Icases3 - shift(Icases3, 1L, type = "lag"))
    ][, RSV_incidence_4 := (Rcases4 - shift(Rcases4 ,1L, type = "lag"))
    ][, INF_incidence_4 := (Icases4 - shift(Icases4, 1L, type = "lag"))
    ][, Dual_incidence_3 := (Dual_cases3 - shift(Dual_cases3, 1L, type = "lag"))
    ][, Dual_incidence_4 := (Dual_cases4 - shift(Dual_cases4, 1L, type = "lag"))
      ][, RSV_incidence_5 := (Rcases5 - shift(Rcases5, 1L, type = "lag"))
      ][, INF_incidence_5 := (Icases5 - shift(Icases5, 1L, type = "lag"))
      ][, Dual_incidence_5 := (Dual_cases5 - shift(Dual_cases5, 1L, type = "lag"))]
  
  outall[1, c("RSV_incidence_1", "INF_incidence_1",
              "RSV_incidence_2", "INF_incidence_2", 
              "Dual_incidence_1", "Dual_incidence_2", 
              "RSV_incidence_3", "INF_incidence_3",
              "RSV_incidence_4", "INF_incidence_4", 
              "Dual_incidence_3", "Dual_incidence_4", 
              "RSV_incidence_5", "INF_incidence_5",
              "Dual_incidence_5"
              )] = 0
  
  
  #print(get_attack_rates(outall))
  
  #remove all but the incidence columsn
  outall[, setdiff(names(outall), c("RSV_incidence_1", "INF_incidence_1",
                                    "RSV_incidence_2", "INF_incidence_2", 
                                    "Dual_incidence_1", "Dual_incidence_2", 
                                    "RSV_incidence_3", "INF_incidence_3",
                                    "RSV_incidence_4", "INF_incidence_4", 
                                    "Dual_incidence_3", "Dual_incidence_4", 
                                    "RSV_incidence_5", "INF_incidence_5",
                                    "Dual_incidence_5")) := NULL][] 
  #Convert to weekly
  outall$week <- c(sapply(1:(nrow(outall)/7), function(x){return(rep(x,7))}))
  outall_week <- outall[, lapply(.SD, sum), by=list(week)]
  
  
  return(as.data.table(outall_week))
}