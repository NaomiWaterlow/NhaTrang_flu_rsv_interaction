quantiles_noflu_rel <- function(total_trace, fixed_sim, sample_size,
                                together = F,
                                rolling_mean = 1, from, to){

  # select a sample of 1000 parameter combinations
  # Need to change numbers if number of parameters changes.
  if (together == T){
    trace_burned <- as.data.frame(total_trace[c(from:to),])
  }else {
    trace_burned <- as.data.frame(total_trace[ c(from:to),2:11])
  }
  
  sample_trace <- Momocs::sample_n(tbl= trace_burned, size = sample_size, replace = T)
  # for eaach sampl
  all_output_store <- as.data.frame(matrix(nrow = sample_size , ncol = 13))
  for (r in 1:nrow(sample_trace)){
    # for each season
    noflu_flu_0_store <- 0
    noflu_flu_1_store <- 0
    noflu_rsv_0_store <- 0
    noflu_rsv_1_store <- 0
    norsv_flu_0_store <- 0
    norsv_flu_1_store <- 0
    norsv_rsv_0_store <- 0
    norsv_rsv_1_store <- 0
    original_flu_0_store <- 0
    original_flu_1_store <- 0
    original_rsv_0_store <- 0
    original_rsv_1_store <- 0
    
   for(season in 1:11){

      # run the model
      sample_trace$RSV_default <- 1
      pars_out <-unlist(subset_parameters(sample_trace[r,], season))
      pars_df <- matrix(pars_out,nrow=1)
      colnames(pars_df) <- names(pars_out)
      pars_df[,"sig"] = 1-pars_df[,"sig"] 
      sim_out <- run_simulation_quantiles(pars_df,season)
      
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

    
    noflu_flu_0_store <- noflu_flu_0_store + sum(sim_out_noflu$Hos_INF_0) + sum(sim_out_noflu$Dual_0)
    noflu_flu_1_store <- noflu_flu_1_store + sum(sim_out_noflu$Hos_INF_1) + sum(sim_out_noflu$Dual_1)
    noflu_rsv_0_store <- noflu_rsv_0_store + sum(sim_out_noflu$Hos_RSV_0) + sum(sim_out_noflu$Dual_0)
    noflu_rsv_1_store <- noflu_rsv_1_store + sum(sim_out_noflu$Hos_RSV_1) + sum(sim_out_noflu$Dual_1)
    
    norsv_flu_0_store <- norsv_flu_0_store + sum(sim_out_norsv$Hos_INF_0) + sum(sim_out_norsv$Dual_0)
    norsv_flu_1_store <- norsv_flu_1_store + sum(sim_out_norsv$Hos_INF_1) + sum(sim_out_norsv$Dual_1)
    norsv_rsv_0_store <- norsv_rsv_0_store + sum(sim_out_norsv$Hos_RSV_0) + sum(sim_out_norsv$Dual_0)
    norsv_rsv_1_store <- norsv_rsv_1_store + sum(sim_out_norsv$Hos_RSV_1) + sum(sim_out_norsv$Dual_1)
                                       
    original_flu_0_store <- original_flu_0_store + sum(sim_out$Hos_INF_0) + sum(sim_out$Dual_0)
    original_flu_1_store <- original_flu_1_store + sum(sim_out$Hos_INF_1) + sum(sim_out$Dual_1)
    original_rsv_0_store <- original_rsv_0_store + sum(sim_out$Hos_RSV_0) + sum(sim_out$Dual_0)
    original_rsv_1_store <- original_rsv_1_store + sum(sim_out$Hos_RSV_1) + sum(sim_out$Dual_1)
  
   }
    
    all_output_store[r,] <- c(noflu_flu_0_store, 
                              noflu_flu_1_store, 
                              noflu_rsv_0_store, 
                              noflu_rsv_1_store, 
                              norsv_flu_0_store, 
                              norsv_flu_1_store, 
                              norsv_rsv_0_store, 
                              norsv_rsv_1_store,
                              original_flu_0_store,
                              original_flu_1_store,
                              original_rsv_0_store,
                              original_rsv_1_store,
                              r) 
    
  }

  colnames(all_output_store) <- c("noflu_flu_0",
                                  "noflu_flu_1",
                                  "noflu_rsv_0",
                                  "noflu_rsv_1",
                                  "norsv_flu_0",
                                  "norsv_flu_1",
                                  "norsv_rsv_0",
                                  "norsv_rsv_1",
                                  "original_flu_0",
                                  "original_flu_1",
                                  "original_rsv_0",
                                  "original_rsv_1", 
                                  "sample")
  
  all_output_store$noflu_change_rsv <- ((all_output_store$original_rsv_0 + all_output_store$original_rsv_1 )-
    (all_output_store$noflu_rsv_0 +all_output_store$noflu_rsv_1)) / (all_output_store$original_rsv_0 + all_output_store$original_rsv_1 )
  
  all_output_store$norsv_change_flu <- ((all_output_store$original_flu_0 +  all_output_store$original_flu_1 )  - 
    (all_output_store$norsv_flu_0 +all_output_store$norsv_flu_1)) / (all_output_store$original_flu_0 +  all_output_store$original_flu_1 )
  
  
 quantile_table <- apply(all_output_store, 2, function(x)quantile(x, probs = c(0.025,0.5,0.975), na.rm = T))
 
 print(paste0("vaccination against flu increase rsv by ", -quantile_table[2,"noflu_change_rsv"]))
 print(paste0("CIs are ", -quantile_table[1,"noflu_change_rsv"], " to ", -quantile_table[3,"noflu_change_rsv"]))
 
 print(paste0("vaccination against rsv increase flu by ", -quantile_table[2,"norsv_change_flu"]))
 print(paste0("CIs are ", -quantile_table[1,"norsv_change_flu"], " to ", -quantile_table[3,"norsv_change_flu"]))
 
 quantile_subset <- data.frame(matrix(ncol = 5, nrow =12 ))
 quantile_subset[1,] <- c(quantile_table[,1],"Flu", "INF < 2")
 quantile_subset[2,] <- c(quantile_table[,2],"Flu", "INF 2-4")
 quantile_subset[3,] <- c(quantile_table[,3],"Flu", "RSV < 2")
 quantile_subset[4,] <- c(quantile_table[,4],"Flu", "RSV 2-4")
 quantile_subset[5,] <- c(quantile_table[,5],"RSV", "INF < 2")
 quantile_subset[6,] <- c(quantile_table[,6],"RSV", "INF 2-4")
 quantile_subset[7,] <- c(quantile_table[,7],"RSV", "RSV < 2")
 quantile_subset[8,] <- c(quantile_table[,8],"RSV", "RSV 2-4")
 quantile_subset[9,] <- c(quantile_table[,9],"None", "INF < 2")
 quantile_subset[10,] <- c(quantile_table[,10],"None", "INF 2-4")
 quantile_subset[11,] <- c(quantile_table[,11],"None", "RSV < 2")
 quantile_subset[12,] <- c(quantile_table[,12],"None", "RSV 2-4")
colnames(quantile_subset) <- c("p0.025", "p0.5", "p0.975", "Vaccination", "Type")
 
quantile_subset[,1] <- as.numeric(as.character(quantile_subset[,1] ))
quantile_subset[,2] <- as.numeric(as.character(quantile_subset[,2] ))
quantile_subset[,3] <- as.numeric(as.character(quantile_subset[,3] ))

PLOT <- ggplot(quantile_subset, aes(x = Type, y = p0.5, colour = Vaccination)) + 
  geom_point(position = position_dodge(width = 0.2)) + 
  geom_errorbar(aes(ymin = p0.025, ymax = p0.975),
                width = 0.1, position = position_dodge(width = 0.2)) +
  theme_linedraw() + 
  scale_color_manual(values = c("#fc8d62","black","#66c2a5"))+
  labs(x = "group", y = "Number of cases over all seasons", colour = "Vaccination", 
       title = "E: Vaccination scenarios") 
 return(PLOT)
  }
