
#population_numbers_season <- population_numbers[[season_number]]

label_trace <- function(trace_in, detectI = F){

  if(detectI == F){
  colnames(trace_in) <- c("liklihood", names(init_theta)[-44], "Step_accepted","Swap_accepted_to",
                          "swap_attempted")
  } else if(detectI==T){
    colnames(trace_in) <- c("liklihood", names(init_theta), "Step_accepted","Swap_accepted_to",
                                                    "swap_attempted")
                            
                          }
  
  return(trace_in)
}

format_trace <- function(trace_in, thin = 1, burnin=0, keep_all = F){

  trace_dt <- data.table(trace_in)
  trace_dt[,"time_steps"] = timesteps
  trace_dt <- trace_dt[burnin:dim(trace_dt)[1],]
  trace_dt <- trace_dt[seq(1,dim(trace_dt)[1], by=thin),]
  if(keep_all == F){
  trace_dt[,c("Step_accepted", "Swap_accepted_to", "swap_attempted") := NULL]}

  return(trace_dt)
}

plot_trace <- function(trace_in, multi=T){
  if (multi == T){
  TRACE_PLOTS <- ggplot(as.data.frame(trace_in)) +
    geom_line(aes(y=value, x=time_steps, colour=factor(temp))) +
    facet_wrap(variable~., scales="free", ncol=2 ) +
    theme(
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 10),
      axis.text = element_text(size=10),
      axis.title = element_text(size=10, margin = margin(l=0,r=100,b=0,t=0)),
      title = element_text(size=10),
      legend.position = "none")
  } else {
    TRACE_PLOTS <- ggplot(as.data.frame(trace_in)) +
      geom_line(aes(y=value, x=time_steps)) +
      facet_wrap(variable~., scales="free", ncol=2) +
      theme(
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 10),
      axis.text = element_text(size=10),
      axis.title = element_text(size=10, margin = margin(l=0,r=100,b=0,t=0)),
      title = element_text(size=10))
  }

  return(TRACE_PLOTS)
}

swaps_accepted <- function(trace_in){

  num <- sum(trace_in[,"Swap_accepted_to"]!=0)
  denom <- sum(trace_in[,"swap_attempted"])

  rate <- num / denom
  return(rate)
}

plot_density <- function(trace_in, multi = T){
if (multi == T){
  DENS_PLOTS <- ggplot(as.data.frame(trace_combined)) +
    geom_density(aes(x=value, colour=factor(temp))) +
    facet_wrap(variable~., scales="free")
} else {
  DENS_PLOTS <- ggplot(as.data.frame(trace_combined)) +
    geom_density(aes(x=value)) +
    facet_wrap(variable~., scales="free")
}
  return(DENS_PLOTS)
}


transform_back <- function(parameters){
  
  # parameters[,sig := (exp(sig)/(1+exp(sig)))]
 

     # parameters[,betaR_1 := exp(betaR_1)]
     # parameters[,betaR_2 := exp(betaR_2)]
     # parameters[,betaI_1 := exp(betaI_1)]
     # parameters[,betaI_2 := exp(betaI_2)]
     
     # parameters[,Rdetect2 := exp(Rdetect2)]
     # parameters[,Rdetect5 := exp(Rdetect5)]
     # parameters[,Idetect := exp(Idetect)]
    
     parameters[,propR_1 := propR_1/10000]
     parameters[,propR_2 := propR_2/10000]
     parameters[,propR_3 := propR_3/10000]
     parameters[,propR_4 := propR_4/10000]
     parameters[,propR_5 := propR_5/10000]
     parameters[,propR_6 := propR_6/10000]
     parameters[,propR_7 := propR_7/10000]
     parameters[,propR_8 := propR_8/10000]
     parameters[,propR_9 := propR_9/10000]
     parameters[,propR_10 := propR_10/10000]
     parameters[,propR_11 := propR_11/10000]

     parameters[,propI_1 := propI_1/10000]
     parameters[,propI_2 := propI_2/10000]
     parameters[,propI_3 := propI_3/10000]
     parameters[,propI_4 := propI_4/10000]
     parameters[,propI_5 := propI_5/10000]
     parameters[,propI_6 := propI_6/10000]
     parameters[,propI_7 := propI_7/10000]
     parameters[,propI_8 := propI_8/10000]
     parameters[,propI_9 := propI_9/10000]
     parameters[,propI_10 := propI_10/10000]
     parameters[,propI_11 := propI_11/10000]
     

     
     parameters[,rho := exp(rho)]
     # parameters[,detect_multiplier :=exp(rho)]
  
  return(parameters)
}
