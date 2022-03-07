# profile likelihoods

# list of sigma values
sigma_values <- seq(0,1, by= 0.02)
profile_storage <- data.frame()
sample_size <- 10

for(mode in c("low", "moderate")){
  if(mode == "low"){input_trace <- trace_edited[[1]][sig <0.2]}
  if(mode == "moderate"){input_trace <- trace_edited[[1]][sig >= 0.2]}
  
  for(sig_val in sigma_values){
    
    # select a sample of 1000 parameter combinations
    sample_trace <- Momocs::sample_n(tbl= input_trace, size = sample_size, replace = T)
    sample_trace[,"RSV_default"] <- 1
    data_to_fit <- multi_data
    for(sample_num in 1:sample_size){
      
      pars_transformed <- sample_trace[sample_num,]
      pars_transformed[["sig"]] <- sig_val
      
      season_out <-  sum(sapply(1:length(data_to_fit), function(season){
        #subset parameters
        pars_season <- subset_parameters(as.data.frame(pars_transformed), season)
        data_season <- data_to_fit[[season]]
        #calculate initial state
        init_state_season <- calculate_init_state(population_numbers, pars_season)
        #calculate likelihood of paraeters
        liklihood_init <- calc_likelihood_overall(pars_season,
                                                  init_state_season,
                                                  data_season,
                                                  log = T)
        
        return(if(is.finite(liklihood_init)){liklihood_init}else{NaN})
      }))
      
      profile_storage <-  rbind(profile_storage, c(sig_val,season_out, sample_num, mode))
    }
  }
}

profile_storage <- as.data.table(profile_storage)
colnames(profile_storage) <- c("sigma", "ll", "sample", "mode")

profile_storage$ll <- as.numeric(profile_storage$ll)

sum_profile <- profile_storage[,quantile(ll, probs = c(0.025), na.rm=T), by = c("sigma", "mode") ]
colnames(sum_profile)[3] <- "lower"
sum_profile$middle <-  profile_storage[,quantile(ll, probs = c(0.5), na.rm=T), by = c("sigma", "mode") ]$V1
sum_profile$upper <-  profile_storage[,quantile(ll, probs = c(0.975), na.rm=T), by = c("sigma", "mode") ]$V1
sum_profile$sigma <- as.numeric(sum_profile$sigma)

sum_profile[,sigma_lab := 1 -sigma]

# RHO_PROF <- ggplot(sum_profile, aes(x = sigma_lab, y = middle, fill = mode)) +
#   geom_line() + 
#   geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5) + 
#   theme_linedraw() + 
#   geom_hline(yintercept = max(sum_profile$middle, na.rm=T), linetype="dotted") + 
#   labs("loglikelihood", x = "1/duration of protection")

SIGMA_PROF <- ggplot(sum_profile, aes(x = sigma_lab, y = middle, fill = mode)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5) + 
  theme_linedraw() + 
  geom_hline(yintercept = max(sum_profile$middle, na.rm=T), linetype="dotted") + 
  labs(y = "Log likelihood", x = "Interaction parameter", fill= "Mode")


tiff(here("Figures","SIGMA_PROF.tiff"), height = 2000, width = 3200, res = 300)
SIGMA_PROF
dev.off()
