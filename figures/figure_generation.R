# Create figures for paper!
temperatures <- rev(c(1:12))
library(coda)
set.seed(100)
##### READ IN REQUIRED TRACES #####
load(here::here("combined_dualmult_16_2021-08-17.Rdata"))
trace_1 <- total_trace[1:(length(total_trace))]
trace_1 <- lapply(trace_1, function(i) na.omit(i))
timesteps <- 1:dim(trace_1[[1]])[1]
trace_1 <- lapply(trace_1, function(i) label_trace(i, detectI = T))
trace_1_thinned <- trace_1[[1]][seq(1,nrow(trace_1[[1]]), by = 1),]
trace_edited_1 <- lapply(trace_1, function(i) format_trace(i, thin=10, burnin=250000,
                                                       keep_all = F))
trace_edited_1 <- lapply(trace_edited_1, function(i) transform_back(i))
trace_edited_1 <- Map(cbind, trace_edited_1, temp = temperatures)


 trace_to_sample <- trace_edited_1[[1]]
source("figure_generation_functions.R")

####### FIGURE 2 - Data/Fit ######
# creates fit plot
FITS_TOGETHER <- fit_quantiles_together(trace_to_sample, multi_data, 
                                       sample_size = 1000, together=T,
                                       rolling_mean = 1,  
                                       from=1, to=19651)
FITS_TOGETHER <- FITS_TOGETHER + labs(title="A: Model fit")

VOID <- ggplot()+ theme_void() + 
  labs(title="C: Model Diagram") + 
  theme(plot.title = element_text(hjust = 0.05)) 


###### SUPPLEMENTARY FIGURES DENSITY #####
trace_single <- trace_to_sample
trace_single_mt <- melt(trace_single, id=c("time_steps"))

DENSITY_PLOT <- 
  ggplot(trace_single_mt, aes( x = value)) + 
  facet_wrap(variable~., scales = "free", nrow=8) + 
  geom_density() + 
  theme_classic() + 
  theme(legend.position = "none", strip.text.y = element_text(angle = 0), 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(colour = 'black',hjust = 0)) + 
  labs(y = "Density", x = "Value") +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

DENSITY_PLOT

tiff(here("Figures","Sup_density.tiff"), height = 2000, width = 4000, res = 300)
DENSITY_PLOT
dev.off()
 



####### SUPPLEMENTARY Exponential susceptibility #######
sus_profile <- matrix( nrow = 10, ncol = 7)
sus_profile[1,] <- c(1,calc_start_sus(1),1) 
sus_profile[2,] <- c(2,calc_start_sus(0.9),0.9) 
sus_profile[3,] <- c(3,calc_start_sus(0.8),0.8) 
sus_profile[4,] <- c(4,calc_start_sus(0.7),0.7) 
sus_profile[5,] <- c(5,calc_start_sus(0.6),0.6) 
sus_profile[6,] <- c(6,calc_start_sus(0.5),0.5) 
sus_profile[7,] <- c(7,calc_start_sus(0.4),0.4) 
sus_profile[8,] <- c(8,calc_start_sus(0.3),0.3) 
sus_profile[9,] <- c(9,calc_start_sus(0.2),0.2) 
sus_profile[10,] <- c(10,calc_start_sus(0.1),0.1) 
sus_profile[10,] <- c(10,calc_start_sus(0),0) 
colnames(sus_profile) <- c("ID", "0-1", "2-4", "5-15",
                           "16-64" ,"64+", "parameter_value")

sus_profile_2 <- melt(data = data.table(sus_profile), 
                      id.vars=c("ID","parameter_value"))
sus_profile_2$value <- as.numeric(as.character(sus_profile_2$value))
sus_profile_2$parameter_value <- as.factor(sus_profile_2$parameter_value)
SUS_PLOT <- ggplot(sus_profile_2, aes(x=variable, y=value, 
                                      group = ID, colour = parameter_value)) +
  geom_line() +
  theme_linedraw() + 
  labs(y = "Susceptibility", colour = "Parameter value", x ="Age group") 
  
tiff(here("Figures","Sup_susceptibility.tiff"), height = 2000, width = 3200, res = 300)
SUS_PLOT
dev.off()




####### FIGURE 3: Attack rates and relative seasonal susceptibility #######

# take X samples 

attack_store <- data.table()
case_store <- data.table()
n_samples <- 50
for(j in c(1:n_samples)){
  trace_to_sample_conv <- trace_to_sample
    trace_to_sample_conv[,"sig"] = 1-trace_to_sample_conv[,"sig"] 
  sample_num <- sample(x=1:nrow(trace_to_sample_conv), size = 1)
  input_params <- unlist(trace_to_sample_conv[sample_num,2:46])
  input_params <- as.list(input_params)
  # Likelihood by season
  for (season in c(1:11)){
    #subset parameters
  season_lab <- seq(2007,2017, by =1)[season]
    pars_season <- subset_parameters(input_params, season)
    if(season_params[season, "RSV_detect_multiplier"]=="RSV_default"){
      pars_season$RSV_mult <- 1
    } else if(season_params[season, "RSV_detect_multiplier"]=="RSV_detect_muliplier"){
      pars_season$RSV_mult <- input_params$RSV_detect_multiplier
    }
 #   pars_season$RSV_mult <- season_params[season, "RSV_detect_multiplier"]
    data_season <- multi_data[[season]]
    init_state <- calculate_init_state(population_numbers,
                                       pars_season)
    rates <- run_model_attack(theta_season = pars_season,
                                                   init_state = init_state,
                                                   season_data = data_season)
    attack_rates <- as.data.frame(rates[[1]])
    case_totals <- as.data.frame(rates[[2]])
    
    attack_rates$season <- season_lab
    attack_rates$sample <- sample_num
    attack_rates$age <- factor(rownames(attack_rates), 
                               levels = rownames(attack_rates))
    
    case_totals$season <- season_lab
    case_totals$sample <- sample_num
    case_totals$age <- factor(rownames(case_totals), 
                               levels = rownames(case_totals))
    
    attack_store <- rbind(attack_store, attack_rates)
    case_store <- rbind(case_store, case_totals)
  }}

cc <- scales::seq_gradient_pal("darkmagenta",  "darkorange", "Lab")(seq(0,1,length.out=5))

attack_store_m <- melt(attack_store, id.vars = c("season", "sample", "age"))
attack_store_m$season <- as.factor(attack_store_m$season)
ATTACK_RATES <- ggplot(attack_store_m, aes(x = season, y = value, colour = age)) + 
  facet_wrap(variable~., scale = "free") + 
  geom_jitter(width = 0.2, alpha = 0.6) +
  coord_flip() +
  theme_linedraw() + 
  scale_colour_manual(values=cc) +
  labs(x= "Season Number", y = "Season Attack Rates", colour = "Age group", 
       title = "A: Attack Rates")

summarise_season <- function(i){
  a_season <- data.table(i)
  a_season_m <- melt(a_season, id.vars=c("start_week", "section", "week", "year"))
  summary <- a_season_m[, sum(value), by =variable]
  return(summary$V1)
}

summary <- as.data.table(sapply(multi_data, function(i) summarise_season(i)))
summary$variable <- c("FLU", "FLU", "RSV", "RSV", "Dual", "Dual")
summary$age <- c("0-1", "2-4","0-1", "2-4", "0-1", "2-4")  
colnames(summary) <- c(as.character(2007:2017), "variable", "age")
summary_m <- melt.data.table(summary, id.vars = c("variable", "age"))
colnames(summary_m) <- c("variable",  "age", "season", "value")
case_store_m <- melt(case_store, id.vars= c("sample", "age", "season"))
case_store_m$season <- as.factor(case_store_m$season)
combined_table <- merge(case_store_m, summary_m, by = c("variable", "age", "season") )
colnames(combined_table) <- c("virus", "age", "season", "sample", "Model", "Data")
combined_table$virus <- factor(combined_table$virus, 
                               levels = c("RSV", "FLU", "Dual"))

CASE <- ggplot(combined_table, aes(x = Data, y = Model, colour =  virus, shape = age )) + 
  geom_abline(intercept = 0, slope = 1) +
  geom_point() + 
  theme_linedraw() +
  labs(x = "Observed Cases", y = "Modelled Cases", colour = "Virus", shape = "Age group") + 
  scale_colour_manual(values=c("#66c2a5","#fc8d62", "#8da0cb")) +
  labs(title = "D: Goodness of fit")
  


# Susceptibility by age
cols <- as.numeric(grep("*sus*", colnames(trace_to_sample)))
 sus_tab <- trace_to_sample[, ..cols]
 colnames(sus_tab) <- as.character(seq(2007,2017, by =1))
 sus_tab_m <- melt(sus_tab)
 
 SUS <- ggplot(sus_tab_m, aes(x =  value, fill = variable)) +
   geom_density(alpha = 0.3)  + 
   theme_linedraw() + 
   theme(axis.text.y=element_blank(),
         axis.ticks.y=element_blank()) +
   labs(x = "Susceptibility parameter", y = "Density", fill = "Season", 
        title = "B: Influenza Susceptibility") + 
   guides(fill = guide_legend(override.aes = list(size = 0.1))) +
   theme(legend.title = element_text(size = 8), 
         legend.text = element_text(size = 8),
         legend.key.size = unit(0.7,"line"))
 
 
 sus_store <- data.frame(matrix(nrow = 11, ncol = 6))
 sus_store[1,] <-c(2007,calc_start_sus(quantile(trace_to_sample$inf_sus_1, probs = c(0.5))))
 sus_store[2,] <-c(2008,calc_start_sus(quantile(trace_to_sample$inf_sus_2, probs = c(0.5))))
 sus_store[3,] <-c(2009,calc_start_sus(quantile(trace_to_sample$inf_sus_3, probs = c(0.5))))
 sus_store[4,] <-c(2010,calc_start_sus(quantile(trace_to_sample$inf_sus_4, probs = c(0.5))))
 sus_store[5,] <-c(2011,calc_start_sus(quantile(trace_to_sample$inf_sus_5, probs = c(0.5))))
 sus_store[6,] <-c(2012,calc_start_sus(quantile(trace_to_sample$inf_sus_6, probs = c(0.5))))
 sus_store[7,] <-c(2013,calc_start_sus(quantile(trace_to_sample$inf_sus_7, probs = c(0.5))))
 sus_store[8,] <-c(2014,calc_start_sus(quantile(trace_to_sample$inf_sus_8, probs = c(0.5))))
 sus_store[9,] <-c(2015,calc_start_sus(quantile(trace_to_sample$inf_sus_9, probs = c(0.5))))
 sus_store[10,] <-c(2016,calc_start_sus(quantile(trace_to_sample$inf_sus_10, probs = c(0.5))))
 sus_store[11,] <-c(2017,calc_start_sus(quantile(trace_to_sample$inf_sus_11, probs = c(0.5))))
 colnames(sus_store) <- c("Year", "0-1", "2-4", "5-15","16-64", "65+" )
 sus_store_m <- melt(sus_store, id.vars = "Year")
 sus_store_m$Year <- factor(sus_store_m$Year, levels = unique(sus_store_m$Year))
 
 SUS2 <- ggplot(sus_store_m, aes(x = Year, y = value, colour = variable, group =variable)) + 
   geom_point() +
   geom_line() +
   scale_colour_manual(values=cc) + 
   theme_linedraw() + 
   labs(y = "Proportion Susceptible", x="", title = "B: Proportion susceptible to influenza by season") + 
   theme(legend.position = "none")
  
 
 SUS2
 
 g_legend<-function(a.gplot){
   tmp <- ggplot_gtable(ggplot_build(a.gplot))
   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
   legend <- tmp$grobs[[leg]]
   return(legend)}
 
shared_legend <- g_legend(ATTACK_RATES)

p3 <- grid.arrange(arrangeGrob(ATTACK_RATES + theme(legend.position="none"),
                               SUS2, nrow=2),
                   shared_legend, nrow=1,widths=c(10, 1))
 
 tiff(here("Figures","FIG3_Sus"), height = 2000, width = 3200, res = 300)
 grid.arrange(arrangeGrob(ATTACK_RATES + theme(legend.position="none"),
                          SUS2, nrow=2, heights=c(2,1)),
              shared_legend, nrow=1,widths=c(10, 1))
 dev.off()
 
 
 
 trace_to_sample <- data.table(trace_to_sample)
 trace_to_sample <-trace_to_sample[,1:45]
 
 trace_to_sample[(sig >0.2), section := "sigma > 0.2"]
 trace_to_sample[(sig <= 0.2), section := "sigma <= 0.2"]
 
 trace_to_sample_m <- melt.data.table(trace_to_sample)
 trace_to_sample_m[variable == "time_steps",] <- NA
 trace_to_sample_m <- na.omit(trace_to_sample_m)
 
 library(ggstance)
 
 ALL <- ggplot(trace_to_sample_m, aes(x = value, fill=section, colour =section)) + 
  # geom_histogram(bins =1000, position = "stack") + 
  geom_density() +
    facet_wrap(variable~., scales = "free") + 
   labs(x = "Parameter value", y = "Count", title = "A: Posterior estimates") + 
   theme_linedraw(
   )
 
 ALL <- ggplot(trace_to_sample_m, aes(x = value)) + 
   geom_histogram(bins =250, position = "stack", fill = "deepskyblue") + 
   facet_wrap(variable~., scales = "free") + 
   labs(x = "Parameter value", y = "Count", title = "A: Posterior estimates") + 
   theme_linedraw() + 
   geom_pointrangeh(data = estimates_mod, aes(x= p0.5, xmin = p0.025, xmax=p0.975, 
                                              y = 1))+ 
   theme(
         axis.text.y=element_blank(),
         axis.ticks.y=element_blank())
 
 
  
 key_subset <- trace_to_sample_m[variable=="liklihood" |
                                   variable == "sig" |
                                   variable =="rho", ]
 key_subset$variable <- factor(key_subset$variable, levels = c("sig", "liklihood","rho"))
key_subset[variable == "liklihood", variable := "Likelihood"]
key_subset[variable == "sig", variable := "Interaction parameter (sigma)"]
key_subset[variable == "rho", variable := "Duration of cross-protection (rho)"]

 
 SUSBET <- ggplot(key_subset, aes(x = value)) + 
    geom_histogram(bins =30, position = "stack", aes(fill = section)) + 
    #geom_density(aes(fill = section) )+
   facet_wrap(variable~., scales = "free", ncol = 1) + 
   labs(x = "Parameter value", y = "Count", title = "B: Posteriors by mode", 
        fill = "Sigma subset") + 
   theme_linedraw() + 
   scale_fill_manual(values = c("slategray4", "skyblue3")) +
   theme(strip.background = element_rect(colour="white", fill="white"),
         strip.text = element_text(colour = 'black',hjust = 0),
         panel.spacing.x=unit(0.1, "lines"),panel.spacing.y=unit(-0, "lines"))
 
 
 tiff(here("figures","All_interaction.tiff"), height = 2000, width = 3200, res = 300)
 ALL
 dev.off()
 

###### FIGURE 4: Simulations no covid #######
 VACC <- quantiles_noflu_rel(trace_to_sample[sig>0.2,], multi_data, 
                             sample_size = 1000, together=T,
                             rolling_mean = 1,  
                             from=1, to=15953)
 
 
 
 tiff(here::here("Figures","Fig2a_fit.tiff"), height = 2000, width = 3200, res = 300)
 
 grid.arrange(FITS_TOGETHER, SUSBET,
              layout_matrix = rbind(c(1,1,1,1,2,2),
                                    c(1,1,1,1,2,2),
                                    c(1,1,1,1,2,2)))
 
 dev.off()
 
 tiff(here::here("Figures","Fig2b_fit.tiff"), height = 1000, width = 3500, res = 300)
 
 grid.arrange(VOID, CASE, VACC,
              layout_matrix = rbind(c(1,2,3),
                                    c(1,2,3)))
 
 dev.off()
 
 
 tiff(here::here("Figures","Fig2a_fit.tiff"), height = 2000, width = 3200, res = 300)
 
 grid.arrange(FITS_TOGETHER, SUSBET,
              layout_matrix = rbind(c(1,1,1,1,2,2),
                                    c(1,1,1,1,2,2),
                                    c(1,1,1,1,2,2)))
 
 dev.off()
 
 tiff(here::here("Figures","Fig2_fit.tiff"), height = 3000, width = 3500, res = 300)
 
 grid.arrange(FITS_TOGETHER, SUSBET, VOID, CASE, VACC,
              layout_matrix = rbind(c(1,1,1,1,2,2),
                                    c(1,1,1,1,2,2),
                                    c(1,1,1,1,2,2),
                                    c(4,4,5,5,6,6),
                                    c(4,4,5,5,6,6)))
 
 dev.off()
 
 
 
 ###### POSTERIORS ########
 # Posteriors for text
 quants <- c(0.025,0.5,0.975)
 quantile(trace_to_sample$betaR, probs = quants)*R0ratR
 quantile(trace_to_sample$betaI, probs = quants)*R0ratI
 
 quantile(trace_to_sample$sig, probs = quants) 
 quantile(trace_to_sample[sig>0.2,]$sig, probs = quants) 
 quantile(trace_to_sample[sig<=0.2,]$sig, probs = quants) 
 
 1/quantile(trace_to_sample[sig>0.2,]$rho, probs = quants)
 1/quantile(trace_to_sample$rho, probs = quants) 
 
 quantile(trace_to_sample[sig>0.2,]$sig, probs = quants) 
 
 # Probability of reporting dual
 trace_to_sample_sub <- trace_to_sample[sig>0.2,]
 #P(reporting| flu + rsv)
trace_to_sample_sub[,"rep_dual_0"] <- (trace_to_sample_sub$Dual_mult*trace_to_sample_sub$Rdetect2) /
# P(reporting|flu) + P(reporting|RSV)
  (trace_to_sample_sub$Rdetect2 + trace_to_sample_sub$Idetect)

quantile(trace_to_sample_sub$rep_dual_0, probs = quants)
 
 #P(reporting| flu + rsv)
trace_to_sample_sub[,"rep_dual_1"] <- (trace_to_sample_sub$Dual_mult*trace_to_sample_sub$Rdetect5) /
   # P(reporting|flu) + P(reporting|RSV)
   (trace_to_sample_sub$Rdetect5 + (trace_to_sample_sub$Idetect*
      trace_to_sample_sub$detect_multiplier))  
 
quantile(trace_to_sample_sub$rep_dual_1, probs = quants)
 
 # Probability of reporting dual
 trace_to_sample_sub_other <- trace_to_sample[sig<=0.2,]
 #P(reporting| flu + rsv)

   trace_to_sample_sub_other[,"rep_dual_0"] <- (trace_to_sample_sub_other$Dual_mult*trace_to_sample_sub_other$Rdetect2) /
   # P(reporting|flu) + P(reporting|RSV)
   (trace_to_sample_sub_other$Rdetect2 + trace_to_sample_sub_other$Idetect)
 
 quantile(trace_to_sample_sub_other$rep_dual_0, probs = quants)
 
 #P(reporting| flu + rsv)
 trace_to_sample_sub_other[,"rep_dual_1"] <- (trace_to_sample_sub_other$Dual_mult*trace_to_sample_sub_other$Rdetect5) /
   # P(reporting|flu) + P(reporting|RSV)
   (trace_to_sample_sub_other$Rdetect5 + (trace_to_sample_sub_other$Idetect*
                                            trace_to_sample_sub_other$detect_multiplier))  
 
 quantile(trace_to_sample_sub_other$rep_dual_1, probs = quants)
 