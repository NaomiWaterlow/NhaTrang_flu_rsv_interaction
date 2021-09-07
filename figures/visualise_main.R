
load(here("AWS_out","AWS_mid_detect_5_2020-07-30.Rdata"))

# get all the chains but not the temperature
trace <- total_trace[1:(length(total_trace))]
trace <- lapply(trace, function(i) na.omit(i))

######### TRACES COMBINED ANALYSIS ###########

trace_temps <- trace[[length(trace)]]
#temperatures <- 1/trace_temps[dim(trace_temps)[1],]
temperatures <- rev(c(1:12))
# add the timesteps to the trace outputs
timesteps <- 1:dim(trace[[1]])[1]
trace <- lapply(trace, function(i) label_trace(i, detectI = T))

trace_edited <- lapply(trace, function(i) format_trace(i, thin=10, burnin=250000,
                                                       keep_all = F))

# transform the parameters back to the origional values
trace_edited <- lapply(trace_edited, function(i) transform_back(i))

trace_edited <- Map(cbind, trace_edited, temp = temperatures)

trace_combined <- do.call(rbind.data.frame, trace_edited)

trace_combined <- melt(trace_combined, id=c("time_steps", "temp"))
trace_combined$temp <- 1/trace_combined$temp

trace_combined$temp <- factor(trace_combined$temp)
trace_combined$temp <- fct_rev(trace_combined$temp)

nb.cols <- length(unique(trace_combined$temp))
mycolors <- colorRampPalette(brewer.pal(8, "YlGnBu"))(nb.cols)

DENS_PLOTS_TRANS <- list()
TRACE_PLOTS_TRANS <- list()

trace_combined <- as.data.frame(trace_combined)
for (i in c(1:nlevels(trace_combined$variable))){
  
  parm <- levels(trace_combined$variable)[i]
  parm <- as.character(parm)
  #denisty plo
  DENS_PLOTS_TRANS[[parm]]<-(ggplot(subset(
    trace_combined
    , variable==parm), aes(x=value)) +
      geom_density()+
      ggtitle(parm)+
      theme(plot.title = element_text(face = "bold"))
  )
  # trace plot
  TRACE_PLOTS_TRANS[[parm]] <- (ggplot(subset(
    trace_combined
    , variable==parm)) +
      geom_path(aes(y=value, x=time_steps, colour = temp)) +
      labs(title =parm)+
      #  scale_color_manual(values=mycolors) + theme_bw()+
      theme(plot.title = element_text(face = "bold"),
            legend.position = "none")
  )
}



do.call("grid.arrange", c(TRACE_PLOTS_SINGLE[c( "liklihood",#"betaR", "betaI", "Rdetect2", "Rdetect5",
                                               "Idetect", #"propR_1", "propI_1", 
                                               "sig", "rho", 
                                               "inf_sus_1", "overdispersion")], ncol=1))

tiff(here("Figures","FIGSX2_TRACE_MAIN"), height = 2000, width = 3200, res = 300)
do.call("grid.arrange", c(TRACE_PLOTS_TRANS[c( "liklihood","betaR", "betaI", "Rdetect5",
                                                "propR_1", "sig", "rho",
                                                "inf_sus_1","overdispersion", "Dual_mult" )],ncol=2))
do.call("grid.arrange", c(TRACE_PLOTS_SINGLE[c( "liklihood","betaR", "betaI", "Rdetect5",
                                                "propR_1", "sig", "rho",
                                                "inf_sus_1","overdispersion", "Dual_mult" )], ncol=2))
dev.off()

trace_temps <- as.data.frame(trace_temps)
trace_temps$step <- c(1:dim(trace_temps)[1])
trace_temps_2 <- melt(trace_temps, id="step")

TEMPS <- ggplot(trace_temps_2, aes(x=step, y=value, colour = variable)) +
  geom_line()#+ ylim(c(0,500))
TEMPS

# calculate acceptance rate of parameters in chain 1
AR <- sum(trace[[1]][,"Step_accepted"])/dim(trace[[1]])[1]
AR
# calculate the number of swaps done by chain 1
swap_rate <- sapply(trace, function(i) swaps_accepted(i))
swap_rate



trace_single <- trace[[1]]

# format the trace
colnames(trace_single) <- c("liklihood", names(init_theta), "Step_accepted","Swap_accepted_to",
                            "swap_attempted")
trace_single <- format_trace(trace_single, thin=1, burnin = 1, keep_all = F)
trace_single <- transform_back(trace_single)
trace_single_mt <- melt(trace_single[,], id=c("time_steps"))

DENS_PLOTS_SINGLE <- list()
TRACE_PLOTS_SINGLE <- list()

for (i in c(1:nlevels(trace_single_mt$variable))){
  
  parm <- levels(trace_single_mt$variable)[i]
  parm <- as.character(parm)
  #denisty plo
  DENS_PLOTS_SINGLE[[parm]]<-(ggplot(subset(
    trace_single_mt
    , variable==parm), aes(x=value)) +
      geom_density()+
      ggtitle(parm)+
      theme(plot.title = element_text(face = "bold"))
  )
  # trace plot
  TRACE_PLOTS_SINGLE[[parm]] <- (ggplot(subset(
    trace_single_mt
    , variable==parm)) +
      geom_path(aes(y=value, x=time_steps)) +
      labs(title =parm)+
      # scale_color_manual(values=mycolors) + theme_bw()+
      theme(plot.title = element_text(face = "bold"),
            legend.position = "none")
  )
}

do.call("grid.arrange", c(TRACE_PLOTS_SINGLE, ncol=1))

# # plot the quantiles from the main trace
set.seed(100)
FITS_TOGETHER <-fit_quantiles_together(trace_edited[[1]], multi_data, 
                                             sample_size = 5, together=T,
                                             rolling_mean = 1,  
                                             from=700, to=7500)
FITS_TOGETHER


####### Likelihoods ######
lik_store <- data.table()
trace_tosamp <- trace_edited[[1]]
trace_tosamp$RSV_default <- 1
for(j in c(1:2)){
  sample_num <- sample(x=1:nrow(trace_tosamp), size = 1)
input_params <- unlist(trace_tosamp[sample_num,])
input_params <- as.list(input_params)
# Likelihood by season
 for (season in c(1:11)){
#subset parameters
  
pars_season <- subset_parameters(input_params, season)
data_season <- multi_data[[season]]

init_state_proposed <- calculate_init_state(population_numbers,
                                            pars_season)
likelihood <- calc_likelihood_overall(pars_season,
                                               init_state_proposed,
                                               data_season,
                                               log = T)

lik_store <- rbind(lik_store, data.frame(j, season, likelihood))
}}
lik_store$season <- as.factor(lik_store$season)
ggplot(lik_store, aes(x = season, y = likelihood)) + 
  geom_point(alpha=0.5)

# pairs plot
PAIRS_PLOT <- ggpairs(as.data.frame(trace_single[seq(1,dim(trace_single)[1], by=10),
                                                 c(2:45)]))
PAIRS_PLOT

# create a correlation matrix
covmat_out <- cov(trace[[1]][,c(2:ncol(trace[[1]]))])
saveRDS(covmat_out[1:43,1:43], "20210628_multiinf_comb.RDS", version=2)
covmat_out <- cov(traces_cov[,2:43])
covmat <- readRDS("20210518_nos9.RDS")
covmat[,""]

cov_plot <- melt(covmat_out)
ggplot(cov_plot, aes(x = Var1, y = Var2, fill = abs(value))) + 
  geom_tile()

cor_mat <- cor(trace_edited[[1]][,2:45])

cor_mat <- data.frame(cor_mat)
cor_mat$rowname <- rownames(cor_mat)
cor_mat_2 <- melt(cor_mat, id = "rowname")
cor_mat_2$rowname <- factor(cor_mat_2$rowname, levels = levels(cor_mat_2$variable))
 
ggplot(cor_mat_2, aes(x = rowname, y=variable, fill = value)) + 
  geom_tile() + 
  scale_fill_distiller(palette = "Spectral")

ggplot(cor_mat_2, aes(x = rowname, y=variable, fill = abs(value))) + 
  geom_tile() #+ 
#scale_fill_distiller(palette = "Spectral")

cor_mat_2 <- data.table(cor_mat_2)
cor_mat_2[abs(value) < 0.85,cutoff := 0 ]
cor_mat_2[abs(value) >= 0.85,cutoff := 1 ]

ggplot(cor_mat_2, aes(x = rowname, y=variable, fill = cutoff)) + 
  geom_tile()


trace_1 <- trace_combined
trace_1$run <- 1
trace_2 <- trace_combined
trace_2$run <- 2
traces_together <- data.table(rbind(trace_1, trace_2))
traces_plotter <- traces_together[(temp == 1)]
traces_plotter <- traces_plotter[time_steps > 43000]

traces_plotter$run <- as.factor(traces_plotter$run)
#traces_plotter$run <-fct_rev(traces_plotter$run)
ggplot(traces_plotter) + 
  geom_density(aes(x = value, fill =run)) + 
  facet_wrap(variable~., scales = "free", ncol=5)


######## Combined plots

testing <- as.data.table(trace_1[[1]])[10000:32500,]
testing[sig > 0.25, set := 1]
testing[sig <= 0.25, set := 2]

testing_m <- melt(testing, id.vars= c("set") )

ggplot(testing_m, aes(x = value, group=set,
                    fill =set, colour  = set)) + 
  geom_density() + 
  facet_wrap(variable~., scales="free")

####### Combine traces

load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/multinom_infections_covmat_17_2021-07-01.Rdata")
trace <- total_trace[1:(length(total_trace))-1]
trace <- lapply(trace, function(i) na.omit(i))
first_trace <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/multinomial_infections_cont_17_2021-07-05.Rdata")
trace <- total_trace[1:(length(total_trace))]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/multinomial_infections_cont2_17_2021-07-09.Rdata")
trace <- total_trace[1:(length(total_trace))]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace1 <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/multinomial_infections_cont3_17_2021-07-12.Rdata")
trace <- total_trace[1:(length(total_trace))]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace2 <- trace


total_trace <- lapply(seq_along(first_trace), function(x) rbind(first_trace[[x]], addon_trace[[x]]))
total_trace <- lapply(seq_along(total_trace), function(x) rbind(total_trace[[x]], addon_trace1[[x]]))
total_trace <- lapply(seq_along(total_trace), function(x) rbind(total_trace[[x]], addon_trace2[[x]]))

name <- "combined_trace_17"
save(total_trace, file=paste0(name,"_", Sys.Date() , ".Rdata"))




load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/multinomial_infections_covmat_16_2021-07-01.Rdata")
trace <- total_trace[1:(length(total_trace))-1]
trace <- lapply(trace, function(i) na.omit(i))
first_trace <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/multinomial_infections_cont_16_2021-07-05.Rdata")
trace <- total_trace[1:(length(total_trace))]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/multinomial_infections_cont2_16_2021-07-09.Rdata")
trace <- total_trace[1:(length(total_trace))]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace1 <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/multinomial_infections_cont3_16_2021-07-12.Rdata")
trace <- total_trace[1:(length(total_trace))]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace2 <- trace


total_trace <- lapply(seq_along(first_trace), function(x) rbind(first_trace[[x]], addon_trace[[x]]))
total_trace <- lapply(seq_along(total_trace), function(x) rbind(total_trace[[x]], addon_trace1[[x]]))
total_trace <- lapply(seq_along(total_trace), function(x) rbind(total_trace[[x]], addon_trace2[[x]]))

name <- "combined_trace_16"
save(total_trace, file=paste0(name,"_", Sys.Date() , ".Rdata"))


load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_standard_17_2021-08-01.RData")
trace <- total_trace[1:(length(total_trace))-1]
trace <- lapply(trace, function(i) na.omit(i))
first_trace <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_standard_cont_17_2021-08-08.RData")
trace <- total_trace[1:(length(total_trace))-1]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_standard_cont2_17_2021-08-15.RData")
trace <- total_trace[1:(length(total_trace))-1]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace2 <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_standard_cont3_17_2021-08-23.RData")
trace <- total_trace[1:(length(total_trace))-1]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace3 <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_standard_cont4_17_2021-08-25.RData")
trace <- total_trace[1:(length(total_trace))]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace4 <- trace

total_trace <- lapply(seq_along(first_trace), function(x) rbind(first_trace[[x]], addon_trace[[x]]))
total_trace <- lapply(seq_along(first_trace), function(x) rbind(total_trace[[x]], addon_trace2[[x]]))
total_trace <- lapply(seq_along(first_trace), function(x) rbind(total_trace[[x]], addon_trace3[[x]]))
total_trace <- lapply(seq_along(first_trace), function(x) rbind(total_trace[[x]], addon_trace4[[x]]))

name <- "combined_standard_17"
save(total_trace, file=paste0("~/Documents/GitHub/NhaTrang_fitting/AWS_out/",name,"_", Sys.Date() , ".Rdata"))




load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_dualmult_16_2021-08-01.RData")
trace <- total_trace[1:(length(total_trace))]
trace <- lapply(trace, function(i) na.omit(i))
first_trace <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_dualmult_cont_16_2021-08-08.RData")
trace <- total_trace[1:(length(total_trace))-1]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_dualmult_cont2_16_2021-08-15.RData")
trace <- total_trace[1:(length(total_trace))-1]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace2 <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_dualmult_cont3_16_2021-08-24.RData")
trace <- total_trace[1:(length(total_trace))-1]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace3 <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_dualmult_cont4_16_2021-08-31.RData")
trace <- total_trace[1:(length(total_trace))-1]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace4 <- trace

total_trace <- lapply(seq_along(first_trace), function(x) rbind(first_trace[[x]], addon_trace[[x]]))
total_trace <- lapply(seq_along(first_trace), function(x) rbind(total_trace[[x]], addon_trace2[[x]]))
total_trace <- lapply(seq_along(first_trace), function(x) rbind(total_trace[[x]], addon_trace3[[x]]))
total_trace <- lapply(seq_along(first_trace), function(x) rbind(total_trace[[x]], addon_trace4[[x]]))

name <- "combined_dualmult_16"
save(total_trace, file=paste0("~/Documents/GitHub/NhaTrang_fitting/AWS_out/",name,"_", Sys.Date() , ".Rdata"))


load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_dualmult_17_2021-08-14.RData")
trace <- total_trace[1:(length(total_trace))-1]
trace <- lapply(trace, function(i) na.omit(i))
first_trace <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_dualmult_cont_17_2021-08-28.RData")
trace <- total_trace[1:(length(total_trace))-1]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_dualmult_cont2_17_2021-08-31.RData")
trace <- total_trace[1:(length(total_trace))]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace2<- trace

total_trace <- lapply(seq_along(first_trace), function(x) rbind(first_trace[[x]], addon_trace[[x]]))
total_trace <- lapply(seq_along(first_trace), function(x) rbind(total_trace[[x]], addon_trace2[[x]]))
name <- "combined_dualmult_17"
save(total_trace, file=paste0("~/Documents/GitHub/NhaTrang_fitting/AWS_out/",name,"_", Sys.Date() , ".Rdata"))



load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_interactionprior_16_2021-08-01.RData")
trace <- total_trace[1:(length(total_trace))-1]
trace <- lapply(trace, function(i) na.omit(i))
first_trace <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_interaction_cont_16_2021-08-08.RData")
trace <- total_trace[1:(length(total_trace))-1]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_interaction_cont2_16_2021-08-15.RData")
trace <- total_trace[1:(length(total_trace))-1]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace2 <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_interaction_cont3_16_2021-08-23.RData")
trace <- total_trace[1:(length(total_trace))-1]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace3 <- trace
load("~/Documents/GitHub/NhaTrang_fitting/AWS_out/Pop_interaction_cont4_16_2021-08-25.RData")
trace <- total_trace[1:(length(total_trace))]
trace <- lapply(trace, function(i) na.omit(i))
addon_trace4 <- trace

total_trace <- lapply(seq_along(first_trace), function(x) rbind(first_trace[[x]], addon_trace[[x]]))
total_trace <- lapply(seq_along(first_trace), function(x) rbind(total_trace[[x]], addon_trace2[[x]]))
total_trace <- lapply(seq_along(first_trace), function(x) rbind(total_trace[[x]], addon_trace3[[x]]))
total_trace <- lapply(seq_along(first_trace), function(x) rbind(total_trace[[x]], addon_trace4[[x]]))

name <- "combined_interaction_16"
save(total_trace, file=paste0("~/Documents/GitHub/NhaTrang_fitting/AWS_out/",name,"_", Sys.Date() , ".Rdata"))


addon_trace4

addon_trace4 <- lapply(addon_trace4, function(x) cut(x))

cut <- function(x){
  x <- x[1:50000,]
  return(x)
}


