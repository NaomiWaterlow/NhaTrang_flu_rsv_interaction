# Extra plots following review.

source(here("revision_extras","reviewer_plots_functions.R"))

# 1. Dynamics in the other ages. 
# the input of the posteriors is "trace_edited[[1]]"

ALL_AGES_PLOT_PROP <- fit_all_ages_spaghetti(trace_edited[[1]], multi_data, 
                        sample_size = 100, together=T,
                        rolling_mean = 1,  
                        from=1, to=24651)

# 2. Flu -> RSV sigma make any difference?

trace_subset_mode0 <- trace_edited[[1]][sig < 0.2]

FITS_MODE_0 <-fit_quantiles_together(trace_subset_mode0 , multi_data, 
                                       sample_size = 5, together=T,
                                       rolling_mean = 1,  
                                       from=1, to=nrow(trace_subset_mode0))

trace_subset_mode1 <- trace_edited[[1]][sig >= 0.2]

FITS_MODE_1 <-fit_quantiles_together(trace_subset_mode1 , multi_data, 
                                     sample_size = 5, together=T,
                                     rolling_mean = 1,  
                                     from=1, to=nrow(trace_subset_mode0))

tiff(here("Figures","ALL_AGES_PLOT.tiff"), height = 2000, width = 3200, res = 300)
ALL_AGES_PLOT
dev.off()

tiff(here("Figures","ALL_AGES_PLOT_PROP.tiff"), height = 2000, width = 3200, res = 300)
ALL_AGES_PLOT_PROP
dev.off()

tiff(here("Figures","FITS_MODE_0.tiff"), height = 2000, width = 3200, res = 300)
FITS_MODE_0
dev.off()

tiff(here("Figures","FITS_MODE_1.tiff"), height = 2000, width = 3200, res = 300)
FITS_MODE_1
dev.off()

compileModel(here::here("sigma_oneway.cpp"), "build2")
dyn.load("build2/sigma_oneway.so")
source(here::here("visualise fit","quantile_plot_functions.R"))
source(here::here("visualise fit","quantile_plots_dual.R"))

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
    overdispersion = unname(param_values[,"overdispersion"]), 
    RSV_mult = unname(param_values[,"RSV_mult"]), 
    Dual_mult = unname(param_values[,"Dual_mult"]), 
    sig_fixed = 1
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
  ][, Dual_0 := Dual_incidence_1*Parameters[["Rdetect2"]]*Parameters[["RSV_mult"]] *Parameters[["Dual_mult"]]
  ][, Dual_1 := Dual_incidence_2*Parameters[["Rdetect5"]]*Parameters[["RSV_mult"]]*Parameters[["Dual_mult"]]]
  return(trajectory)
  
}

## Simulate the SIPR model, with detections
simulate_SIPR <- function (theta_season, init_state, times)
{
  outall <- data.table(ode(y = init_state, 
                           t = times, 
                           initfunc = "initmod",
                           dllname = "sigma_oneway",
                           func = "derivatives", 
                           parms = theta_season, 
                           method = "ode23"))
  
  #Calculate the incidence
  trajectory <- summary_stats(outall)

  return(trajectory)
}



FITS_ONEWAY_0 <-fit_quantiles_together(trace_edited[[1]], multi_data, 
                                       sample_size = 5, together=T,
                                       rolling_mean = 1,  
                                       from=700, to=7500)
# Change the value above in the run simulations quantiles function
FITS_ONEWAY_1 <-fit_quantiles_together(trace_edited[[1]], multi_data, 
                                       sample_size = 5, together=T,
                                       rolling_mean = 1,  
                                       from=700, to=7500)

tiff(here("Figures","FITS_ONEWAY_1.tiff"), height = 2000, width = 3200, res = 300)
FITS_ONEWAY_1
dev.off()

tiff(here("Figures","FITS_ONEWAY_0.tiff"), height = 2000, width = 3200, res = 300)
FITS_ONEWAY_0
dev.off()
