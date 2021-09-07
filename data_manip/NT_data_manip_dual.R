#Manipulate the data!

##### Dates to cut #####

# These were used as initial tests
# dates_to_cut <- c(as.Date("2007-02-07"),
#                   as.Date("2008-01-28"),
#                   as.Date("2009-03-16"),
#                   as.Date("2010-03-15"),
#                   as.Date("2011-05-09"),
#                   as.Date("2012-02-27"),
#                   as.Date("2013-01-28"),
#                   as.Date("2014-01-06"),
#                   as.Date("2015-02-02"),
#                   as.Date("2016-03-28"),
#                   as.Date("2017-04-17"), 
#                   as.Date("2018-01-01"))

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

##### Load and manipulate#####

#read in the data
setwd("/Users/lsh1402815/Documents/Data/NhaTrang")
case_table <- as.data.table(read.csv("infa_rsv_age.csv"))

# Calculate average age of cases ######

mean(case_table[(ageyear < 5) & 
           (infa ==1),]$ageyear)*12
mean(case_table[(ageyear < 5) & 
                  (rsv ==1),]$ageyear)*12

#remove case with no date (formid == 4482)
case_table[npsdate==""] <- NA
# split into age groups, remove older ages and remove unneccesary column
case_table[ageyear==0 | ageyear==1, age_group := 1]
case_table[ageyear==2 | ageyear==3 | ageyear==4, age_group := 2]
case_table <- case_table[!is.na(age_group), ]
case_table[,sex:=NULL]
case_table[,formid:=NULL]
case_table[,ageyear:=NULL]
case_table[,both := 0]
# add the dual cases to the both column
case_table[infa > 0 & rsv > 0 ,both := 1]
# remove the dual cases from the individual columns
case_table[both ==1 , infa := 0]
case_table[both ==1 & rsv > 0 , rsv := 0]
#sum for each day rather than individual cases
case_table_summed <- case_table[, lapply(.SD, sum), by=list(npsdate, age_group)]
# make date format
case_table_summed$npsdate <- as.Date(as.character(case_table_summed$npsdate),
                                               format="%d-%b-%y")
#rearrange to desired column titles (e.g. Hos_RSV_0)
case_table_summed[which(infa>0 & age_group==1), "Hos_INF_0" :=infa]
case_table_summed[which(infa>0 & age_group==2), "Hos_INF_1" :=infa]
case_table_summed[which(rsv>0 & age_group==1), "Hos_RSV_0" :=rsv]
case_table_summed[which(rsv>0 & age_group==2), "Hos_RSV_1" :=rsv]
case_table_summed[which(both>0 & age_group==1), "both_0" :=both]
case_table_summed[which(both>0 & age_group==2), "both_1" :=both]

case_table_summed[,c("age_group", "infa", "rsv", "both") :=NULL]
# convert NAs to 0s. 
case_table_summed[is.na(case_table_summed)]<-0

ctm <- melt.data.table(case_table_summed, id.vars = c("npsdate"))
ctc <- dcast(ctm, formula = npsdate ~ variable, fun.aggregate = sum)

# Create an empty incidence table to fill in gaps
empty_incidence <- data.table(npsdate = as.Date(c(1:3962), origin = "2007-01-30"), 
                              Hos_INF_0 = as.numeric(0), 
                              Hos_INF_1 = as.numeric(0), 
                              Hos_RSV_0 = as.numeric(0), 
                              Hos_RSV_1 = as.numeric(0), 
                              both_0 = as.numeric(0), 
                              both_1 = as.numeric(0))

#replace the 0 values with real values where exist
empty_incidence[npsdate %in% ctc$npsdate, Hos_INF_0 := unlist(ctc[,"Hos_INF_0"])]
empty_incidence[npsdate %in% ctc$npsdate, Hos_INF_1 := unlist(ctc[,"Hos_INF_1"])]
empty_incidence[npsdate %in% ctc$npsdate, Hos_RSV_0 := unlist(ctc[,"Hos_RSV_0"])]
empty_incidence[npsdate %in% ctc$npsdate, Hos_RSV_1 := unlist(ctc[,"Hos_RSV_1"])]
empty_incidence[npsdate %in% ctc$npsdate, both_0 := unlist(ctc[,"both_0"])]
empty_incidence[npsdate %in% ctc$npsdate, both_1 := unlist(ctc[,"both_1"])]
dat_c <- empty_incidence

# convert to week / year
dat_c[, week := isoweek(npsdate)]
dat_c[, year := isoyear(npsdate)]

start_weeks <- dat_c[,min(npsdate), by =c("week", "year")]

for(i in 1:nrow(start_weeks)){
  target_week <- unlist(start_weeks[i,"week"])
  target_year <- unlist(start_weeks[i,"year"])
  target_date <- as.Date(unlist(start_weeks[i,"V1"]))
  dat_c[week == target_week & year == target_year, start_week := 
         as.character(target_date) ]
}
# label the sections/ seasons
for(i in 1:(length(dates_to_cut)-1)){
dat_c[((npsdate >= dates_to_cut[i]) & (npsdate < dates_to_cut[i+1])),section := i]
  }
#format
dat_c[,c("npsdate") := NULL]
dat_cm <- melt.data.table(dat_c, id.vars = c("section","week", "year", "start_week"))
dat_cmc <- na.omit(dcast.data.table(dat_cm,
                            formula = start_week + section + week + year~variable,
                            fun.aggregate = sum))

# Gove a quick plot to test
plot_m <- melt.data.table(dat_cmc, id.vars = c("start_week", "section", "week", "year"))
plot_m$start_week <- as.Date(plot_m$start_week)
ggplot(plot_m, aes(x = start_week, y = value )) + 
  geom_line()+ 
  facet_grid(variable~.) +
  geom_vline(xintercept = dates_to_cut, colour = "darkred", alpha = 0.5) 

#Split into individual seasons
  y_1 <- as.data.frame(dat_cmc[section ==1,])
  y_2 <- as.data.frame(dat_cmc[section ==2,])
  y_3 <- as.data.frame(dat_cmc[section ==3,])
  y_4 <- as.data.frame(dat_cmc[section ==4,])
  y_5 <- as.data.frame(dat_cmc[section ==5,])
  y_6 <- as.data.frame(dat_cmc[section ==6,])
  y_7 <- as.data.frame(dat_cmc[section ==7,])
  y_8 <- as.data.frame(dat_cmc[section ==8,])
  y_9 <- as.data.frame(dat_cmc[section ==9,])
  y_10 <- as.data.frame(dat_cmc[section ==10,])
  y_11 <- as.data.frame(dat_cmc[section ==11,])

  multi_data <- list(observed1 = y_1,
                     observed2 = y_2,
                     observed3 = y_3,
                     observed4 = y_4,
                     observed5 = y_5,
                     observed6 = y_6,
                     observed7 = y_7,
                     observed8 = y_8, 
                     observed9 = y_9,
                     observed10 = y_10,
                     observed11 = y_11)
# Save the file  
# saveRDS(multi_data, file = "multi_data_dual_JULY.RDS", version = 2)
  
  ######### Work out when most cases occur ######
  when <- copy(dat_cmc)
  when[, influenza := Hos_INF_0 + Hos_INF_1 + both_0 + both_1]
  when[, rsv := Hos_RSV_0 + Hos_RSV_1 + both_0 + both_1]
  when[,c("Hos_INF_0", "Hos_INF_1", "Hos_RSV_0", "Hos_RSV_1", "both_0", "both_1") := NULL]
  when[, sum(influenza), by = week]
  
  when$inf_roll <- rollmean(when$influenza, k=4, fill ="NA")
  when$rsv_roll <- rollmean(when$rsv, k=4, fill ="NA")
  when[, sum(influenza), by=year]
  
  when[,year_sum_inf := sum(influenza), by=year]
  when[,year_sum_rsv := sum(rsv), by=year]
  
 
  
  both_when <- when[, sum(rsv), by = week]
  both_when$influenza <- when[, sum(influenza), by = week]$V1
  colnames(both_when) <- c("Week", "RSV", "influenza")
  both_when$RSV <- (both_when$RSV / (rsv_0 +rsv_1 + dual_0 + dual_1))*100
  both_when$influenza <- (both_when$influenza / (inf_0 +inf_1 + dual_0 + dual_1))*100
  both_when_m <- melt(both_when, id.vars=c("Week"))
  
  
  when[, influenza := (inf_roll/year_sum_inf)*100]
  when[, RSV := (rsv_roll/year_sum_rsv)*100]
  when_sub <- copy(when[,c("week", "year", "influenza", "RSV")])
  when_sub_m <- melt.data.table(when_sub, id.vars = c("week", "year"))
  
  
  when_sub_m <- as.data.frame(when_sub_m)
  SEASONALITY <- ggplot(both_when_m, aes(x = Week, y = value, colour = variable)) + 
    geom_line(data = when_sub_m, aes(x = week, y = value, group=interaction(variable, year)), alpha = 0.4) +
    theme(legend.position = "none")+
     geom_line(size = 1.5) +
    theme_linedraw() +
    scale_colour_manual(values = c("#fc8d62","#66c2a5")) + 
    labs(title = "C: Seasonality",x = "Week of the year", y = "Percentage of cases reported", colour = "Virus")  
  SEASONALITY
  
# Look at total cases by year (NOT season)
total_year <- data.table(matrix(nrow =11, ncol = 5))
total_year[,1:2] <-dat_cmc[,sum(Hos_INF_0), by=year]
total_year[,3] <-dat_cmc[,sum(Hos_INF_1), by=year][,2] 
total_year[,4] <-dat_cmc[,sum(Hos_RSV_0), by=year][,2] 
total_year[,5] <-dat_cmc[,sum(Hos_RSV_1), by=year][,2] 
colnames(total_year) <- c("year", "flu_infants","flu_children",
                          "rsv_infants", "rsv_children")  
  
temp <- melt.data.table(total_year, id.vars = "year")  
ggplot(temp, aes(x = year, y = value, colour = variable)) +
  geom_line()

##### CALCulate total number of reports by age group / virus #####
dat_cmc_m <- melt.data.table(dat_cmc, id.vars = c("start_week", 
                                                  "section", 
                                                  "week", 
                                                  "year"))

summary_tab <- dat_cmc_m[, sum(value), by = variable]
summary_tab[,virus := c("Influenza", "Influenza", "RSV", "RSV", "Both", "Both")]
summary_tab[,age_group := c("0-1", "2-4","0-1", "2-4","0-1", "2-4")]
summary_tab$virus <- factor(summary_tab$virus, 
                            levels = c("RSV", "Influenza","Both"))
cc <- c("#66c2a5","#fc8d62", "#8da0cb")
TOTAL_REPORTS <-ggplot(summary_tab, aes(x = age_group, y = V1, fill = virus)) + 
  geom_bar(stat= "identity", position = "dodge") + 
  theme_minimal()  + 
  scale_fill_manual(values= cc) + 
  labs(x = "Age group", y = "Count", fill = "Virus", 
       title = "B: Total cases reported by age and infection")

##### correlations ####
nice_names <- dat_cmc
nice_names$flu <- nice_names$Hos_INF_0 + nice_names$Hos_INF_1 + 
  nice_names$both_0 + nice_names$both_1
nice_names$rsv <- nice_names$Hos_RSV_0 + nice_names$Hos_RSV_1 + 
  nice_names$both_0 + nice_names$both_1
colnames(nice_names) <- c("start_week", "section", "week", "year", 
                          "Influenza 0-1", "Influenza 2-4", 
                          "RSV 0-1", "RSV 2-4", "Both 0-1", 
                          "Both 2-4", "influenza", "RSV")
cor_tab <- as.data.frame(cor(nice_names[,5:10]))

cor_tab$type <- rownames(cor_tab)
cor_tab_m <- melt.data.table(as.data.table(cor_tab), id.vars = "type")
cor_tab_m$type <- factor(cor_tab_m$type, levels = levels(cor_tab_m$variable))
CORRELATION_PLOT <- ggplot(cor_tab_m, aes(x = variable, y = type)) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient2(low = "#f1a340", mid = "#f7f7f7", high = "#998ec3", 
                       midpoint = 0) + 
  theme_minimal() + 
  labs(x = "", y = "", fill = "Correlation", 
       title = "C: Correlation between time series") + 
  theme(axis.text.x = element_text(angle =45, 
                                   hjust = 1))
CORRELATION_PLOT

# do a correlation test on flu vs RSV
cor.test(x=nice_names$influenza, y = nice_names$RSV)


SCATTER_COR <- ggplot(nice_names, aes(x = influenza, y = RSV)) + geom_point(alpha = 0.4) + 
  theme_linedraw() + 
  labs(x = "Weekly influenza cases reported", 
       y = "Weekly RSV cases reported")

tiff(here::here("Figures","SUP_scatter"), height = 2000, width = 2000, res = 300)
SCATTER_COR
dev.off()


#### Create the time series plot #####

dat_cmc[, Influenza := Hos_INF_0 + Hos_INF_1 + both_0 +both_1]
dat_cmc[,RSV := Hos_RSV_0 + Hos_RSV_1 + both_0 +both_1]

last_graph <- dat_cmc[,c(1,11,12)]
last_graph_m <- melt.data.table(last_graph, id.vars= c("start_week"))
last_graph_m$start_week <- as.Date(last_graph_m$start_week)
last_graph_m$variable <- factor(last_graph_m$variable, levels = c("RSV", "Influenza"))
TIMESERIES <- ggplot(last_graph_m, aes(x = start_week, y = value, colour = variable)) + 
  geom_line() +# facet_grid(variable~.) + 
  theme_linedraw() + 
  labs(x = "Time", y = "Cases Reported", colour = "Virus", 
       title = "A: Weekly time series of Influenza and RSV reported") + 
  theme(legend.position = "NULL") + 
  scale_colour_manual(values = cc[c(1,2)]) #+ 
 # geom_vline(xintercept = dates_to_cut)

#### SAVE FIGURE 1 #####

tiff(here::here("Figures","FIGX_description"), height = 2000, width = 3200, res = 300)

grid.arrange(TOTAL_REPORTS, TIMESERIES, SEASONALITY,
             layout_matrix = rbind(c(2,2),
                                   c(1,3)) )
dev.off()



###### Expected casses #####


dat_cmc$INF_0_incidence <- (dat_cmc$Hos_INF_0 + dat_cmc$both_0)/ population_numbers[1]

dat_cmc$RSV_0_incidence <- (dat_cmc$Hos_RSV_0+ dat_cmc$both_0) / population_numbers[1]

dat_cmc$Dual_0_incidence <- dat_cmc$RSV_0_incidence*dat_cmc$INF_0_incidence 

dual_e_0 <-sum(dat_cmc$Dual_0_incidence)*population_numbers[1]


dat_cmc$INF_1_incidence <- (dat_cmc$Hos_INF_1+ dat_cmc$both_1) / population_numbers[2]

dat_cmc$RSV_1_incidence <- (dat_cmc$Hos_RSV_1+ dat_cmc$both_1) / population_numbers[2]

dat_cmc$Dual_1_incidence <- dat_cmc$RSV_1_incidence*dat_cmc$INF_1_incidence 

dual_e_1 <- sum(dat_cmc$Dual_1_incidence)*population_numbers[2]


##### z statistic #####

inf_0 <- sum(dat_cmc$Hos_INF_0)
inf_1 <- sum(dat_cmc$Hos_INF_1)
rsv_0<- sum(dat_cmc$Hos_RSV_0)
rsv_1<- sum(dat_cmc$Hos_RSV_1)
dual_0<- sum(dat_cmc$both_0)
dual_1<- sum(dat_cmc$both_1)
all_rsv_0 <- rsv_0 + dual_0
all_rsv_1 <- rsv_1 + dual_1
all_inf_0 <- inf_0 + dual_0
all_inf_1 <- inf_1 + dual_1
tot_0 <- sum(inf_0, rsv_0, dual_0)
tot_1 <- sum(inf_1, rsv_1, dual_1)
dual_e_0
dual_e_1

# youngest 
((dual_0/tot_0) - (dual_e_0/tot_0))/
  sqrt(((dual_0/tot_0) *(1 - dual_e_0/tot_0))/ tot_0)

# older
((dual_1/tot_1) - (dual_e_1/tot_1))/
  sqrt(((dual_1/tot_1) *(1 - dual_e_1/tot_1))/ tot_1)

########
dual_0/dual_e_0
dual_1/dual_e_1



####### Calculating the Attack Rate RSV #######


prev_RSV_0 <- (dual_0/tot_0) /(all_inf_0/tot_0 )

RSV_reporting <-  all_rsv_0/tot_0  * ((1/gammaR)/7) / prev_RSV_0 


dat_cmc[, weekly_total_0 := Hos_INF_0  + Hos_RSV_0 + both_0 ]
dat_cmc[, prop_dual := both_0 / weekly_total_0]
dat_cmc$prop_dual[is.nan(dat_cmc$prop_dual)] <- 0


calc_dual_ARI_0 <- function(RSV_reporting, input_table){
  
  #prevalence of RSV
  input_table[, prev_RSV_0 := ((((Hos_RSV_0 + both_0)/population_numbers[1]) * ((1/gammaR)/7) )/ RSV_reporting)]
  #incidence of dual in ARI
  input_table[, dual_ARI_0 := ((Hos_INF_0 + both_0)/population_numbers[1]) * prev_RSV_0 ]
  # remove NaNs
  input_table$dual_ARI_0[is.nan(input_table$dual_ARI_0)] <- 0
  #likelihood
  # input_table[, ll := dnorm(x = dual_ARI_0, mean = prop_dual, log =T)]
  print(RSV_reporting)
  input_table[, ll := dnbinom(x = both_0,
                              mu = dual_ARI_0*population_numbers[1], 
                              size = 1,  log=T)]
  
  return(-sum(input_table$ll))
}

annual_rates <- dat_cmc[, sum(Hos_RSV_0+both_0), by = year]

reporting_rate_0_all <- optim(par = 0.1, fn = calc_dual_ARI_0, method = "Brent", 
                              lower = c(0), upper = c(10),hessian = T, input_table = dat_cmc)
reporting_rate_0 <- reporting_rate_0_all$par

Attack_rate_0 <- ((annual_rates$V1)/population_numbers[1])/reporting_rate_0

standard_error_0 <- sqrt(diag(solve(reporting_rate_0_all$hessian)))
CI_upper_0 <- reporting_rate_0 + 1.96*standard_error_0
CI_lower_0 <-reporting_rate_0 - 1.96*standard_error_0
Attack_rate_upper <- ((annual_rates$V1)/population_numbers[1])/CI_upper_0
Attack_rate_lower <-((annual_rates$V1)/population_numbers[1])/CI_lower_0

print(paste0( "the average annual attack rate is ", mean(Attack_rate_0), " (", mean(Attack_rate_lower), " - ", mean(Attack_rate_upper), ")"))


# second age group

dat_cmc[, weekly_total_1 := Hos_INF_1  + Hos_RSV_1 + both_1 ]
dat_cmc[, prop_dual_1 := both_1 / weekly_total_1]
dat_cmc$prop_dual_1[is.nan(dat_cmc$prop_dual_1)] <- 0

calc_dual_ARI_1 <- function(RSV_reporting, input_table){
   RSV_reporting <- RSV_reporting/100
  #prevalence of RSV
  input_table[, prev_RSV_1 := ((((Hos_RSV_1 + both_1)/population_numbers[2]) * ((1/gammaR)/7) )/ RSV_reporting)]
  #incidence of dual in ARI
  input_table[, dual_ARI_1 := ((Hos_INF_1 + both_1)/population_numbers[2]) * prev_RSV_1 ]
  # remove NaNs
  input_table$dual_ARI_1[is.nan(input_table$dual_ARI_1)] <- 0
  #likelihood
  # input_table[, ll := dnorm(x = dual_ARI_1, mean = prop_dual, log =T)]
  print(RSV_reporting)
  input_table[, ll := dnbinom(x = both_1,
                              mu = dual_ARI_1*population_numbers[2], 
                              size = 1,  log=T)]

  return(-sum(input_table$ll))
}

annual_rates <- dat_cmc[, sum(Hos_RSV_1 + both_1), by = year]

reporting_rate_1_all <- optim(par = 100, fn = calc_dual_ARI_1, method = "Brent", 
                              lower = c(0), upper = c(10000),hessian = T, input_table = dat_cmc)
reporting_rate_1 <- reporting_rate_1_all$par

Attack_rate_1 <- ((annual_rates$V1)/population_numbers[2])/reporting_rate_1 *100

standard_error_1 <- sqrt(diag(solve(reporting_rate_1_all$hessian)))
CI_upper_1 <- reporting_rate_1 + 1.96*standard_error_1
CI_lower_1 <-reporting_rate_1 - 1.96*standard_error_1
Attack_rate_upper <- ((annual_rates$V1)/population_numbers[2])/CI_upper_1 *100
Attack_rate_lower <-((annual_rates$V1)/population_numbers[2])/CI_lower_1 *100

print(paste0( "the average annual attack rate is ", mean(Attack_rate_1), " (", mean(Attack_rate_lower), " - ", mean(Attack_rate_upper), ")"))



####### Calculating the Attack Rate Influenza #######


prev_INF_0 <- (dual_0/tot_0) /(all_inf_0/tot_0 )

INF_reporting <-  all_inf_0/tot_0  * ((1/gammaI)/7) / prev_INF_0 


dat_cmc[, weekly_total_0 := Hos_INF_0  + Hos_INF_0 + both_0 ]
dat_cmc[, prop_dual := both_0 / weekly_total_0]
dat_cmc$prop_dual[is.nan(dat_cmc$prop_dual)] <- 0


calc_dual_ARI_0 <- function(INF_reporting, input_table){
  INF_reporting <- INF_reporting/100
  #prevalence of INF
  input_table[, prev_INF_0 := ((((Hos_INF_0 + both_0)/population_numbers[1]) * ((1/gammaI)/7) )/ INF_reporting)]
  #incidence of dual in ARI
  input_table[, dual_ARI_0 := ((Hos_INF_0 + both_0)/population_numbers[1]) * prev_INF_0 ]
  # remove NaNs
  input_table$dual_ARI_0[is.nan(input_table$dual_ARI_0)] <- 0
  #likelihood
  # input_table[, ll := dnorm(x = dual_ARI_0, mean = prop_dual, log =T)]
  print(INF_reporting)
  input_table[, ll := dnbinom(x = both_0,
                              mu = dual_ARI_0*population_numbers[1], 
                              size = 1,  log=T)]
  
  return(-sum(input_table$ll))
}

annual_rates <- dat_cmc[, sum(Hos_INF_0+both_0), by = year]

reporting_rate_0_all <- optim(par = 100, fn = calc_dual_ARI_0, method = "Brent", 
                              lower = c(0), upper = c(10000),hessian = T, input_table = dat_cmc)
reporting_rate_0 <- reporting_rate_0_all$par

Attack_rate_0 <- ((annual_rates$V1)/population_numbers[1])/reporting_rate_0 *100

standard_error_0 <- sqrt(diag(solve(reporting_rate_0_all$hessian)))
CI_upper_0 <- reporting_rate_0 + 1.96*standard_error_0
CI_lower_0 <-reporting_rate_0 - 1.96*standard_error_0
Attack_rate_upper <- ((annual_rates$V1)/population_numbers[1])/CI_upper_0 *100
Attack_rate_lower <-((annual_rates$V1)/population_numbers[1])/CI_lower_0 *100

print(paste0( "the average annual attack rate is ", mean(Attack_rate_0), " (", mean(Attack_rate_lower), " - ", mean(Attack_rate_upper), ")"))


# second age group

dat_cmc[, weekly_total_1 := Hos_INF_1  + Hos_INF_1 + both_1 ]
dat_cmc[, prop_dual_1 := both_1 / weekly_total_1]
dat_cmc$prop_dual_1[is.nan(dat_cmc$prop_dual_1)] <- 0

calc_dual_ARI_1 <- function(INF_reporting, input_table){
  INF_reporting <- INF_reporting/100
  #prevalence of INF
  input_table[, prev_INF_1 := ((((Hos_INF_1 + both_1)/population_numbers[2]) * ((1/gammaI)/7) )/ INF_reporting)]
  #incidence of dual in ARI
  input_table[, dual_ARI_1 := ((Hos_INF_1 + both_1)/population_numbers[2]) * prev_INF_1 ]
  # remove NaNs
  input_table$dual_ARI_1[is.nan(input_table$dual_ARI_1)] <- 0
  #likelihood
  # input_table[, ll := dnorm(x = dual_ARI_1, mean = prop_dual, log =T)]
  print(INF_reporting)
  input_table[, ll := dnbinom(x = both_1,
                              mu = dual_ARI_1*population_numbers[2], 
                              size = 1,  log=T)]
  
  return(-sum(input_table$ll))
}

annual_rates <- dat_cmc[, sum(Hos_INF_1 + both_1), by = year]

reporting_rate_1_all <- optim(par = 100, fn = calc_dual_ARI_1, method = "Brent", 
                              lower = c(0), upper = c(10000),hessian = T, input_table = dat_cmc)
reporting_rate_1 <- reporting_rate_1_all$par

Attack_rate_1 <- ((annual_rates$V1)/population_numbers[2])/reporting_rate_1 *100

standard_error_1 <- sqrt(diag(solve(reporting_rate_1_all$hessian)))
CI_upper_1 <- reporting_rate_1 + 1.96*standard_error_1
CI_lower_1 <-reporting_rate_1 - 1.96*standard_error_1
Attack_rate_upper <- ((annual_rates$V1)/population_numbers[2])/CI_upper_1 *100
Attack_rate_lower <-((annual_rates$V1)/population_numbers[2])/CI_lower_1 *100

print(paste0( "the average annual attack rate is ", mean(Attack_rate_1), " (", mean(Attack_rate_lower), " - ", mean(Attack_rate_upper), ")"))
