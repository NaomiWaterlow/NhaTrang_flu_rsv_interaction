# This contains the non-chainging parameter inputs for Model_for_data. 
# Contains combined inf and batch updating

#Define number of agegroups
age_groups <- c(0,2,5,16,65)

#This will need altering. 
#source("PHE_Model_for_data_Mixing_UK_Ages.R") #UK Polymod contact matrix
#contactmatrix <- population_contacts


#fixed RSV susceptibility
rsv_sus <- c(1, 0.75, rep(0.65, (length(age_groups) - 2)))
# fixed recovery rates - also written inside calculating R0 rat functions
gammaR <- 0.1111
gammaI <- 0.2632

#calculate priors

season_params <- data.frame()
season1 <- c("y_07_08",
             "betaI",
             "Idetect",
             "propR_1",
             "inf_sus_1", 
             "Rdetect2", 
             "Rdetect5", 
             "sig", 
             "rho", 
             "betaR", 
             "propI_1",
             "detect_multiplier", 
             "RSV_default")
season2 <- c("y_08_09",
             "betaI",
             "Idetect",
             "propR_2",
             "inf_sus_2", 
             "Rdetect2", 
             "Rdetect5", 
             "sig", 
             "rho", 
             "betaR", 
             "propI_2",
             "detect_multiplier", 
             "RSV_default")
season3 <- c("y_09_10",
             "betaI",
             "Idetect",
             "propR_3",
             "inf_sus_3", 
             "Rdetect2", 
             "Rdetect5", 
             "sig", 
             "rho", 
             "betaR", 
             "propI_3",
             "detect_multiplier", 
             "RSV_default")
season4 <- c("y_10_11",
             "betaI",
             "Idetect",
             "propR_4",
             "inf_sus_4", 
             "Rdetect2", 
             "Rdetect5", 
             "sig", 
             "rho", 
             "betaR", 
             "propI_4",
             "detect_multiplier", 
             "RSV_default")
season5 <- c("y_11_12",
             "betaI",
             "Idetect",
             "propR_5",
             "inf_sus_5", 
             "Rdetect2", 
             "Rdetect5", 
             "sig", 
             "rho", 
             "betaR", 
             "propI_5",
             "detect_multiplier", 
             "RSV_default")
season6 <- c("y_12_13",
             "betaI",
             "Idetect",
             "propR_6",
             "inf_sus_6", 
             "Rdetect2", 
             "Rdetect5", 
             "sig", 
             "rho", 
             "betaR", 
             "propI_6",
             "detect_multiplier", 
             "RSV_detect_multiplier")
season7 <- c("y_13_14",
             "betaI",
             "Idetect",
             "propR_7",
             "inf_sus_7", 
             "Rdetect2", 
             "Rdetect5", 
             "sig", 
             "rho", 
             "betaR", 
             "propI_7",
             "detect_multiplier", 
             "RSV_detect_multiplier")
season8 <- c("y_14_15",
             "betaI",
             "Idetect",
             "propR_8",
             "inf_sus_8", 
             "Rdetect2", 
             "Rdetect5", 
             "sig", 
             "rho", 
             "betaR", 
             "propI_8",
             "detect_multiplier", 
             "RSV_detect_multiplier")
season9 <- c("y_15_16",
             "betaI",
             "Idetect",
             "propR_9",
             "inf_sus_9", 
             "Rdetect2", 
             "Rdetect5", 
             "sig", 
             "rho", 
             "betaR", 
             "propI_9",
             "detect_multiplier", 
             "RSV_detect_multiplier")
season10 <- c("y_16_17",
              "betaI",
              "Idetect",
              "propR_10",
              "inf_sus_10", 
              "Rdetect2", 
              "Rdetect5", 
              "sig", 
              "rho", 
              "betaR", 
              "propI_10",
              "detect_multiplier", 
              "RSV_detect_multiplier")
season11 <- c("y_17_18",
              "betaI",
              "Idetect",
              "propR_11",
              "inf_sus_11", 
              "Rdetect2", 
              "Rdetect5", 
              "sig", 
              "rho", 
              "betaR", 
              "propI_11",
              "detect_multiplier", 
              "RSV_detect_multiplier")


season_params <- rbind(season1, season2, season3, season4, season5, season6,season7, season8,
                       season9,
                       season10,season11)
colnames(season_params) <- c("year","betaI", "Idetect","propR", "inf_sus", 
                             "Rdetect2", 
                             "Rdetect5", 
                             "sig", 
                             "rho", 
                             "betaR", 
                             "propI",
                             "multiplier", 
                             "RSV_detect_multiplier")




#year <- multi_data[[season_number]][1,1]
