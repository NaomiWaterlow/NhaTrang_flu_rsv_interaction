#R0 Calculations
# Influenza only, no RSV in model. (I.e. assume model has only 4 compartments)
#Transmission Matrix

calculate_R0_rat_I <- function(theta_season, season_pop){
  
  gammaI <- 0.2632
  bI_calc_R0 <- 1
  pop_sizes <- season_pop
  # currenntly using year 2014 data to calculate
  contactmatrix <- population_contacts
  ITransmission<- matrix(nrow=length(pop_sizes), ncol=length(pop_sizes))
  for (i in 1:length(pop_sizes)){
    for (j in 1:length(pop_sizes)) {
      ITransmission[i,j] <- bI_calc_R0*contactmatrix[i,j]*pop_sizes[j]
    }
  }
  #Transition matrix
  ITransition <-matrix(0,nrow = length(pop_sizes), ncol = length(pop_sizes) )
  diag(ITransition)<- gammaI
  #Inverse 
  ITransition_inverse<- ginv(ITransition)
  INGM <- ITransmission%*%ITransition_inverse
  IEigen<- unlist((eigen(INGM)[1]))
  R0ratI <-  max(IEigen)
  return(R0ratI)
}

# Now for RSV only
calculate_R0_rat_R <- function(theta_season, season_pop){
  
  gammaR <- 0.1111
  bR_calc_R0 <- 1
  pop_sizes <- season_pop
  # currently using 2014 data to calculate
  contactmatrix <- population_contacts
  rsv_sus <- rsv_sus
  
  RTransmission<- matrix(nrow=length(pop_sizes), ncol=length(pop_sizes))
  
  for (i in 1:length(pop_sizes)){
    for (j in 1:length(pop_sizes)) {
      RTransmission[i,j] <- rsv_sus[j]*bR_calc_R0*contactmatrix[i,j]*pop_sizes[j]
    }
  }
  #Transition matrix
  RTransition <-matrix(0,nrow = length(pop_sizes), ncol = length(pop_sizes) )
  diag(RTransition) <- gammaR
  #Inverse 
  RTransition_inverse<- ginv(RTransition)
  RNGM <- RTransmission%*%RTransition_inverse
  REigen<- unlist((eigen(RNGM)[1]))
  R0ratR <- max(REigen)
  return(R0ratR)
}