# NhaTrang_flu_rsv_interaction
Repository accompanying the paper "Evidence for Influenza and RSV interaction from 10 years of enhanced surveillance in Nha Trang, Vietnam, a modelling study."

All the code for generating the research and figures for the paper is included in this repository. The summarised data is available in the file Nha_Trang_flu_rsv.csv

See LICENSE file for licensing details.

Corresponding author: Naomi R Waterlow, naomi.waterlow1@lshtm.ac.uk

Package versions are:
lubridate_1.7.9.2  GGally_1.5.0       Momocs_1.3.0       forcats_0.5.0     
zoo_1.8-8          gridExtra_2.3      RColorBrewer_1.1-2 Rcpp_1.0.6        
tmvtnorm_1.4-10    gmm_1.6-5          sandwich_2.5-1     Matrix_1.2-18     
mvtnorm_1.1-1      data.table_1.13.4  ggplot2_3.3.2      deSolve_1.28      
MASS_7.3-51.6      reshape2_1.4.4     doParallel_1.0.15  iterators_1.0.12  
foreach_1.5.0      here_0.1           rvest_0.3.5        xml2_1.3.2  

All analysis except for the fitting was done in R version 4.0.0 on macOS Catalina 10.15.7. The fitting was done on R 3.4 on AWS ec2 machines, running a linux AMI. On AWS the package mvtnorm was version 1.0.8.
