#### define season cuts #####

start_dates <-as.Date(c("2007-11-01", "2008-11-01","2009-11-01","2010-11-01","2011-11-01",
                        "2012-11-01","2013-11-01","2014-11-01","2015-11-01","2016-11-01"))
end_dates <- as.Date(c( "2008-05-01","2009-05-01","2010-05-01","2011-05-01",
                        "2012-05-01","2013-05-01","2014-05-01","2015-05-01","2016-05-01", 
                        "2017-05-01"))

dat_cmc[, all_cases := Hos_INF_0+ Hos_INF_1+ Hos_RSV_0+ Hos_RSV_1+ both_0+ both_1]
dat_cmc[, smoothed_cases := rollmean(all_cases, k=4, fill = NA)]

dat_cmc$start_week <- as.Date(dat_cmc$start_week)
dates_out <- c()

for(i in 1:length(start_dates)){
  table_temp <- dat_cmc[((start_week > start_dates[i]) & (start_week < end_dates[i]))]
  date_temp <- (table_temp[which.min(table_temp$smoothed_cases),start_week])
  dates_out <- c(dates_out, as.character(date_temp))
}