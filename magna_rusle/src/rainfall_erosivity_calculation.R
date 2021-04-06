# Script to calculate R-factor, monthly R-factor and past yearly 
# rainfall erosivitys for UK precipiation stations (from hourly records)

# load required libraries
library(readxl)
library(tidyverse)
#install.packages("devtools")
#devtools::install_github("kvantas/hyetor")
library(hyetor)
library(tibble)
library(dplyr)
library(lubridate)

# read in the list of uk rain gauge stations 
station_names <- read_csv("~/magna_rusle/data/input/station_names.csv")
# create empty dataframe that will contain the processed data
station_rainfall_erosivity <- data.frame()
# set the directory where the data is located
dir <- file.path("~", "magna_rusle","data","input","prec_stations")

# create function calculate_erosivities that calculates required erosivities
calculate_erosivities = function(name) {
  
  # define file path where station rainfall data are located
  file_path = file.path(dir,name)
  setwd(file_path)
  
  # create a list of files for each station
  file_list <- list.files(path=file_path)
  
  # create empty dataframe
  station <- data.frame()
  
  # process and combine files for each year
  for (i in 1:length(file_list)){
    data = read.csv(file_list[i], skip = 61) # skip metadata
    data = data[-nrow(data),] # remove last row
    data = data[data$ob_hour_count==1,] # keep only hourly summaries
    data = data[,c(1,7,9,10)] # remove unnescessary rows
    data$file = i # record the file number data came from
    station = rbind(station,data) # append rows
  }
  
  # ensure all observation times are unique
  station = station %>% 
    group_by(ob_end_time) %>%
    summarise(src_id=unique(src_id),prcp_amt=mean(prcp_amt,na.rm=TRUE),prc_dur=mean(prcp_dur,na.rm=TRUE),
              file=unique(file))
  
  # obtain first and last date/time
  start = station[1,1]$ob_end_time
  end = station[nrow(station),1]$ob_end_time
  
  # create variables for each timestamp
  station$ob_end_time = as.POSIXct(station$ob_end_time, format = "%Y-%m-%d %H:%M")
  station$year = year(station$ob_end_time) 
  station$month = month(station$ob_end_time)
  station$day = day(station$ob_end_time) 
  station$hour = hour(station$ob_end_time)
  
  #  create dataframe with all hours
  NoOfHours <- as.numeric(ymd_hms(end) - ymd_hms(start))*24 
  ob_end_time = ymd_hms(start) + hours(0:NoOfHours)
  ob_end_time = data.frame(ob_end_time)
  
  # join station and all hours
  imputed = left_join(ob_end_time,station)
  imputed$year = year(imputed$ob_end_time) 
  imputed$month = month(imputed$ob_end_time)
  imputed$day = day(imputed$ob_end_time) 
  imputed$hour = hour(imputed$ob_end_time)
  
  # calculate percentage of missing values for each year
  missing = imputed %>% 
    group_by(year) %>%
    summarise(n=n(),missing=sum(is.na(prcp_amt)),perc_missing=100*missing/n)
  # get a list of years with less that 10% missing data - can change this value
  # for higher quality results (less missingness)
  bad_quality_years = missing$year[missing$perc_missing>10]
  # save info on percentage of missing data for each year for each station
  x =paste("missing",name,sep="_")
  missing_path = paste("~/magna_rusle/reports/output_data/station_missing_data/",x,".csv",sep="")
  write.csv(missing,missing_path, row.names = FALSE)
  
  # remove bad quality years
  station = station[!(station$year %in% bad_quality_years),]
  station = station[!is.na(station$ob_end_time),]
  
  # create and process prec time series for erosivity calculation
  prec60min = station[,c(1,3)]
  prec60min$ob_end_time = as.POSIXct(prec60min$ob_end_time, format = "%Y-%m-%d %H:%M")
  colnames(prec60min) = c("date","prec")
  
  # breakdown or individual storm and erosivities
  ei_values <- prec60min %>%
    hyet_fill(time_step = 60, ts_unit = "mins") %>%
    hyet_erosivity(time_step = 60)
  
  # add years and months variables
  ei_values <- ei_values %>%
    mutate( year = year(begin)) %>%
    mutate( month = month(begin))
  
  # calculate long term R coeff
  R_coeff <- ei_values %>%
    group_by(year) %>%
    summarise(R = sum(erosivity)*1.5597) %>%
    ungroup() %>%
    summarise(R = mean(R)) %>%
    unlist()
  
  # summarise yearly erosivities
  yearly = ei_values %>% 
    group_by(year) %>%
    summarise(sum(erosivity)*1.5597)
  
  # summarise monthly erosivities
  monthly = ei_values %>%
    group_by(month) %>%
    summarise((sum(erosivity)*1.5597)/length(unique(ei_values$year)))
  
  # make column names the same and join together
  colnames(yearly) = c("time","erosivity")
  colnames(monthly) = c("time","erosivity")
  monthly$time = c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
  all = rbind(yearly,monthly)
  
  # transpose 
  all <- as.data.frame(t(as.matrix(all)))
  colnames(all) = all[1,]
  all = all[-1,]
  all$R = R_coeff
  # get station information
  station_id = unique(station$src_id)
  station =station_names[station_names$src_id==station_id,]
  all = cbind(station,all)
  rownames(all)= c()
  all$years_of_hourly = length(unique(ei_values$year))
  return(all)
}

all_stations = list.files(path="~/magna_rusle/data/input/prec_stations")
for (i in 1:length(all_stations)){
  station_rainfall_erosivity = bind_rows(station_rainfall_erosivity,
                                         calculate_erosivities(all_stations[i]))
}

# remove years with less than 15 years of hourly data

station_re = station_rainfall_erosivity[station_rainfall_erosivity$years_of_hourly>=15,]

# format the data a bit

col_order <- c("src_id","station_name","station_file_name",
               "historic_county","station_latitude",
               "station_longitude","station_elevation",
               "first_year","last_year","1995","1996",
               "1997","1998","1999","2000","2001",
               "2002","2003","2004","2005","2006",
               "2007","2008","2009","2010","2011",
               "2012","2013","2014","2015","2016",
               "2017","2018","2019","jan","feb","mar",
               "apr","may","jun","jul","aug","sep","oct","nov",
               "dec","R","years_of_hourly")

station_re2 = station_re[, col_order]
colnames(station_re2) = c("src_id","station_name","station_file_name",
                                         "historic_county","station_latitude",
                                         "station_longitude","station_elevation",
                                         "first_year","last_year","R_1995","R_1996",
                                         "R_1997","R_1998","R_1999","R_2000","R_2001",
                                         "R_2002","R_2003","R_2004","R_2005","R_2006",
                                         "R_2007","R_2008","R_2009","R_2010","R_2011",
                                         "R_2012","R_2013","R_2014","R_2015","R_2016",
                                         "R_2017","R_2018","R_2019","jan","feb","mar",
                                         "apr","may","jun","jul","aug","sep","oct","nov",
                                         "dec","R","years_of_hourly")

# create dataframe with monthly r factors and save
monthly_r_factors = station_re2[,c("src_id","station_name","historic_county",
                                  "station_latitude","station_longitude",
                                  "station_elevation","years_of_hourly","R",
                                  "jan","feb","mar","apr","may","jun","jul",
                                  "aug","sep","oct","nov","dec")]
write.csv(monthly_r_factors, file ="~/magna_rusle/reports/output_data/monthly_station_r_factors.csv" , 
          row.names = FALSE)

# create dataframe with past erosivities and save
past_erosivity = station_re2[,c("src_id","station_name","historic_county",
                                "station_latitude","station_longitude",
                                "station_elevation","years_of_hourly","R",
                                "R_1995","R_1996",
                                "R_1997","R_1998","R_1999","R_2000","R_2001",
                                "R_2002","R_2003","R_2004","R_2005","R_2006",
                                "R_2007","R_2008","R_2009","R_2010","R_2011",
                                "R_2012","R_2013","R_2014","R_2015","R_2016",
                                "R_2017","R_2018","R_2019")]
write.csv(past_erosivity, file ="~/magna_rusle/reports/output_data/station_past_erosivity_1995-2019.csv" , 
          row.names = FALSE)
