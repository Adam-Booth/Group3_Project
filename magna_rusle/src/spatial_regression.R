library(mgcv)
library(tidyverse)
library(caret)
library(raster)
library(ModelMetrics)
library(spaMM)
library(sf)
library(GGally)

monthly_r_factors <- read_csv("~/magna_rusle/reports/output_data/monthly_station_r_factors.csv")

# get extent of the fort
fort_boundary <- st_read(
  "~/magna_rusle/data/input/fort_shapefile/magna_fort_boundary.shp")
plot(fort_boundary)

# get extent of UK 
gbr <- raster::getData("GADM", country="GBR", level=1) 

# crop monthly mean precipiation data
prec_monthly_crop = crop(prec_monthly,gbr)

# crop monthly mean tmin
tmin_monthly_crop = crop(tmin_monthly,gbr)

# crop monthly mean tmax
tmax_monthly_crop = crop(tmax_monthly,gbr)

# crop monthly mean bioclim vars
biovars_crop = crop(biovars_current,gbr)

# get elevation data 
elev = stack("~/magna_rusle/data/input/wc2.1_30s_elev.tif")
elev = crop(elev,gbr)
names(elev) = "station_elevation"

# make a copy of the R factor database
r_data = monthly_r_factors

# create new variables by extracting monthly mean prec, tmin, tmax and biovars
# at each site

# create list of months in order
months = c("jan", "feb", "mar", "apr" ,"may", "jun", "jul",
           "aug", "sep", "oct" ,"nov", "dec")

# create columns for mean monthly prec
for (i in 1:12) {
  x=paste("prec",months[i],sep="_")
  r_data$new = extract(prec_monthly_crop[[i]], r_data[,c("station_longitude", "station_latitude")])
  colnames(r_data)[which(names(r_data) == "new")] <- x
}

# create columns for mean monthly tmin
for (i in 1:12) {
  y=paste("tmin",months[i],sep="_")
  r_data$new = extract(tmin_monthly_crop[[i]], r_data[,c("station_longitude", "station_latitude")])
  colnames(r_data)[which(names(r_data) == "new")] <- y
}

# create columns for mean monthly tmax
for (i in 1:12) {
  z=paste("tmax",months[i],sep="_")
  r_data$new = extract(tmax_monthly_crop[[i]], r_data[,c("station_longitude", "station_latitude")])
  colnames(r_data)[which(names(r_data) == "new")] <- z
}

# create columns for biovars
for (i in 1:19) {
  b=paste("bio",i,sep="")
  r_data$new = extract(biovars_crop[[i]], r_data[,c("station_longitude", "station_latitude")])
  colnames(r_data)[which(names(r_data) == "new")] <- b
}

# Create an empty raster with the same extent and resolution as climate layers
latitude_raster = longitude_raster = raster(nrows = nrow(prec_monthly_crop),
                                            ncols = ncol(prec_monthly_crop),
                                            ext = extent(prec_monthly_crop))

# Change the values to be latitude and longitude respectively
longitude_raster[] <- coordinates(longitude_raster)[,1]
latitude_raster[] <- coordinates(latitude_raster)[,2]
names(longitude_raster) = c("station_longitude")
names(latitude_raster) = c("station_latitude")

# Create a prediction stack of all variables

# first ensure all names are the same
names(prec_monthly_crop) = c("prec_jan","prec_feb","prec_mar","prec_apr",
                             "prec_may","prec_jun","prec_jul","prec_aug",
                             "prec_sep","prec_oct","prec_nov","prec_dec")
names(tmin_monthly_crop) = c("tmin_jan","tmin_feb","tmin_mar","tmin_apr",
                             "tmin_may","tmin_jun","tmin_jul","tmin_aug",
                             "tmin_sep","tmin_oct","tmin_nov","tmin_dec")
names(tmax_monthly_crop) = c("tmax_jan","tmax_feb","tmax_mar","tmax_apr",
                             "tmax_may","tmax_jun","tmax_jul","tmax_aug",
                             "tmax_sep","tmax_oct","tmax_nov","tmax_dec")

plot(elev)
elev2 <- resample(elev, biovars_crop, method='bilinear')
plot(elev2)
pred_stack <- stack(prec_monthly_crop,
                    tmin_monthly_crop,
                    tmax_monthly_crop,
                    biovars_crop,
                    longitude_raster,
                    latitude_raster,
                    elev2)
names(pred_stack)

# validate the model
# split dataset into train and test dataset (80% split)
# and perform k-fold cross validation

# create empty dataframes to store cv stats
rmse = data.frame()
R2 = data.frame()
mae = data.frame()

# remove one column that has NA for all climate vars
r_data2 = na.omit(r_data)
r_data2_monthly = r_data2[,1:20]
#write.csv(r_data2_monthly, file ="~/magna_climate/data/future_climate/predictions/R_monthly_database.csv" , 
#          row.names = FALSE)


# plot pairs and correlation between vars - this was used to help decide which covariates to use in models
#ggpairs(r_data2[c(4:6,8)])
#ggpairs(r_data2[c(8,57:60)])
#ggpairs(r_data2[c(8,60:65)])


# k fold cross validation
for (i in 1:10) {
  # set seed to ensure reproducible
  set.seed(666+9*i)
  training_samples <- r_data2$src_id %>%
    createDataPartition(p = 0.8, list = FALSE)
  train_data  <- r_data2[training_samples, ]
  valid_data <- r_data2[-training_samples, ]
  model <- fitme(R ~ prec_mar + prec_apr + prec_jan + prec_feb + prec_may +
                   prec_jun + prec_jul + prec_aug + prec_sep + prec_oct +
                   prec_nov + prec_dec + station_elevation + bio2 + bio3 + bio4+
                   bio15 + bio18 + bio19 +
                   Matern(1|station_latitude+station_longitude), 
                 family=gaussian(link="log"), data = train_data)
  predictions_validation <- predict(model, valid_data)
  ggplot() + geom_point(aes(as.vector(predictions_validation), valid_data$R))
  # Calculate rmse and r^2
  rmse = rbind(rmse,rmse(predictions_validation, valid_data$R))
  R2 = rbind(R2,R2(predictions_validation, valid_data$R))
  mae = rbind(mae,mae(predictions_validation, valid_data$R))
}

# create dataframe containing all cross validation stats for the models
model = "R"
cross_validation = cbind(R2,rmse,mae)
colnames(cross_validation)
cv_stats = colMeans(cross_validation)
cross_validation_stats = data.frame(model,cv_stats[1],cv_stats[2],cv_stats[3])
colnames(cross_validation_stats) = c("model","R2","rmse","mae")

# retrain model with all the data
model <- fitme(R ~ prec_mar + prec_apr + prec_jan + prec_feb + prec_may +
                 prec_jun + prec_jul + prec_aug + prec_sep + prec_oct +
                 prec_nov + prec_dec + station_elevation +  bio2 + bio3 + bio4+
                 bio15 + bio18 + bio19 +
                 Matern(1|station_latitude+station_longitude), 
               family=gaussian(link="log"), data = r_data2)

# interpolate current map of R factors
current_r_pred <- predict(pred_stack, model)
current_r_pred <- mask(current_r_pred, gbr)
plot(current_r_pred)
setwd("~/magna_rusle/reports/output_data/r_factor_predictions")
writeRaster(current_r_pred,filename="current_r_factor_pred",format="GTiff",overwrite=TRUE)
current_magna = crop(current_r_pred,fort_boundary)
mean_magna = cellStats(current_magna, stat='mean', na.rm=TRUE)
magna_years = c("2001-2020")
magna_ssps = c(NA)
magna_rs = c(mean_magna)

magna_R = data.frame(magna_years,magna_ssps,magna_rs)

# function which predicts future R factor for years and ssp specified
dir =  "~/magna_rusle/data/input/future_climate"

# create empty dataframe for storing means

predict_future_R = function(years,ssp) {
  
  yrs = paste("years",years,sep="_")
  path = file.path(dir,"CNRM-CM6-1",yrs,ssp)
  setwd(path)
  
  bi = paste("wc2.1_2.5m_bioc_CNRM-CM6-1_",ssp,"_",years,".tif",sep="")
  pr = paste("wc2.1_2.5m_prec_CNRM-CM6-1_",ssp,"_",years,".tif",sep="")
  tn = paste("wc2.1_2.5m_tmin_CNRM-CM6-1_",ssp,"_",years,".tif",sep="")
  tx = paste("wc2.1_2.5m_tmax_CNRM-CM6-1_",ssp,"_",years,".tif",sep="")
  
  bioc_stack2 = stack(bi)
  bioc_stack2 = crop(bioc_stack2,gbr)
  prec_stack2 = stack(pr)
  prec_stack2 = crop(prec_stack2,gbr)
  tmin_stack2 = stack(tn)
  tmin_stack2 = crop(tmin_stack2,gbr)
  tmax_stack2 = stack(tx)
  tmax_stack2 = crop(tmax_stack2,gbr)
  
  names(prec_stack2) = c("prec_jan","prec_feb","prec_mar","prec_apr",
                         "prec_may","prec_jun","prec_jul","prec_aug",
                         "prec_sep","prec_oct","prec_nov","prec_dec")
  names(tmin_stack2) = c("tmin_jan","tmin_feb","tmin_mar","tmin_apr",
                         "tmin_may","tmin_jun","tmin_jul","tmin_aug",
                         "tmin_sep","tmin_oct","tmin_nov","tmin_dec")
  names(tmax_stack2) = c("tmax_jan","tmax_feb","tmax_mar","tmax_apr",
                         "tmax_may","tmax_jun","tmax_jul","tmax_aug",
                         "tmax_sep","tmax_oct","tmax_nov","tmax_dec")
  names(bioc_stack2) = c("bio1","bio2","bio3","bio4","bio5","bio6","bio7",
                         "bio8","bio9","bio10","bio11","bio12","bio13","bio14",
                         "bio15","bio16","bio17","bio18","bio19")
  pred_stack2 <- stack(prec_stack2,
                       tmin_stack2,
                       tmax_stack2,
                       bioc_stack2,
                       longitude_raster,
                       latitude_raster,
                       elev2)
  names(pred_stack2)
  
  predicted_raster = predict(pred_stack2,model)
  predicted_raster = mask(predicted_raster, gbr)
  plot(predicted_raster)
  x = paste("future_r_factor_pred",years,ssp,sep="_")
  setwd("~/magna_rusle/reports/output_data/r_factor_predictions")
  writeRaster(predicted_raster,filename=x,format="GTiff",overwrite=TRUE)
  
  rel_change = 100*((predicted_raster - current_r_pred)/current_r_pred)
  y = paste(x,"relative_change")
  writeRaster(rel_change,filename=y,format="GTiff",overwrite=TRUE)
  
  at_magna = crop(predicted_raster,fort_boundary)
  avgr = cellStats(at_magna, stat='mean', na.rm=TRUE)
  new_row = c(years,ssp,avgr)
  return(new_row)
}

years_list=c("2021-2040","2041-2060","2061-2080","2081-2100")
ssp_list=c("ssp126","ssp245","ssp585")

# calculate predictions for every ssp and all years
for (i in 1:length(years_list)) {
  for (j in 1:length(ssp_list)) {
    magna_R[nrow(magna_R) +1,] = predict_future_R(years_list[i],ssp_list[j])
  }
}

# change column names and save
colnames(magna_R) = c("Years","SSP","R_factor_Magna")
write.csv(magna_R, file ="~/magna_rusle/reports/output_data/r_factor_predictions/R_factor_magna.csv" , 
          row.names = FALSE)

# Next calculate predicted monthly R factors for each SSP 

# create empty dataframes to store cv stats
rmse_jan = data.frame()
R2_jan = data.frame()
mae_jan = data.frame()

rmse_feb = data.frame()
R2_feb = data.frame()
mae_feb = data.frame()

rmse_mar = data.frame()
R2_mar = data.frame()
mae_mar = data.frame()

rmse_apr = data.frame()
R2_apr = data.frame()
mae_apr = data.frame()

rmse_may = data.frame()
R2_may = data.frame()
mae_may = data.frame()

rmse_jun = data.frame()
R2_jun = data.frame()
mae_jun = data.frame()

rmse_jul = data.frame()
R2_jul = data.frame()
mae_jul = data.frame()

rmse_aug = data.frame()
R2_aug = data.frame()
mae_aug = data.frame()

rmse_sep = data.frame()
R2_sep = data.frame()
mae_sep = data.frame()

rmse_oct = data.frame()
R2_oct = data.frame()
mae_oct = data.frame()

rmse_nov = data.frame()
R2_nov = data.frame()
mae_nov = data.frame()

rmse_dec = data.frame()
R2_dec = data.frame()
mae_dec = data.frame()

r_data2$jan = as.numeric(r_data2$jan)
r_data2$feb = as.numeric(r_data2$feb)
r_data2$mar = as.numeric(r_data2$mar)
r_data2$apr = as.numeric(r_data2$apr)
r_data2$may = as.numeric(r_data2$may)
r_data2$jun = as.numeric(r_data2$jun)
r_data2$jul = as.numeric(r_data2$jul)
r_data2$aug = as.numeric(r_data2$aug)
r_data2$sep = as.numeric(r_data2$sep)
r_data2$oct = as.numeric(r_data2$oct)
r_data2$nov = as.numeric(r_data2$nov)
r_data2$dec = as.numeric(r_data2$dec)

for (i in 1:10) {
  # set seed to ensure repeatable
  set.seed(666+9*i)
  training_samples <- r_data2$src_id %>%
    createDataPartition(p = 0.8, list = FALSE)
  train_data  <- r_data2[training_samples, ]
  valid_data <- r_data2[-training_samples, ]
  model_jan <- fitme(jan ~ prec_jan + station_elevation + tmin_jan + tmax_jan +
                       Matern(1|station_latitude+station_longitude),
                     family=gaussian(link="log"), data = train_data)
  model_feb <- fitme(feb ~ prec_feb + station_elevation + tmin_feb + tmax_feb +
                       Matern(1|station_latitude+station_longitude),
                     family=gaussian(link="log"), data = train_data)
  model_mar <- fitme(mar ~ prec_mar + station_elevation + tmin_mar + tmax_mar +
                       Matern(1|station_latitude+station_longitude),
                     family=gaussian(link="log"), data = train_data)
  model_apr <- fitme(apr ~ prec_apr + station_elevation + tmin_apr + tmax_apr +
                       Matern(1|station_latitude+station_longitude),
                     family=gaussian(link="log"), data = train_data)
  model_may <- fitme(may ~ prec_may + station_elevation + tmin_may + tmax_may +
                       Matern(1|station_latitude+station_longitude),
                     family=gaussian(link="log"), data = train_data)
  model_jun <- fitme(jun ~ prec_jun + station_elevation + tmin_jun + tmax_jun +
                       Matern(1|station_latitude+station_longitude),
                     family=gaussian(link="log"), data = train_data)
  model_jul <- fitme(jul ~ prec_jul + station_elevation + tmin_jul + tmax_jul +
                       Matern(1|station_latitude+station_longitude),
                     family=gaussian(link="log"), data = train_data)
  model_aug <- fitme(aug ~ prec_aug + station_elevation + tmin_aug + tmax_aug + 
                       Matern(1|station_latitude+station_longitude),
                     family=gaussian(link="log"), data = train_data)
  model_sep <- fitme(sep ~ prec_sep + station_elevation + tmin_sep + tmax_sep +
                       Matern(1|station_latitude+station_longitude),
                     family=gaussian(link="log"), data = train_data)
  model_oct <- fitme(oct ~ prec_oct + station_elevation + tmin_oct + tmax_oct +
                       Matern(1|station_latitude+station_longitude),
                     family=gaussian(link="log"), data = train_data)
  model_nov <- fitme(nov ~ prec_nov + station_elevation + tmin_nov + tmax_nov +
                       Matern(1|station_latitude+station_longitude),
                     family=gaussian(link="log"), data = train_data)
  model_dec <- fitme(dec ~ prec_dec + station_elevation + tmin_dec + tmax_dec +
                       Matern(1|station_latitude+station_longitude),
                     family=gaussian(link="log"), data = train_data)
  
  
  
  
  predictions_jan <- predict(model_jan, valid_data)
  predictions_feb <- predict(model_feb, valid_data)
  predictions_mar <- predict(model_mar, valid_data)
  predictions_apr <- predict(model_apr, valid_data)
  predictions_may <- predict(model_may, valid_data)
  predictions_jun <- predict(model_jun, valid_data)
  predictions_jul <- predict(model_jul, valid_data)
  predictions_aug <- predict(model_aug, valid_data)
  predictions_sep <- predict(model_sep, valid_data)
  predictions_oct <- predict(model_oct, valid_data)
  predictions_nov <- predict(model_nov, valid_data)
  predictions_dec <- predict(model_dec, valid_data)
  
  # Calculate rmse and r^2 for each model
  rmse_jan = rbind(rmse_jan,rmse(predictions_jan, valid_data$jan))
  R2_jan = rbind(R2_jan,R2(predictions_jan, valid_data$jan))
  mae_jan = rbind(mae_jan,mae(predictions_jan, valid_data$jan))
  
  rmse_feb = rbind(rmse_feb,rmse(predictions_feb, valid_data$feb))
  R2_feb = rbind(R2_feb,R2(predictions_feb, valid_data$feb))
  mae_feb = rbind(mae_feb,mae(predictions_feb, valid_data$feb))
  
  rmse_mar = rbind(rmse_mar,rmse(predictions_mar, valid_data$mar))
  R2_mar = rbind(R2_mar,R2(predictions_mar, valid_data$mar))
  mae_mar = rbind(mae_mar,mae(predictions_mar, valid_data$mar))
  
  rmse_apr = rbind(rmse_apr,rmse(predictions_apr, valid_data$apr))
  R2_apr = rbind(R2_apr,R2(predictions_apr, valid_data$apr))
  mae_apr = rbind(mae_apr,mae(predictions_apr, valid_data$apr))
  
  rmse_may = rbind(rmse_may,rmse(predictions_may, valid_data$may))
  R2_may = rbind(R2_may,R2(predictions_may, valid_data$may))
  mae_may = rbind(mae_may,mae(predictions_may, valid_data$may))
  
  rmse_jun = rbind(rmse_jun,rmse(predictions_jun, valid_data$jun))
  R2_jun = rbind(R2_jun,R2(predictions_jun, valid_data$jun))
  mae_jun = rbind(mae_jun,mae(predictions_jun, valid_data$jun))
  
  rmse_jul = rbind(rmse_jul,rmse(predictions_jul, valid_data$jul))
  R2_jul = rbind(R2_jul,R2(predictions_jul, valid_data$jul))
  mae_jul = rbind(mae_jul,mae(predictions_jul, valid_data$jul))
  
  rmse_aug = rbind(rmse_aug,rmse(predictions_aug, valid_data$aug))
  R2_aug = rbind(R2_aug,R2(predictions_aug, valid_data$aug))
  mae_aug = rbind(mae_aug,mae(predictions_aug, valid_data$aug))
  
  rmse_sep = rbind(rmse_sep,rmse(predictions_sep, valid_data$sep))
  R2_sep = rbind(R2_sep,R2(predictions_sep, valid_data$sep))
  mae_sep = rbind(mae_sep,mae(predictions_sep, valid_data$sep))
  
  rmse_oct = rbind(rmse_oct,rmse(predictions_oct, valid_data$oct))
  R2_oct = rbind(R2_oct,R2(predictions_oct, valid_data$oct))
  mae_oct = rbind(mae_oct,mae(predictions_oct, valid_data$oct))
  
  rmse_nov = rbind(rmse_nov,rmse(predictions_nov, valid_data$nov))
  R2_nov = rbind(R2_nov,R2(predictions_nov, valid_data$nov))
  mae_nov = rbind(mae_nov,mae(predictions_nov, valid_data$nov))
  
  rmse_dec = rbind(rmse_dec,rmse(predictions_dec, valid_data$dec))
  R2_dec = rbind(R2_dec,R2(predictions_dec, valid_data$dec))
  mae_dec = rbind(mae_dec,mae(predictions_dec, valid_data$dec))
}

# add to dataframe containing cross validation stats 
cv_jan = cbind(R2_jan,rmse_jan,mae_jan)
cv_jan = colMeans(cv_jan)
cv_row_jan = data.frame("jan",cv_jan[1],cv_jan[2],cv_jan[3])

cv_feb = cbind(R2_feb,rmse_feb,mae_feb)
cv_feb = colMeans(cv_feb)
cv_row_feb = data.frame("feb",cv_feb[1],cv_feb[2],cv_feb[3])

cv_mar = cbind(R2_mar,rmse_mar,mae_mar)
cv_mar = colMeans(cv_mar)
cv_row_mar = data.frame("mar",cv_mar[1],cv_mar[2],cv_mar[3])

cv_apr = cbind(R2_apr,rmse_apr,mae_apr)
cv_apr = colMeans(cv_apr)
cv_row_apr = data.frame("apr",cv_apr[1],cv_apr[2],cv_apr[3])

cv_may = cbind(R2_may,rmse_may,mae_may)
cv_may = colMeans(cv_may)
cv_row_may = data.frame("may",cv_may[1],cv_may[2],cv_may[3])

cv_jun = cbind(R2_jun,rmse_jun,mae_jun)
cv_jun = colMeans(cv_jun)
cv_row_jun = data.frame("jun",cv_jun[1],cv_jun[2],cv_jun[3])

cv_jul = cbind(R2_jul,rmse_jul,mae_jul)
cv_jul = colMeans(cv_jul)
cv_row_jul = data.frame("jul",cv_jul[1],cv_jul[2],cv_jul[3])

cv_aug = cbind(R2_aug,rmse_aug,mae_aug)
cv_aug = colMeans(cv_aug)
cv_row_aug = data.frame("aug",cv_aug[1],cv_aug[2],cv_aug[3])

cv_sep = cbind(R2_sep,rmse_sep,mae_sep)
cv_sep = colMeans(cv_sep)
cv_row_sep = data.frame("sep",cv_sep[1],cv_sep[2],cv_sep[3])

cv_oct = cbind(R2_oct,rmse_oct,mae_oct)
cv_oct = colMeans(cv_oct)
cv_row_oct = data.frame("oct",cv_oct[1],cv_oct[2],cv_oct[3])

cv_nov = cbind(R2_nov,rmse_nov,mae_nov)
cv_nov = colMeans(cv_nov)
cv_row_nov = data.frame("nov",cv_nov[1],cv_nov[2],cv_nov[3])

cv_dec = cbind(R2_dec,rmse_dec,mae_dec)
cv_dec = colMeans(cv_dec)
cv_row_dec = data.frame("dec",cv_dec[1],cv_dec[2],cv_dec[3])

colnames(cv_row_jan) = c("model","R2","rmse","mae")
colnames(cv_row_feb) = c("model","R2","rmse","mae")
colnames(cv_row_mar) = c("model","R2","rmse","mae")
colnames(cv_row_apr) = c("model","R2","rmse","mae")
colnames(cv_row_may) = c("model","R2","rmse","mae")
colnames(cv_row_jun) = c("model","R2","rmse","mae")
colnames(cv_row_jul) = c("model","R2","rmse","mae")
colnames(cv_row_aug) = c("model","R2","rmse","mae")
colnames(cv_row_sep) = c("model","R2","rmse","mae")
colnames(cv_row_oct) = c("model","R2","rmse","mae")
colnames(cv_row_nov) = c("model","R2","rmse","mae")
colnames(cv_row_dec) = c("model","R2","rmse","mae")

cv_months = rbind(cv_row_jan,cv_row_feb,cv_row_mar,cv_row_apr,cv_row_may,
                  cv_row_jun,cv_row_jul,cv_row_aug,cv_row_sep,cv_row_oct,
                  cv_row_nov,cv_row_dec)
cross_validation_stats = rbind(cross_validation_stats,cv_months)
write.csv(cross_validation_stats, 
          file ="~/magna_rusle/reports/output_data/r_factor_predictions/model_cross_validation_stats.csv" , 
          row.names = FALSE)

# retrain models using all data available

model_jan <- fitme(jan ~ prec_jan + station_elevation + tmin_jan + tmax_jan +
                     Matern(1|station_latitude+station_longitude),
                   family=gaussian(link="log"), data = r_data2)
model_feb <- fitme(feb ~ prec_feb + station_elevation + tmin_feb + tmax_feb +
                     Matern(1|station_latitude+station_longitude),
                   family=gaussian(link="log"), data = r_data2)
model_mar <- fitme(mar ~ prec_mar + station_elevation + tmin_mar + tmax_mar +
                     Matern(1|station_latitude+station_longitude),
                   family=gaussian(link="log"), data = r_data2)
model_apr <- fitme(apr ~ prec_apr + station_elevation + tmin_apr + tmax_apr +
                     Matern(1|station_latitude+station_longitude),
                   family=gaussian(link="log"), data = r_data2)
model_may <- fitme(may ~ prec_may + station_elevation + tmin_may + tmax_may +
                     Matern(1|station_latitude+station_longitude),
                   family=gaussian(link="log"), data = r_data2)
model_jun <- fitme(jun ~ prec_jun + station_elevation + tmin_jun + tmax_jun +
                     Matern(1|station_latitude+station_longitude),
                   family=gaussian(link="log"), data = r_data2)
model_jul <- fitme(jul ~ prec_jul + station_elevation + tmin_jul + tmax_jul +
                     Matern(1|station_latitude+station_longitude),
                   family=gaussian(link="log"), data = r_data2)
model_aug <- fitme(aug ~ prec_aug + station_elevation + tmin_aug + tmax_aug + 
                     Matern(1|station_latitude+station_longitude),
                   family=gaussian(link="log"), data = r_data2)
model_sep <- fitme(sep ~ prec_sep + station_elevation + tmin_sep + tmax_sep +
                     Matern(1|station_latitude+station_longitude),
                   family=gaussian(link="log"), data = r_data2)
model_oct <- fitme(oct ~ prec_oct + station_elevation + tmin_oct + tmax_oct +
                     Matern(1|station_latitude+station_longitude),
                   family=gaussian(link="log"), data = r_data2)
model_nov <- fitme(nov ~ prec_nov + station_elevation + tmin_nov + tmax_nov + 
                     Matern(1|station_latitude+station_longitude),
                   family=gaussian(link="log"), data = r_data2)
model_dec <- fitme(dec ~ prec_dec + station_elevation + tmin_dec + tmax_dec +
                     Matern(1|station_latitude+station_longitude),
                   family=gaussian(link="log"), data = r_data2)

# predict/interpolate current monthly erosivities using models

jan_r_pred <- predict(pred_stack, model_jan)
feb_r_pred <- predict(pred_stack, model_feb)
mar_r_pred <- predict(pred_stack, model_mar)
apr_r_pred <- predict(pred_stack, model_apr)
may_r_pred <- predict(pred_stack, model_may)
jun_r_pred <- predict(pred_stack, model_jun)
jul_r_pred <- predict(pred_stack, model_jul)
aug_r_pred <- predict(pred_stack, model_aug)
sep_r_pred <- predict(pred_stack, model_sep)
oct_r_pred <- predict(pred_stack, model_oct)
nov_r_pred <- predict(pred_stack, model_nov)
dec_r_pred <- predict(pred_stack, model_dec)


jan_r_pred <- mask(jan_r_pred, gbr)
feb_r_pred <- mask(feb_r_pred, gbr)
mar_r_pred <- mask(mar_r_pred, gbr)
apr_r_pred <- mask(apr_r_pred, gbr)
may_r_pred <- mask(may_r_pred, gbr)
jun_r_pred <- mask(jun_r_pred, gbr)
jul_r_pred <- mask(jul_r_pred, gbr)
aug_r_pred <- mask(aug_r_pred, gbr)
sep_r_pred <- mask(sep_r_pred, gbr)
oct_r_pred <- mask(oct_r_pred, gbr)
nov_r_pred <- mask(nov_r_pred, gbr)
dec_r_pred <- mask(dec_r_pred, gbr)

setwd("~/magna_rusle/reports/output_data/r_factor_predictions/monthly_maps/2001-2020")
writeRaster(jan_r_pred,filename="jan_r_factor_pred_2001-2020",format="GTiff",overwrite=TRUE)
writeRaster(feb_r_pred,filename="feb_r_factor_pred_2001-2020",format="GTiff",overwrite=TRUE)
writeRaster(mar_r_pred,filename="mar_r_factor_pred_2001-2020",format="GTiff",overwrite=TRUE)
writeRaster(apr_r_pred,filename="apr_r_factor_pred_2001-2020",format="GTiff",overwrite=TRUE)
writeRaster(may_r_pred,filename="may_r_factor_pred_2001-2020",format="GTiff",overwrite=TRUE)
writeRaster(jun_r_pred,filename="jun_r_factor_pred_2001-2020",format="GTiff",overwrite=TRUE)
writeRaster(jul_r_pred,filename="jul_r_factor_pred_2001-2020",format="GTiff",overwrite=TRUE)
writeRaster(aug_r_pred,filename="aug_r_factor_pred_2001-2020",format="GTiff",overwrite=TRUE)
writeRaster(sep_r_pred,filename="sep_r_factor_pred_2001-2020",format="GTiff",overwrite=TRUE)
writeRaster(oct_r_pred,filename="oct_r_factor_pred_2001-2020",format="GTiff",overwrite=TRUE)
writeRaster(nov_r_pred,filename="nov_r_factor_pred_2001-2020",format="GTiff",overwrite=TRUE)
writeRaster(dec_r_pred,filename="dec_r_factor_pred_2001-2020",format="GTiff",overwrite=TRUE)

jan_magna = crop(jan_r_pred,fort_boundary)
feb_magna = crop(feb_r_pred,fort_boundary)
mar_magna = crop(mar_r_pred,fort_boundary)
apr_magna = crop(apr_r_pred,fort_boundary)
may_magna = crop(may_r_pred,fort_boundary)
jun_magna = crop(jun_r_pred,fort_boundary)
jul_magna = crop(jul_r_pred,fort_boundary)
aug_magna = crop(aug_r_pred,fort_boundary)
sep_magna = crop(sep_r_pred,fort_boundary)
oct_magna = crop(oct_r_pred,fort_boundary)
nov_magna = crop(nov_r_pred,fort_boundary)
dec_magna = crop(dec_r_pred,fort_boundary)

mean_jan = cellStats(jan_magna, stat='mean', na.rm=TRUE)
mean_feb = cellStats(feb_magna, stat='mean', na.rm=TRUE)
mean_mar = cellStats(mar_magna, stat='mean', na.rm=TRUE)
mean_apr = cellStats(apr_magna, stat='mean', na.rm=TRUE)
mean_may = cellStats(may_magna, stat='mean', na.rm=TRUE)
mean_jun = cellStats(jun_magna, stat='mean', na.rm=TRUE)
mean_jul = cellStats(jul_magna, stat='mean', na.rm=TRUE)
mean_aug = cellStats(aug_magna, stat='mean', na.rm=TRUE)
mean_sep = cellStats(sep_magna, stat='mean', na.rm=TRUE)
mean_oct = cellStats(oct_magna, stat='mean', na.rm=TRUE)
mean_nov = cellStats(nov_magna, stat='mean', na.rm=TRUE)
mean_dec = cellStats(dec_magna, stat='mean', na.rm=TRUE)

magna_years = c("2001-2020")
magna_ssps = c(NA)
magna_monthly_r = data.frame(magna_years,magna_ssps,mean_jan,mean_feb,mean_mar,
                             mean_apr,mean_may,mean_jun,mean_jul,mean_aug,mean_sep,
                             mean_oct,mean_nov,mean_dec)

# next predict future monthly r factors for different ssps

predict_future_monthly_R = function(years,ssp) {
  
  yrs = paste("years",years,sep="_")
  path = file.path(dir,"CNRM-CM6-1",yrs,ssp)
  setwd(path)
  
  bi = paste("wc2.1_2.5m_bioc_CNRM-CM6-1_",ssp,"_",years,".tif",sep="")
  pr = paste("wc2.1_2.5m_prec_CNRM-CM6-1_",ssp,"_",years,".tif",sep="")
  tn = paste("wc2.1_2.5m_tmin_CNRM-CM6-1_",ssp,"_",years,".tif",sep="")
  tx = paste("wc2.1_2.5m_tmax_CNRM-CM6-1_",ssp,"_",years,".tif",sep="")
  
  bioc_stack2 = stack(bi)
  bioc_stack2 = crop(bioc_stack2,gbr)
  prec_stack2 = stack(pr)
  prec_stack2 = crop(prec_stack2,gbr)
  tmin_stack2 = stack(tn)
  tmin_stack2 = crop(tmin_stack2,gbr)
  tmax_stack2 = stack(tx)
  tmax_stack2 = crop(tmax_stack2,gbr)
  
  names(prec_stack2) = c("prec_jan","prec_feb","prec_mar","prec_apr",
                         "prec_may","prec_jun","prec_jul","prec_aug",
                         "prec_sep","prec_oct","prec_nov","prec_dec")
  names(tmin_stack2) = c("tmin_jan","tmin_feb","tmin_mar","tmin_apr",
                         "tmin_may","tmin_jun","tmin_jul","tmin_aug",
                         "tmin_sep","tmin_oct","tmin_nov","tmin_dec")
  names(tmax_stack2) = c("tmax_jan","tmax_feb","tmax_mar","tmax_apr",
                         "tmax_may","tmax_jun","tmax_jul","tmax_aug",
                         "tmax_sep","tmax_oct","tmax_nov","tmax_dec")
  names(bioc_stack2) = c("bio1","bio2","bio3","bio4","bio5","bio6","bio7",
                         "bio8","bio9","bio10","bio11","bio12","bio13","bio14",
                         "bio15","bio16","bio17","bio18","bio19")
  pred_stack2 <- stack(prec_stack2,
                       tmin_stack2,
                       tmax_stack2,
                       bioc_stack2,
                       longitude_raster,
                       latitude_raster,
                       elev2)
  names(pred_stack2)
  
  jan_r_pred2 <- predict(pred_stack2, model_jan)
  feb_r_pred2 <- predict(pred_stack2, model_feb)
  mar_r_pred2 <- predict(pred_stack2, model_mar)
  apr_r_pred2 <- predict(pred_stack2, model_apr)
  may_r_pred2 <- predict(pred_stack2, model_may)
  jun_r_pred2 <- predict(pred_stack2, model_jun)
  jul_r_pred2 <- predict(pred_stack2, model_jul)
  aug_r_pred2 <- predict(pred_stack2, model_aug)
  sep_r_pred2 <- predict(pred_stack2, model_sep)
  oct_r_pred2 <- predict(pred_stack2, model_oct)
  nov_r_pred2 <- predict(pred_stack2, model_nov)
  dec_r_pred2 <- predict(pred_stack2, model_dec)
  
  jan_r_pred2 <- mask(jan_r_pred2, gbr)
  feb_r_pred2 <- mask(feb_r_pred2, gbr)
  mar_r_pred2 <- mask(mar_r_pred2, gbr)
  apr_r_pred2 <- mask(apr_r_pred2, gbr)
  may_r_pred2 <- mask(may_r_pred2, gbr)
  jun_r_pred2 <- mask(jun_r_pred2, gbr)
  jul_r_pred2 <- mask(jul_r_pred2, gbr)
  aug_r_pred2 <- mask(aug_r_pred2, gbr)
  sep_r_pred2 <- mask(sep_r_pred2, gbr)
  oct_r_pred2 <- mask(oct_r_pred2, gbr)
  nov_r_pred2 <- mask(nov_r_pred2, gbr)
  dec_r_pred2 <- mask(dec_r_pred2, gbr)
  
  fp = file.path("~/magna_rusle/reports/output_data/r_factor_predictions/monthly_maps",
                 years)
  setwd(fp)
  fname = paste("future_r",years,ssp,sep="_")
  writeRaster(jan_r_pred2,filename=paste("jan",fname,sep="_"),format="GTiff",overwrite=TRUE)
  writeRaster(feb_r_pred2,filename=paste("feb",fname,sep="_"),format="GTiff",overwrite=TRUE)
  writeRaster(mar_r_pred2,filename=paste("mar",fname,sep="_"),format="GTiff",overwrite=TRUE)
  writeRaster(apr_r_pred2,filename=paste("apr",fname,sep="_"),format="GTiff",overwrite=TRUE)
  writeRaster(may_r_pred2,filename=paste("may",fname,sep="_"),format="GTiff",overwrite=TRUE)
  writeRaster(jun_r_pred2,filename=paste("jun",fname,sep="_"),format="GTiff",overwrite=TRUE)
  writeRaster(jul_r_pred2,filename=paste("jul",fname,sep="_"),format="GTiff",overwrite=TRUE)
  writeRaster(aug_r_pred2,filename=paste("aug",fname,sep="_"),format="GTiff",overwrite=TRUE)
  writeRaster(sep_r_pred2,filename=paste("sep",fname,sep="_"),format="GTiff",overwrite=TRUE)
  writeRaster(oct_r_pred2,filename=paste("oct",fname,sep="_"),format="GTiff",overwrite=TRUE)
  writeRaster(nov_r_pred2,filename=paste("nov",fname,sep="_"),format="GTiff",overwrite=TRUE)
  writeRaster(dec_r_pred2,filename=paste("dec",fname,sep="_"),format="GTiff",overwrite=TRUE)
  
  jan_magna2 = crop(jan_r_pred2,fort_boundary)
  feb_magna2 = crop(feb_r_pred2,fort_boundary)
  mar_magna2 = crop(mar_r_pred2,fort_boundary)
  apr_magna2 = crop(apr_r_pred2,fort_boundary)
  may_magna2 = crop(may_r_pred2,fort_boundary)
  jun_magna2 = crop(jun_r_pred2,fort_boundary)
  jul_magna2 = crop(jul_r_pred2,fort_boundary)
  aug_magna2 = crop(aug_r_pred2,fort_boundary)
  sep_magna2 = crop(sep_r_pred2,fort_boundary)
  oct_magna2 = crop(oct_r_pred2,fort_boundary)
  nov_magna2 = crop(nov_r_pred2,fort_boundary)
  dec_magna2 = crop(dec_r_pred2,fort_boundary)
  
  mean_jan2 = cellStats(jan_magna2, stat='mean', na.rm=TRUE)
  mean_feb2 = cellStats(feb_magna2, stat='mean', na.rm=TRUE)
  mean_mar2 = cellStats(mar_magna2, stat='mean', na.rm=TRUE)
  mean_apr2 = cellStats(apr_magna2, stat='mean', na.rm=TRUE)
  mean_may2 = cellStats(may_magna2, stat='mean', na.rm=TRUE)
  mean_jun2 = cellStats(jun_magna2, stat='mean', na.rm=TRUE)
  mean_jul2 = cellStats(jul_magna2, stat='mean', na.rm=TRUE)
  mean_aug2 = cellStats(aug_magna2, stat='mean', na.rm=TRUE)
  mean_sep2 = cellStats(sep_magna2, stat='mean', na.rm=TRUE)
  mean_oct2 = cellStats(oct_magna2, stat='mean', na.rm=TRUE)
  mean_nov2 = cellStats(nov_magna2, stat='mean', na.rm=TRUE)
  mean_dec2 = cellStats(dec_magna2, stat='mean', na.rm=TRUE)
  
  new_row = c(years,ssp,mean_jan2,mean_feb2,mean_mar2,
              mean_apr2,mean_may2,mean_jun2,mean_jul2,
              mean_aug2,mean_sep2,mean_oct2,mean_nov2,mean_dec2)
  return(new_row)
}


# calculate predictions for every ssp and all years
for (i in 1:length(years_list)) {
  for (j in 1:length(ssp_list)) {
    magna_monthly_r[nrow(magna_monthly_r) +1,] = 
      predict_future_monthly_R(years_list[i],ssp_list[j])
  }
}

write.csv(magna_monthly_r, file ="~/magna_rusle/reports/output_data/r_factor_predictions/magna_monthly_r_factors.csv",
          row.names = FALSE)

