# Script to get the mean prec, tmin, tmax and calculate biovars (rasters) 
# for years 2000-2018 from worldclim data

library(rgdal)
library(raster)
library(sf)
library(maps)
library(mapdata) 
library(dismo)

# get extent of UK 
gbr <- raster::getData("GADM", country="GBR", level=1) 

# define an array of months
months = c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")

# get average tmin for 2000-2018
dir="~/magna_rusle/data/input/monthly_climate/tmin"
for (i in 1:length(months)){
  file_path = file.path(dir,months[i])
  setwd(file_path)
  tmin_files <- list.files(path=file_path)
  tmin_files
  tmin_stack = stack(tmin_files)
  tmin_stack = crop(tmin_stack,gbr)
  mean <- stackApply(tmin_stack, indices =  rep(1,nlayers(tmin_stack)), fun = "mean", na.rm = T)
  x <- paste("tmin_mean",months[i], sep="_")
  print(x)
  eval(call("<-", as.name(x), mean))
}

tmin_monthly = stack(tmin_mean_jan,tmin_mean_feb,tmin_mean_mar,tmin_mean_apr,
                     tmin_mean_may,tmin_mean_jun,tmin_mean_jul,tmin_mean_aug,
                     tmin_mean_sep,tmin_mean_oct,tmin_mean_nov,tmin_mean_dec)

# get average tmax for 2000-2018
dir="~/magna_rusle/data/input/monthly_climate/tmax"
for (i in 1:length(months)){
  file_path = file.path(dir,months[i])
  setwd(file_path)
  tmax_files <- list.files(path=file_path)
  tmax_files
  tmax_stack = stack(tmax_files)
  tmax_stack = crop(tmax_stack,gbr)
  mean <- stackApply(tmax_stack, indices =  rep(1,nlayers(tmax_stack)), fun = "mean", na.rm = T)
  x <- paste("tmax_mean",months[i], sep="_")
  print(x)
  eval(call("<-", as.name(x), mean))
}

tmax_monthly = stack(tmax_mean_jan,tmax_mean_feb,tmax_mean_mar,tmax_mean_apr,
                     tmax_mean_may,tmax_mean_jun,tmax_mean_jul,tmax_mean_aug,
                     tmax_mean_sep,tmax_mean_oct,tmax_mean_nov,tmax_mean_dec)

# get average prec for 2000-2018
dir="~/magna_rusle/data/input/monthly_climate/prec"
for (i in 1:length(months)){
  file_path = file.path(dir,months[i])
  setwd(file_path)
  prec_files <- list.files(path=file_path)
  prec_files
  prec_stack = stack(prec_files)
  prec_stack = crop(prec_stack,gbr)
  mean <- stackApply(prec_stack, indices =  rep(1,nlayers(prec_stack)), fun = "mean", na.rm = T)
  x <- paste("prec_mean",months[i], sep="_")
  print(x)
  eval(call("<-", as.name(x), mean))
}

prec_monthly = stack(prec_mean_jan,prec_mean_feb,prec_mean_mar,prec_mean_apr,
                     prec_mean_may,prec_mean_jun,prec_mean_jul,prec_mean_aug,
                     prec_mean_sep,prec_mean_oct,prec_mean_nov,prec_mean_dec)


# calculate biovars for 2000-2018
biovars_current = biovars(prec_monthly,tmin_monthly,tmax_monthly)

# save raster stacks
writeRaster(biovars_current, filename="biovars_current.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(prec_monthly, filename="prec_monthly.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(tmin_monthly, filename="tmin_monthly.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
writeRaster(tmax_monthly, filename="tmax_monthly.tif", options="INTERLEAVE=BAND", overwrite=TRUE)


# or save individual raster layers
dir = "~/magna_rusle/reports/output_data/climate_means"
setwd(dir)

for (i in 1:nlayers(biovars_current)) {
  band = subset(biovars_current,subset=i)
  x=paste("bio",i,".tif",sep="")
  writeRaster(band,filename=x,format="GTiff",overwrite=TRUE)
}

for (i in 1:nlayers(prec_monthly)) {
  band = subset(prec_monthly,subset=i)
  x=paste("prec_",i,".tif",sep="")
  writeRaster(band,filename=x,format="GTiff",overwrite=TRUE)
}

for (i in 1:nlayers(tmin_monthly)) {
  band = subset(tmin_monthly,subset=i)
  x=paste("tmin_",i,".tif",sep="")
  writeRaster(band,filename=x,format="GTiff",overwrite=TRUE)
}

for (i in 1:nlayers(tmax_monthly)) {
  band = subset(tmax_monthly,subset=i)
  x=paste("tmax_",i,".tif",sep="")
  writeRaster(band,filename=x,format="GTiff",overwrite=TRUE)
}


