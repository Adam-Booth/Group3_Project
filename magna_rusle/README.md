# magna_rusle

This folder contains the ProjectTemplate directory for calculating RUSLE paramter Rainfall erosivity for UK and at the site of Magna Fort using R. Source codes can be found in 'src' folder. Input data can be found in 'data' folder. Various input data have not been included within the ProjectTemplate directory mainly due to size of files. However, detailed instructions can be found on how to obtain the data and where to store within the directory in order to run the scripts provided. In addition, any outputs from the scripts have been saved under "reports\output_data" to allow these to be checked and compared to initial results and to facilitate running of later scripts.  
There are 3 scripts involved in this section of the project, which are meant to be ran in the corresponding orders as they will use ouputs from previous scripts. However, outputs have been stored, as mentioned above, so not all scripts need to be ran. (For example, script 1 - rainfall_erosivity_calculation.R takes over an hour to run and thus providing output data allows this to be avoided.)  

## 1. Script 1: 'rainfall_erosivity_calculation.R'

**Purpose:** Script calculates station-wise annual and monthly rainfall erosivity based on high temporal precipitation records. 

**Input data to include:**

* 'station_names.csv' - list of precipitation stations in the UK including information such as name, location and elevation. This file is already included within the directory (obtained openly from CEDA catelogue).
* 'prec_stations'- folder containing individual folders for each precipitation station included to calculate rainfall erosivity. In each station's folder csv files for each year of hourly rainfall data available should be found. These have not been included within this directory due to file size. However, they can be openly downloaded from the CEDA catelogue after registration from the following link: https://catalogue.ceda.ac.uk/uuid/77187ac1e0a341ca993c3366f8c59c3c. Years will only be included within the calculations if they have under 10% missing hourly data - the code will filter for this criteria automatically. 

**Output data created:** (can be found in '~/magna_rusle/reports/output_data')

* 'station_missing_data' - csv for each station detailing percentage of missing data in each year of hourly precipitation records contained within input files
* 'monthly_station_r_factors.csv' - csv containing data on annual and monthly rainfall erosivitys of each precipitation station (also include number of years of hourly records used in the calculations)
* 'station_past_erosivity_1995-2019.csv'- csv containing data on the past yearly erosivitys for each station

# 2. Script 2: 'mean_climate_vars.R'

**Purpose:** Script calculates mean monthly precipitation, minimum and maximum monthly temperatures and bioclimatic variables (see https://www.worldclim.org/data/bioclim.html), from WorldClim raster data for the UK over years of data included. 

**Input data to include:**

* 'monthly_climate' - contains folders for precipitation ('prec'), minimum and maximum temperatures ('tmin' and 'tmax'). In each folder there are 12 folders - one for each month. In these folders raster datasets for the corresponding monthly climate data should be obtained from WorldClim via the link: https://www.worldclim.org/data/monthlywth.html#. In this project years 2000-2018 were included. (i.e. there should then be 19 files in each subfolder i.e. variable/month - this script is used to take an average over how many years of data are included.) 

**Output data created:** (can be found in '~/magna_rusle/reports/output_data')

* 'climate_means' - folder contains mean monthly climate variables (prec, tmin, tmax and biovars) for years of WorldClim data contained in input, for the UK. 

## 3. Script 3: 'spatial_regression.R'

**Purpose:** Script interpolates station R factors (and monthly R factors) to create a map of UK r factors and extracts these values for the area containing Magna Fort. It then goes on to use the same regression model used to interpolate R factor to predict future r factors based on future climate projections.

**Input data to include:**

* 'future_climate' - future climate data (2.5min spatial resolution) obtained from WorldClim via the link: (https://www.worldclim.org/data/cmip6/cmip6_clim2.5m.html). In this project projections from CNRM-CM6-1 GCM were used however the script can be edited to use other models. In each subdirectory (years/ssp) there should be 4 files corresponding to tmin, tmax, prec and bioc obtained for that yearly period and ssp from WorldClim. These have not been included due to large file sizes.
* 'magna_fort_boundary.shp' - shapefile of Magna Fort boundary used to extract mean value of R factors over this area (this was created in ArcGIS Pro and is already included in the directory.) 
* script 2 is required to be ran before this as 'prec/tmin/tmax_monthly' and 'biovars_current' are required in the global environment. Alternatively, these files have been saved as outputs from script 2 and an extra command can be added in to load them. 
* 'wc2.1_30s_elev.tif' - global elevation raster obtained from worldclim (30s spatial resolution).

**Output data created:** (can be found in '~/magna_rusle/reports/output_data/r_factor_predictions')

* 'R_factor_magna.csv' - csv file containing R factor predictions at Magna Fort for current and future years for different SSPs
* 'current_r_factor_pred.tif' - raster of UK current R factors
* 'future_r_factor_pred_[years]_[ssp].tif' - raster of UK predicted R factors for future years (and different SSPs)
* 'future_r_factor_pred_[years]_[ssp] relative_change.tif' - raster of relative change since baseline of UK predicted R factors for future years (and different SSPs)
* 'monthly_r_factors' - rasters of future monthly r factors 
* 'magna_monthly_r_factors.csv' - csv file containing data of current and future monthly R factors at Magna Fort. 
* 'model_cross_validation_stats.csv' - csv file containing cross validation stats for each model (yearly and 12 monthly models), includes R2, MAE and RMSE
