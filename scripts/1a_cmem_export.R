#set libraries
library(tidyverse)
library(here)
library(terra)

# get parameters to inform domain ####
bw_locs <- read.csv(here("data/locs/blwh_full_dataset_ps.csv")) #blue whale location data
ymin <- min(bw_locs$lat) - 2
ymax <- max(bw_locs$lat) + 2
xmin <- min(bw_locs$lon) - 2
xmax <- max(bw_locs$lon) + 2

# connect to full cmem catalog
path_copernicusmarine <- here("data/cmem_library/copernicusmarine.exe") 

# function to extract data from CMEM catalog ####
cmem_nc <- function(cmem_id, out_directory, out_file, cmem_var, start_date, end_date, lon_min = xmin, lat_min = ymin, lon_max = xmax, lat_max = ymax){

    command <- paste(
      shQuote(path_copernicusmarine),
      "subset",
      paste("--dataset-id", cmem_id),
      paste("--variable", cmem_var), 
      paste("--start-datetime", start_date),
      paste("--end-datetime", end_date),
      paste("--minimum-longitude", round(lon_min, 3)),
      paste("--maximum-longitude", round(lon_max, 3)),
      paste("--minimum-latitude", round(lat_min, 3)),
      paste("--maximum-latitude", round(lat_max, 3)),
      "--minimum-depth 0", #just downloads surface layer
      "--maximum-depth 1",
      "-o", shQuote(out_directory), #save directory
      paste("--output-filename", out_file), #save filename within directory
      sep = " "
    )

  system(command)
}

# ocean physics subset extract: 0m, full domain, Jan 2016 ####
phys_vars <- c("thetao", "so", "mlotst", "uo", "vo") #SST, salinity, MLD, eastward ocean current velocity, northward ocean current velocity

#loop through physics vars and downloads one .nc file per variable across the specified time window
for(i in 1:length(phys_vars)){
  curr_var <- phys_vars[i]
  print(curr_var)

  cmem_nc(cmem_id = "cmems_mod_glo_phy_my_0.083deg_P1D-m",
          out_directory = here("data/physics/test"), 
          out_file = paste0(curr_var,"_test.nc"),
          cmem_var = curr_var,
          start_date = "2016-01-01T00:00:00",
          end_date = "2016-02-01T00:00:00") 
  
  #regrid phys vars to 0.25 res
  curr_rast <- rast(here(paste0(out_directory,"/", out_file))) 

  template_rast <- rast(
      crs = crs(curr_rast),
      extent = ext(-166, -81, -10, 61), 
      resolution = 0.25 
    )

  rast_clean <- resample(curr_rast, template_rast)
  time(rast_clean) <- time(curr_rast)
  varnames(rast_clean) <- varnames(curr_rast)

  writeCDF(rast_clean, here(paste0("data/physics/test_clean/", curr_var, ".nc")))
}

# ocean biogeochem subset extract: 0m, full domain, Jan 2016 ####
biogeo_vars <- c("chl", "o2", "nppv") #chlorophyll, dissolved oxygen, net primary productivity

#loop through biogeochem vars 
for(i in 1:length(biogeo_vars)){
  curr_var <- biogeo_vars[i]
  print(curr_var)

  cmem_nc(cmem_id = "cmems_mod_glo_bgc_my_0.25deg_P1D-m",
          out_directory = here("data/biogeo/test"), 
          out_file = paste0(curr_var,"_test.nc"),
          cmem_var = curr_var,
          start_date = "2016-01-01T00:00:00",
          end_date = "2016-02-01T00:00:00") 
  
  #crop biogeo vars to match extent
  curr_rast <- rast(here(paste0(out_directory,"/", out_file))) 

  template_rast <- rast(
      crs = crs(curr_rast),
      extent = ext(-166, -81, -10, 61), 
      resolution = 0.25 
    )

  rast_clean <- resample(curr_rast, template_rast)
  time(rast_clean) <- time(curr_rast)
  varnames(rast_clean) <- varnames(curr_rast)

  writeCDF(rast_clean, here(paste0("data/biogeo/test_clean/", curr_var, ".nc")))
}

# static variables ####
#bathymetry
static_directory <- here("data/static")

command <- paste(
  shQuote(path_copernicusmarine),
  "subset",
  "--dataset-id cmems_mod_glo_bgc_my_0.25deg_static",
  "--dataset-part mask",
  "--variable deptho", #bathy (m)
  "--minimum-longitude -166.049",
  "--maximum-longitude -81.314",
  "--minimum-latitude -10.053",
  "--maximum-latitude 60.881",
  "-o", shQuote(static_directory), #save directory
  "--output-filename bathy.nc", #save filename within directory
  sep = " "
)

 #crop bathy vars to match extent
  curr_rast <- rast(here("data/static/bathy.nc")) 

  template_rast <- rast(
      crs = crs(curr_rast),
      extent = ext(-166, -81, -10, 61), 
      resolution = 0.25 
    )

  rast_clean <- resample(curr_rast, template_rast)
  varnames(rast_clean) <- varnames(curr_rast)

  writeCDF(rast_clean, here("data/static/bathy_clean.nc"))