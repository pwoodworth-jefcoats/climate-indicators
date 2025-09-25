library(ggplot2)
library(dplyr)
library(plotly)
library(raster)
library(sf)
library(scales)
library(nmfspalette)
library(terra)
library(here)

#function to calculate anomaly, recenter data over Pacific, crop to size, downsample, and trim whitespace
prep_map_data <- function(indicator, raster_ann, raster_clim, min_x, max_x, min_y, max_y, target_res){
  
  #calculate anomaly raster
  raster_anom <- raster_ann - raster_clim
  names(raster_ann) <- 'layer'
  
  #convert to points for plotting
  raster_df <- raster_ann %>%
    rasterToPoints(.) %>% #convert to points   
    as.data.frame() %>% #convert to dataframe
    mutate(x_180 = ifelse(x > 180, x-360, x)) %>%
    rename(x_disp=x, x=x_180) %>%  
    mutate(labels = paste0(signif(x, 3), ", ", signif(y, 3),
                           "<br>", indicator,
                           ": ", signif(layer, 3))) %>%
    mutate(ID = indicator)
  
  #repeat for anomaly data
  raster_anom_df <- raster_anom %>%
    rasterToPoints(.) %>% #convert to points   
    as.data.frame() %>% #convert to dataframe
    mutate(x_180 = ifelse(x > 180, x-360, x)) %>%
    rename(x_disp=x, x=x_180) %>% 
    mutate(layer_disp = ifelse(layer > 0, paste0("+ ", layer), layer)) %>% #format hovertext
    mutate(layer_disp = ifelse(layer == 0, 0, layer)) %>%
    mutate(labels = paste0(signif(x, 3), ", ", signif(y, 3),
                           "<br>", indicator,
                           " Anomaly: ", signif(layer, 3))) %>%
    mutate(ID = paste0(indicator, "_anom"))
  
  #combine datasets
  raster_df <- bind_rows(raster_df, raster_anom_df)
  
  #return data
  return(list(raster_ann, raster_anom, raster_df))
  
}

#read in raster data
tatd_2024 <- raster(here("Temperature_at_Depth","T_at_200300_yr2024.nc"))
tatd_clim <- raster(here("Temperature_at_Depth", "T_at_200300_climo_1980thru2009.nc"))
sst_2024 <- raster(here("Sea_Surface_Temperature", "sst_yr2024_dashboard.nc"))
sst_clim <- raster(here("Sea_Surface_Temperature", "sst_climo_dashboard.nc"))
chl_2024 <- raster(here("Ocean_Color", "chl_yr2024_dashboard.nc")) 
chl_clim <- raster(here("Ocean_Color", "chl_climo_dashboard.nc")) 
md50_2024 <- raster(here("Median_Phytoplankton_Size", "medphyto_yr2024_dashboard.nc"))
md50_clim <- raster(here("Median_Phytoplankton_Size", "medphyto_climo_dashboard.nc"))

# get lat coordinates from raw data
# Different data sources have different bounds
reanalysis_bbox <- bbox(tatd_2024) 
satellite_bbox <- bbox(sst_2024)
r_min_x <- reanalysis_bbox[1,1] 
r_max_x <- reanalysis_bbox[1,2]
r_min_y <- reanalysis_bbox[2,1]
r_max_y <- reanalysis_bbox[2,2]
s_min_x <- satellite_bbox[1,1] 
s_max_x <- satellite_bbox[1,2]
s_min_y <- satellite_bbox[2,1]
s_max_y <- satellite_bbox[2,2]

#set target resolution in degrees
target_res <- 0.5

# Workaround for missing information
srs <- tatd_2024@srs
sst_2024@srs <- srs
chl_2024@srs <- srs
md50_2024@srs <- srs

#get prepped data 
#you will get warnings here - they are ignoreable!
tatd_out <- prep_map_data("TatD", tatd_2024, tatd_clim, r_min_x, r_max_x,
                          r_min_y, r_max_y, target_res)
sst_out <- prep_map_data("SST", sst_2024, sst_clim, s_min_x, s_max_x,
                         s_min_y, s_max_y, target_res)
chl_out <- prep_map_data("Chl", chl_2024, chl_clim, s_min_x, s_max_x,
                         s_min_y, s_max_y, target_res)
md50_out <- prep_map_data("MD50", md50_2024, md50_clim, s_min_x, s_max_x,
                          s_min_y, s_max_y, target_res) #note - different max y 

#quick tests - annual maps
plot(tatd_out[[1]])
plot(sst_out[[1]])
plot(chl_out[[1]])
plot(md50_out[[1]])

#quick tests - anomaly maps
plot(tatd_out[[2]])
plot(sst_out[[2]])
plot(chl_out[[2]])
plot(md50_out[[2]])

#combine data and write to file
raster_df_all <- bind_rows(tatd_out[[3]], sst_out[[3]], chl_out[[3]], md50_out[[3]])
write.csv(raster_df_all, here("Indicator_Dashboard", "Data", "Dashboard_Map_Data_2024.csv"))
# write.csv(tatd_out[[3]], here("Temperature_at_Depth", "TempAtDepth_map_data_2024.csv"))
# write.csv(sst_out[[3]], here("Sea_Surface_Temperature", "SST_map_data_2024.csv"))
# write.csv(chl_out[[3]], here("Ocean_Color", "Chl_map_data_2024.csv"))
# write.csv(md50_out[[3]], here("Median_Phytoplankton_Size", "Median_Phyto_map_data_2024.csv"))