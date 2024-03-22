library(ggplot2)
library(dplyr)
library(plotly)
library(raster)
library(rnaturalearth)
library(sf)
library(scales)
library(nmfspalette)
library(terra)

#function to calculate anomaly, recenter data over Pacific, crop to size, downsample, and trim whitespace
prep_map_data <- function(indicator, raster_ann, raster_clim, min_x, max_x, min_y, max_y, target_res){
  
  #calculate anomaly raster
  raster_anom <- raster_ann - raster_clim
  
  #CRS to recenter for plotting
  plot_crs <- "+proj=longlat +x_0=0 +y_0=0 +lat_0=0 +lon_0=180 +datum=WGS84 +no_defs"
  
  #split into two along antimeridian
  raster_w <- raster_ann %>%
    raster::crop(c(xmin = min_x, xmax = 180, ymin = min_y, ymax = max_y))
  
  raster_e <- raster_ann %>%
    raster::crop(c(xmin = 180, xmax = max_x, ymin = min_y, ymax = max_y)) %>%
    raster::rotate(.)
  
  #repeat for anomaly 
  raster_anom_w <- raster_anom %>%
    raster::crop(c(xmin = min_x, xmax = 180, ymin = min_y, ymax = max_y)) 
  
  raster_anom_e <- raster_anom %>%
    raster::crop(c(xmin = 180, xmax = max_x, ymin = min_y, ymax = max_y)) %>%
    raster::rotate(.) 
  
  #downsample if needed, otherwise just reproject
  if (res(raster_ann)[1] < target_res){
    raster_w <- raster_w %>% projectRaster(., crs = plot_crs, res = target_res)
    raster_e <- raster_e %>% projectRaster(., crs = plot_crs, res = target_res)
    raster_anom_w <- raster_anom_w %>% projectRaster(., crs = plot_crs, res = target_res)
    raster_anom_e <- raster_anom_e %>% projectRaster(., crs = plot_crs, res = target_res)
  } else{
    raster_w <- raster_w %>% projectRaster(., crs = plot_crs)
    raster_e <- raster_e %>% projectRaster(., crs = plot_crs)
    raster_anom_w <- raster_anom_w %>% projectRaster(., crs = plot_crs)
    raster_anom_e <- raster_anom_e %>% projectRaster(., crs = plot_crs)
  }
  
  #merge back together
  raster_proj <- raster::merge(raster_w, raster_e, tolerance = 0.5) %>% 
    raster::trim(.) #remove whitespace
  
  #repeat for anomaly map
  raster_anom_proj <- raster::merge(raster_anom_w, raster_anom_e, tolerance = 0.5) %>% 
    raster::trim(.) #remove whitespace
  
  #convert to points for plotting
  raster_df <- raster_proj %>%
    rasterToPoints(.) %>% #convert to points   
    as.data.frame() %>% #convert to dataframe
    mutate(x_disp = ifelse(x <= 0, x + 180, -(180 - x))) %>%
    mutate(labels = paste0(signif(x_disp, 3), ", ", signif(y, 3),
                           "<br>", indicator,
                           ": ", signif(layer, 3))) %>%
    mutate(ID = indicator)
  
  #repeat for anomaly data
  raster_anom_df <- raster_anom_proj %>%
    rasterToPoints(.) %>% #convert to points   
    as.data.frame() %>% #convert to dataframe
    mutate(x_disp = ifelse(x <= 0, x + 180, -(180 - x))) %>% #recenter at 180
    mutate(layer_disp = ifelse(layer > 0, paste0("+ ", layer), layer)) %>% #format hovertext
    mutate(layer_disp = ifelse(layer == 0, 0, layer)) %>%
    mutate(labels = paste0(signif(x_disp, 3), ", ", signif(y, 3),
                           "<br>", indicator,
                           " Anomaly: ", signif(layer, 3))) %>%
    mutate(ID = paste0(indicator, "_anom"))
  
  #combine datasets
  raster_df <- bind_rows(raster_df, raster_anom_df)
  
  #return data
  return(list(raster_proj, raster_anom_proj, raster_df))
  
}

#read in raster data
tatd_2022 <- raster("T_at_200300_yr2022.nc")
tatd_clim <- raster("T_at_200300_climo_1980thru2021.nc")
sst_2022 <- raster("CRW_sst_v3_1_2022-clim_9b86_1fd9_bbd5_U1699483388394.nc")
sst_clim <- raster("CRW_sst_v3_1_1985-2021-clim_437a_63a6_52e7_U1699483412096.nc")
chl_2022 <- raster("esa-cci-chla-2022-clim_v6-0_6f24_0a8d_aa83_U1699486681824.nc",
                   varname = "chlor_a")
chl_clim <- raster("esa-cci-chla-1998-2021-clim-v6-0_5cd0_57c8_f0c8_U1699486666335.nc",
                   varname = "chlor_a")
md50_2022 <- raster("md50_exp-2022-clim_61be_da6b_0507_U1699492344974.nc")
md50_clim <- raster("md50_exp-1998-2021-clim_61be_da6b_0507_U1699492333044.nc")

#get lat coordinates from raw data
original_bbox <- bbox(tatd_2022)
min_x <- original_bbox[1,1] 
max_x <- original_bbox[1,2]
min_y <- original_bbox[2,1]
max_y <- original_bbox[2,2]

#set target resolution in degrees
target_res <- 0.5

#get prepped data 
#you will get warnings here - they are ignoreable!
tatd_out <- prep_map_data("TatD", tatd_2022, tatd_clim, min_x, max_x,
                          min_y, max_y, target_res)
sst_out <- prep_map_data("SST", sst_2022, sst_clim, min_x, max_x,
                         min_y, max_y, target_res)
chl_out <- prep_map_data("Chl", chl_2022, chl_clim, min_x, max_x,
                         min_y, max_y, target_res)
md50_out <- prep_map_data("MD50", md50_2022, md50_clim, min_x, max_x,
                          min_y, max_y = 45, target_res) #note - different max y 

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
write.csv(raster_df_all, "Indicator_Dashboard/Data/Dashboard_Map_Data_2022.csv")

#get lon coordinates from cropped data
cropped_bbox <- bbox(tatd_out[[1]])
min_x <- cropped_bbox[1,1] 
max_x <- cropped_bbox[1,2]

#read in coast data
world_coasts <- rnaturalearth::ne_coastline(scale = 50, returnclass = "sf") %>% st_make_valid()
world_countries <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf") %>% st_make_valid()

#CRS to recenter for plotting - warnings here ok
plot_crs <- "+proj=longlat +x_0=0 +y_0=0 +lat_0=0 +lon_0=180 +datum=WGS84 +no_defs"

#base maps
land_proj <- world_countries %>% 
  st_break_antimeridian(lon_0 = 180) %>% #trim out dateline to avoid artifacts
  st_transform(crs = plot_crs) %>% #adjust CRS for Pacific
  st_crop(x = ., y = c(xmin = min_x, ymin = min_y, 
                       xmax = max_x, ymax = max_y)) %>% #crop to area
  st_cast(., "MULTIPOLYGON") #recast for plotting

coast_proj <- world_coasts %>% 
  st_break_antimeridian(lon_0 = 180) %>% #trim out dateline to avoid artifacts
  st_transform(crs = plot_crs) %>% #adjust CRS for Pacific
  st_crop(x = ., y = c(xmin = min_x, ymin = min_y, 
                       xmax = max_x, ymax = max_y)) %>% #crop to area
  st_cast(., "MULTILINESTRING") #recast for plotting

#base maps for ENSO
land_proj_enso <- world_countries %>% 
  st_break_antimeridian(lon_0 = 180) %>% #trim out dateline to avoid artifacts
  st_transform(crs = plot_crs) %>% #adjust CRS for Pacific
  st_crop(x = ., y = c(xmin = min_x, ymin = -5, 
                       xmax = max_x, ymax = max_y)) %>% #crop to area
  st_cast(., "MULTIPOLYGON") #recast for plotting

coast_proj_enso <- world_coasts %>% 
  st_break_antimeridian(lon_0 = 180) %>% #trim out dateline to avoid artifacts
  st_transform(crs = plot_crs) %>% #adjust CRS for Pacific
  st_crop(x = ., y = c(xmin = min_x, ymin = -5, 
                       xmax = max_x, ymax = max_y)) %>% #crop to area
  st_cast(., "MULTILINESTRING") #recast for plotting

#base maps for Tropical Cyclones
land_proj_tcs <- world_countries %>% 
  st_break_antimeridian(lon_0 = 180) %>% #trim out dateline to avoid artifacts
  st_transform(crs = plot_crs) %>% #adjust CRS for Pacific
  st_crop(x = ., y = c(xmin = 90, ymin = -50, 
                       xmax = 280, ymax = 50)) %>% #crop to area
  st_cast(., "MULTIPOLYGON") #recast for plotting

coast_proj_tcs <- world_coasts %>% 
  st_break_antimeridian(lon_0 = 180) %>% #trim out dateline to avoid artifacts
  st_transform(crs = plot_crs) %>% #adjust CRS for Pacific
  st_crop(x = ., y = c(xmin = 90, ymin = -50, 
                       xmax = 280, ymax = 50)) %>% #crop to area
  st_cast(., "MULTILINESTRING") #recast for plotting

#quick test
ggplot() + geom_sf(data = land_proj, color = "grey", fill = "grey") +
  geom_sf(data = coast_proj, color = "black")

ggplot() + geom_sf(data = land_proj_enso, color = "grey", fill = "grey") +
  geom_sf(data = coast_proj_enso, color = "black")

ggplot() + geom_sf(data = land_proj_tcs, color = "grey", fill = "grey") +
  geom_sf(data = coast_proj_tcs, color = "black")

#save data
saveRDS(land_proj, "Indicator_Dashboard/Data/rnatearth_land.RData")
saveRDS(coast_proj, "Indicator_Dashboard/Data/rnatearth_coast.RData")
 
saveRDS(land_proj_enso, "Indicator_Dashboard/Data/rnatearth_enso_land.RData")
saveRDS(coast_proj_enso, "Indicator_Dashboard/Data/rnatearth_enso_coast.RData")

saveRDS(land_proj_tcs, "Indicator_Dashboard/Data/rnatearth_tcs_land.RData")
saveRDS(coast_proj_tcs, "Indicator_Dashboard/Data/rnatearth_tcs_coast.RData")
