# _common.R - Global parameters and libraries

#------------------------------------------------------
# Set the report year
#------------------------------------------------------
RptYr <- 2025


#------------------------------------------------------
# Install and load libraries
#------------------------------------------------------
# Install pacman and pak for easier library workflow
install.packages(setdiff(c('pacman', 'pak'), rownames(installed.packages())))
# Install github packages
#pak::pak("teunbrand/elementalist")
#pak::pak("nmfs-ost/nmfspalette")

pacman::p_load_gh("teunbrand/elementalist",
                  "nmfs-ost/nmfspalette")

# Install and load all shared packages
pacman::p_load(tidyverse,
               lubridate,
               here,
               grid,
               gridExtra,
               stringr,
               geomtextpath, 
               elementalist,
               nmfspalette,
               scales,
               terra,
               tidyterra,
               metR,
               sf,
               rnaturalearth)


#------------------------------------------------------
# Set color palettes for the chapter
#------------------------------------------------------
oceans <- nmfs_palette("oceans")(3) # 1 = report_year, 3 = previous years
crustacean <- nmfs_palette("crustacean")(4) # 1 = linear trend
coral <- nmfs_palette("coral")(3) # 3 = annual average for Rpt Yr
ann_grey <- "#A5AAAF" # annual means; in NMFS branding guide as dove gray tint 1
waves <- nmfs_palette("waves")(3) # annual means; in NMFS branding guide but not in package
seagrass <- nmfs_palette("seagrass")(3)


#------------------------------------------------------
# Define ggplot theme
#------------------------------------------------------
plotTheme <- theme_minimal() +
  theme(panel.border = element_rect(color='#F1F3F3'),
        panel.grid.major = element_line(color='#F1F3F3'),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        axis.title = element_text(size=14, face="bold"), 
        axis.text = element_text(size=12),
        legend.position = 'bottom', 
        legend.background = element_rect_round(fill = 'gray90', color='transparent', radius = 0.3),
        legend.text = element_text(size=12), 
        legend.key.spacing.y = unit(-0.13,'cm'), 
        legend.key.spacing.x = unit(1.5,'cm'))


#------------------------------------------------------
# Map plotting function
#------------------------------------------------------
plot_chl_map <- function(
    currYr_file,
    climo_file,
    fill_scale,
    show_llbox = TRUE,
    contour_breaks = c(0.1, 0.2),
    min_contour_length = 150,
    proj_crs = "+proj=eck4 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
    # Set manual domain limits if desired; otherwise automatically pulled from raster
    xmin = NULL,
    xmax = NULL,
    ymin = NULL,
    ymax = NULL
) {
  # 1. Calculate Anomaly
  anom <- currYr_file - climo_file
  
  # 2. Determine Domain (Automatic from Raster vs Manual Input)
  r_ext <- ext(currYr)
  dom_xmin <- if (is.null(xmin)) r_ext[1] else xmin
  dom_xmax <- if (is.null(xmax)) r_ext[2] else xmax
  dom_ymin <- if (is.null(ymin)) r_ext[3] else ymin
  dom_ymax <- if (is.null(ymax)) r_ext[4] else ymax
  
  # 3. Reproject Rasters
  anom_proj <- project(anom, proj_crs)
  curr_proj <- project(currYr, proj_crs)
  
  # 4. Create Curved Map Outline from Domain Extent
  map_outline_sf <- st_polygon(list(matrix(c(
    dom_xmin, dom_ymin,  
    dom_xmax, dom_ymin,  
    dom_xmax, dom_ymax,  
    dom_xmin, dom_ymax,  
    dom_xmin, dom_ymin   
  ), ncol = 2, byrow = TRUE))) %>% 
    st_sfc(crs = 4326) %>% 
    st_set_crs(NA) %>% 
    st_segmentize(dfMaxLength = 0.5) %>% 
    st_set_crs(4326) %>% 
    st_sf()
  
  # 5. Land Preparation & Clipping
  world_land <- ne_countries(scale = "medium", returnclass = "sf") %>% 
    st_break_antimeridian(lon_0 = 180) %>% 
    st_make_valid()
  
  map_outline_projected <- st_transform(map_outline_sf, proj_crs)
  world_land_projected  <- st_transform(world_land, proj_crs)
  world_land_clipped   <- st_intersection(world_land_projected, map_outline_projected)
  
  # 6. Dynamic Axis Labels (Generated based on Domain Extent)
  lat_ticks <- pretty(c(dom_ymin, dom_ymax), n = 5)
  lat_ticks <- lat_ticks[lat_ticks >= dom_ymin & lat_ticks <= dom_ymax]
  
  lon_ticks <- pretty(c(dom_xmin, dom_xmax), n = 5)
  lon_ticks <- lon_ticks[lon_ticks >= dom_xmin & lon_ticks <= dom_xmax]
  
  labels_left <- st_as_sf(data.frame(
    lon = dom_xmin - 3.25, 
    lat = lat_ticks, 
    text = ifelse(lat_ticks == 0, "0°", ifelse(lat_ticks > 0, paste0(lat_ticks, "°N"), paste0(abs(lat_ticks), "°S")))
  ), coords = c("lon", "lat"), crs = 4326)
  
  labels_bottom <- st_as_sf(data.frame(
    lon = lon_ticks, 
    lat = dom_ymin - 2.25, 
    text = ifelse(lon_ticks > 180, paste0(360 - lon_ticks, "°W"),
                  ifelse(lon_ticks == 180, "180°",
                         ifelse(lon_ticks == 0, "0°", paste0(lon_ticks, "°E"))))
  ), coords = c("lon", "lat"), crs = 4326)
  
  # 7. Convert raster to data frame for text contours
  curr_df <- as_tibble(curr_proj, xy = TRUE)
  val_col <- names(curr_df)[3]
  
  # 8. Build Base Plot
  p_map <- ggplot() +
    geom_spatraster(data = anom_proj) +
    geom_spatraster_contour(
      data = curr_proj, 
      breaks = contour_breaks, 
      show.legend = FALSE, 
      linewidth = 0.5, 
      linetype = "dashed"
    ) +   
    geom_text_contour(
      data = curr_df,
      aes(x = x, y = y, z = .data[[val_col]]),
      breaks = contour_breaks,
      color = "gray20",
      size = 5,
      min.size = min_contour_length,
      skip = 0,
      label.placer = label_placer_fraction(frac = 0.5)
    ) +
    fill_scale +
    guides(
      fill = guide_colorbar(
        direction = "horizontal",
        barwidth = unit(0.7, "npc"),
        barheight = unit(10, "pt"),
        title.position = "top",
        title.hjust = 0.5,
        label.position = "top",
        label.hjust = 0.5
      )
    ) +
    geom_sf(data = world_land_clipped, fill = "gray80", color = "gray50", linewidth = 0.2) +
    geom_sf(data = map_outline_projected, fill = NA, color = "gray50", linewidth = 0.5)
  
  # 9. Conditionally Add Longline Box
  if (show_llbox) {
    box_sf <- st_polygon(list(matrix(c(
      180, 15, 240, 15, 240, 45, 180, 45, 180, 15
    ), ncol = 2, byrow = TRUE))) %>%
      st_sfc(crs = 4326)
    
    LLbox <- box_sf %>%
      st_set_crs(NA) %>%
      st_segmentize(dfMaxLength = 0.5) %>%
      st_set_crs(4326) %>%
      st_sf()
    
    p_map <- p_map +
      geom_sf(data = LLbox, fill = NA, color = "black", linewidth = 0.75) +
      geom_sf_text(
        data = st_as_sf(data.frame(lon = 180, lat = 45), coords = c("lon", "lat"), crs = 4326),
        aes(label = "Longline fishing grounds"),
        size = 6, vjust = -0.5, hjust = 0, fontface = "bold", color = "black"
      )
  }
  
  # 10. Map Furniture & Theme
  p_map <- p_map +
    geom_sf_text(data = labels_left, aes(label = text), size = 5, hjust = 1, color = "gray30") +
    geom_sf_text(data = labels_bottom, aes(label = text), size = 5, vjust = 1, color = "gray30") +
    coord_sf(
      crs = proj_crs,
      xlim = c(ext(curr_proj)[1], ext(curr_proj)[2]),
      ylim = c(ext(curr_proj)[3], ext(curr_proj)[4]), 
      expand = FALSE, 
      clip = "off"
    ) + 
    theme_minimal() + 
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(), 
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.position = "top",
      legend.justification = c(0.42, 0.5),
      legend.box.margin = margin(b = 0.5),
      plot.margin = margin(10, 10, 30, 40)
    ) +
    xlab(NULL) + 
    ylab(NULL)
  
  return(p_map)
}
