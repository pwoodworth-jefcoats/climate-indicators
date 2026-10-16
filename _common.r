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
               scales)


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
