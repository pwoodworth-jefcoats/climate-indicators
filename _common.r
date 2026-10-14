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




