# Script to make a gif with TZCF and STF for SAFE report
# EDIT: 18 APR 2024 updated the climatology range to 1985:2009 for STF and 1998-2009 for TZCF

# Load libraries
library(raster)
library(rasterVis)
library(mapdata)
library(maptools)
library(cmocean)
library(RColorBrewer)
library(latticeExtra)
library(grid)
library(rerddap)
library(terra)

# Set Last year to include
year=2023

# Set working directory
setwd('~/Documents/SAFEindicators/RAnalysis_SAFE/Output/')
if (!dir.exists(file.path('STFfigs'))) { dir.create(file.path('STFfigs')) }
if (!dir.exists(file.path('TZCFfigs'))) { dir.create(file.path('TZCFfigs')) }

# Area for which to make the maps
lat1 <- 45
lat2 <- 5
lon1 <- 180
lon2 <- 240

#----------------------------------------------------------
#   STF  
#----------------------------------------------------------
## Download data from erddap
# For temp we are using the Coral Reef Watch product https://oceanwatch.pifsc.noaa.gov/erddap/griddap/CRW_sst_v3_1_monthly.html
# for ease I'm downloading one month at a time for all years, so first Jan then Feb and then March
sstInfo <- info(url = 'https://oceanwatch.pifsc.noaa.gov/erddap', datasetid = 'CRW_sst_v3_1_monthly')
# Comment below out if you have already downloaded and saved the .nc file
# Looping through to grab files for Jan, then Feb, then March
for (i in 1:3) {
  sstDat <- griddap(sstInfo, time = c(paste0('1985-0',i,'-28'),paste0(year, '-0',i,'-28')), latitude = c(lat1, lat2), longitude = c(lon1, lon2),
                    fields = sstInfo$variables$variable_name[1],
                    stride = c(12,1,1),
                    store=disk('STFdata'))
}

# Read in the files
STFfiles <- list.files('STFdata/', full.names = T)
# Read all three files into R and combine as stack and sort them after the date (name)
t1 <- brick(STFfiles[1])
t2 <- brick(STFfiles[2])
t3 <- brick(STFfiles[3])
# Make a climatology for all years except the most recent/current one
tclim <- mean(stack(t2[[1:25]], t1[[1:25]], t3[[1:25]]))

## Make the plot
# Get land information and make it into a spatial object
land <- maps::map('world2', fill=TRUE, xlim=c(180,240), ylim=c(20,45), plot=FALSE)
ids <- sapply(strsplit(land$names, ":"), function(x) x[1])
bPols <- map2SpatialPolygons(land, IDs=ids, proj4string=CRS('+proj=longlat +datum=WGS84 +no_defs'))

# make map themes
mapThemeSST <- rasterTheme(region=rev(brewer.pal(11, 'Spectral')))
mapThemeChl <- rasterTheme(region=cmocean('algae')(50))

# make averaging index
idx <- 1:nlayers(t1)
# years for labels
yr <- substr(names(t1), 2, 5) #Ideally you grab this from the names in the file, to make sure 
# Make blank dataframe for 'ghost' lines
datAll <- data.frame()

# Start loop to make plots with 18 degree C highlighted
for (i in idx) {
  # Put the same year into one raster brick and average
  dat <- mean(brick(t1[[i]],t2[[i]],t3[[i]]))
  # Turn the raster into points to isolate the 18 deg C line and save as dataframe for the ghost lines. I couldn't figure out how to get contourline to draw from multiple raster layers so this is my not so elegant fix
  datTemp <- rasterToPoints(dat)
  test <- data.frame(datTemp[near(datTemp[,3], 18, 0.01),], ID=i)
  # Sort every other year in increasing or decreasing long for plotting so you don't have lines crossing the map to connect the beginning of one line with the end of the previous. 
  test <- test[order(test$x, test$y, decreasing = (i %% 2 == 0)),]
  # Add this years line to last years lines
  datAll <- rbind(datAll, test) 
  # Grab the year you are plotting
  names(dat) <- yr[i]
  plot_name <- paste0('STFfigs/STF_', yr[i],'.png')
  png(plot_name, width=7, height=5, res=300, units='in') 
  print(levelplot(dat, pretty=T, margin=F, par.setting=mapThemeSST, at=seq(4,29.5,by=0.3)) + 
          layer(panel.lines(datAll, col='gray', lwd=0.5)) +
          contourplot(dat, pretty=T, margin=F, cuts=1, at=c(0,30,by=18), labels=F) + 
          contourplot(tclim, pretty=T, margin=F, cuts=1, at=c(0,30,by=18), labels=F, lty=2) + 
          layer(sp.polygons(bPols, fill='gray50', col='gray30')) + 
          layer(panel.text(x=183,y=44,labels=yr[i], size=10)))
  dev.off()
}

#----------------------------------------------------------
#   TCZF
#----------------------------------------------------------

## Download data from erddap
# For temp we are using the ESA OC CCI product https://oceanwatch.pifsc.noaa.gov/erddap/griddap/esa-cci-chla-monthly-v6-0.html
# for ease I'm downloading one month at a time for all years, so first Jan then Feb and then March
chlInfo <- info(url = 'https://oceanwatch.pifsc.noaa.gov/erddap', datasetid = 'esa-cci-chla-monthly-v6-0')
# Comment below out if you have already downloaded and saved the .nc file
# Looping through to grab files for Jan, then Feb, then March
for (i in 1:3) {
  chlDat <- griddap(chlInfo, time = c(paste0('1998-0',i, '-01'),paste0(year, '-0',i, '-01')), latitude = c(lat2, lat1), longitude = c(lon1, lon2),
                    fields = chlInfo$variables$variable_name[1],
                    stride = c(12,1,1),
                    store=disk('TZCFdata'))
}

# Read in the files
TZCFfiles <- list.files('TZCFdata/', full.names = T)
# Read all three files into R and combine as stack and sort them after the date (name)
c1 <- brick(TZCFfiles[1])
c2 <- brick(TZCFfiles[2])
c3 <- brick(TZCFfiles[3])
# Make a climatology for all years except the most recent/current one
cclim <- mean(stack(c1[[1:12]], c2[[1:12]], c3[[1:12]]), na.rm=T)
cclim <- aggregate(cclim, fact=7)
cclim_crop <- crop(cclim, extent(c(180, 240, 25, 45))) # to clean up the areas around the front and make figure less busy

## Make the plot
# make averaging index
idx <- 1:nlayers(c1)
# years for labels
yr <- substr(names(c1), 2, 5) #Ideally you grab this from the names in the file, to make sure 
# Make blank dataframe for 'ghost' lines
chlAll <- data.frame()

# Start loop to make plots with 18 degree C highlighted
for (i in idx) {
  # Put the same year into one raster brick and average
  dat <- mean(brick(c1[[i]],c2[[i]],c3[[i]]), na.rm=T)
  dat <- aggregate(dat, fact=7)
  # Grab the year you are plotting
  names(dat) <- yr[i]
  plot_name <- paste0('TZCFfigs/TZCF_coarse_', yr[i],'.png')
  png(plot_name, width=7, height=5, res=300, units='in') 
  print(levelplot(dat, pretty=T, margin=F, par.setting=mapThemeChl, at=seq(0,1,by=0.01)) + 
          contourplot(crop(dat, extent(c(180, 240, 27, 42))), pretty=T, margin=F, cuts=1, at=c(0,30,by=0.2), labels=F) + 
          contourplot(cclim_crop, pretty=T, margin=F, cuts=1, at=c(0,30,by=0.2), labels=F, lty=3) + 
          layer(sp.polygons(bPols, fill='gray50', col='gray30')) + 
          layer(panel.text(x=183,y=44,labels=yr[i], size=10)))
  dev.off()
}



