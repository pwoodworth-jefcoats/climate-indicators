# Indicator script for dealer dataset
# edited 3/21/2017 to include full 2016 and to use departure year instead of landing year for better quarter patterns
# edited 18 APR 2024 to change the climatology range to 2000-2009 for the Weight Density Distribution
#--------------------------------------------------------
# Clear workspace
rm(list=ls())

# Set working directory
mainDir <- '~/Documents/SAFEindicators/RAnalysis_SAFE/'
setwd(mainDir)

#--------------------------------------------------------
# Load Libraries
library(ggplot2)
library(lubridate)
library(scales)
library(gridExtra)
library(dplyr)
library(grid)
library(png)

# Define the last year included in the report
#--------------------------------------------------------
yr <- 2023


# Load Data
#--------------------------------------------------------
setwd(file.path(mainDir, 'Data'))
dealer <- readRDS(paste0('SAFE_dealer_', yr, '.rds'))

#--------------------------------------------------------
# Process data
#--------------------------------------------------------
# Clean up data and calculate new variables
dealer <- dealer[which(dealer$SPECIES_FK != 0),]
dealer <- dealer[which(dealer$SPECIES_FK != 722),]  # Remove squid
dealer <- dealer[which(dealer$SPECIES_FK != 706),]  # Remove spiny lobster
dealer <- dealer[which(dealer$SPECIES_FK != 708),]  # Remove shrimp
dealer <- dealer[which(dealer$TRIP_NUM != 0),] # Remove all records without proper trip numbers
dealer <- dealer[which(dealer$NUM_SOLD > 0),] # Remove all records with <= 0 fish sold
dealer$WeightGG <- dealer$MEAN_SOLD_WEIGHT/2.204624  # convert form lbs to kg

# Convert date to useful format
#dealer$Date <- as.POSIXct(strptime(dealer$REPORT_DATE, '%d-%b-%y'))
dealer$Date <- as.POSIXct(dealer$REPORT_DATE)
dealer$Year <- year(dealer$Date)
dealer$Month <- month(dealer$Date)

# Limit dataset to 2000-2018
dealer <- dealer %>% 
  filter(Year >= 2000 & Year <= yr)

# Convert all gilled/gutted/trunked weights to whole fish weight (using conversion factor in dealer dataset)
dealer$conversion <- unlist(lapply(regmatches(dealer[,'FISH_CONDITION'],gregexpr("[[:digit:]]+\\.[[:digit:]]+\\*?$",dealer[,'FISH_CONDITION'])),"[",1))
dealer$conversionNum <- as.numeric(gsub(".$", "",dealer$conversion))
dealer[which(dealer$FISH_CONDITION==''),'conversionNum'] <- 1
dealer[which(dealer$conversionNum > 1), 'conversionNum'] <- 1
dealer[is.na(dealer$conversionNum),'conversionNum'] <- 1
dealer[which(dealer$conversionNum == 0),'conversionNum'] <- as.numeric(dealer[which(dealer$conversionNum == 0),'conversion'])
dealer$WholeWeight <- dealer$WeightGG/dealer$conversionNum
# For bigeye and yellowfin we are using better published conversion factors (non linear)
dealer[which(dealer$SPECIES_FK==3),'WholeWeight'] <- 1.189346*round(dealer[which(dealer$SPECIES_FK==3),'WeightGG'])^0.972009  # !!! This is for Yellowfin ONLY!!!
dealer[which(dealer$SPECIES_FK==6),'WholeWeight'] <- 1.274959*round(dealer[which(dealer$SPECIES_FK==6),'WeightGG'])^0.960613  # !!! This is for Bigeye ONLY!!!

# Weight - Length relationships for BET and swordfish
# Equation W=a*L^b
# Swordfish: a=1.299x10^(-5) and b=3.0738  Lf=143.6cm and Lm=102cm
SW_L50_kg_f <- (1.299*10^(-5))*(143.6^3.0738)
SW_L50_kg_m <- (1.299*10^(-5))*(102^3.0738)
# Convert from Length to weight using Phoebe's numbers for the Hawaii fishery (gives weight in grams)
# Bigeye: L=103cm from Farley et al. 2018 (WCPFC scientific committee report on age and growth of BET)
BET_L50_kg <- round(0.0187*(103^2.9605))/1000  # these numbers are from phoebe

# Merge dealer with logData so we can get departure dates
dealerAll <- dealer

# Unbin the dealer data for accurate length distributions
nTimes <- dealer$NUM_SOLD
dealerUnbin <- dealer[rep(seq_len(nrow(dealer)), nTimes), ]
dealerUnbin$NUM_SOLD <- 1
dealer <- dealerUnbin

# Parse out bigeye and swordfish datasets (bigeye is deep set only, swordfish is shallow set only)
dealerBE <- dealer[which(dealer$SPECIES_FK == 6 & dealer$SET_TYPE == 'D'),]
dealerSW <- dealer[which(dealer$SPECIES_FK == 11 & dealer$SET_TYPE == 'S'),]

# Functions
# ----------------------------
# Calculate the mean and sd for the violin plots
data_summary <- function(x) {
  m <- median(x)
  ymin <- as.numeric(quantile(x, 0.25))
  ymax <- as.numeric(quantile(x, 0.75))
  return(c(y=m,ymin=ymin,ymax=ymax))
}

dot_summary <- function(x) {
  xyy <- x %>% 
    group_by(Year) %>% 
    summarize(q20=quantile(WholeWeight, 0.20), q80=quantile(WholeWeight, 0.80)) %>% 
    tidyr::gather(Condition, Quartile, q20:q80)
  return(xyy)
}
# ----------------------------


# Read in images as rasters
beye <- readPNG("../../Bigeye_blue.png")
sword <- readPNG("../../Swordfish_yellow.png")
fish <- readPNG("../../Bigeye_white.png")
g <- rasterGrob(beye, interpolate=TRUE)
h <- rasterGrob(sword, interpolate=TRUE)
w <- rasterGrob(fish, interpolate=TRUE)

#--------------------------------------------------------
# Calculate annual box plot distributions for all three
#--------------------------------------------------------
# Function to make violin plot
violfunc <- function(x, col, img, L50, L50col, dotData, imgMin=2022, yMax=max(x$WholeWeight)) {
  numSold <- aggregate(NUM_SOLD ~ Year, data=x, FUN=sum)
  numSold$pretty <- prettyNum(c(numSold$NUM_SOLD), big.mark=",")
  #yMax <- max(x$WholeWeight)#200 # for the shortened all fish figure
  yMin <- round(yMax*1.1)
  int <- ifelse(length(seq(0,yMax,50)) > 7, 100, 50)
  dotData <- dot_summary(x)
  p <- ggplot(data=x, aes(x=Year, y=WholeWeight, group=Year)) + 
    annotation_custom(img, xmin=imgMin, xmax=yr+0.7, ymin=yMax*0.82, ymax=yMax*1.02) +
    geom_violin(fill=col, color=col, alpha=0.6) +
    stat_summary(fun.data=data_summary, size=0.2) +
    geom_point(data=dotData, aes(x=Year, y=Quartile), color='black', size=0.4) +
    geom_hline(yintercept=L50, linetype='dashed', color=L50col) +
    ylab('Weight (kg)') +
    theme_bw() + 
    theme(plot.margin=unit(c(1.2,1,1,1),'lines'), panel.grid=element_blank(), axis.text=element_text(size=9), 
          axis.title.x=element_blank(), axis.title.y=element_text(size=10)) +
    scale_y_continuous(breaks=seq(0,yMax,int), limits=c(0,yMax)) +
    scale_x_continuous(breaks=seq(2000,yr,1), expand = c(0.02,0.02), 
                       sec.axis=sec_axis(~., breaks=seq(2000,yr,1), labels=NULL)) 
  for (i in 1:length(numSold$Year)) {
    p <- p + annotation_custom(grob=textGrob(label=numSold$pretty[i], gp=gpar(cex=0.5, color='gray20')), ymin=yMin, ymax=yMin, xmin=numSold$Year[i], xmax=numSold$Year[i])
  }
  p <- p + annotation_custom(grob=textGrob(label='n=', gp=gpar(cex=0.5, fontface='italic', color='gray20')),ymin=yMin, ymax=yMin, xmin=1999, xmax=1999)
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  return(gt)
}
# Make three part violin plot
a <- violfunc(dealer, '#63a298', w, -1, NA)
a2 <- violfunc(dealer, '#63a298', w, -1, NA, yMax=200)
b <- violfunc(dealerBE, '#5b84a1', g, BET_L50_kg, 'gray50')
c <- violfunc(dealerSW, '#ffba37', h, c(SW_L50_kg_f, SW_L50_kg_m), 'gray50', 2020)
setwd(file.path(mainDir, 'Output'))
ggsave(plot=grid.arrange(a, a2, nrow=2), paste0('SAFEdealerViolinAllFish_', yr, '.png'), height=7, width=9, dpi = 300, units = 'in')
ggsave(plot=grid.arrange(b, c, nrow=2), paste0('SAFEdealerViolinBETSWD_', yr, '.png'), height=7, width=9, dpi = 300, units = 'in')

#--------------------------------------------------------
# Calculate Climatology of weight distribution for all fish, bigeye and swordfish 
#--------------------------------------------------------
# Bin all fish larger than 200kg
# All Fish
dealerSAFE <- dealer
dealer[which(dealer$WholeWeight >= 200), 'WholeWeight'] <- 200
dealerClim <- dealer[which(dealer$Year <=2009),]
dealerYr <- dealer[which(dealer$Year == yr),]
# Bigeye
dealerBE[which(dealerBE$WholeWeight >= 200), 'WholeWeight'] <- 200
dealerBEClim <- dealerBE[which(dealerBE$Year <= 2009),]
dealerBEYr <- dealerBE[which(dealerBE$Year == yr),]
# Swordfish
dealerSW[which(dealerSW$WholeWeight >= 200), 'WholeWeight'] <- 200
dealerSWClim <- dealerSW[which(dealerSW$Year <= 2009),]
dealerSWYr <- dealerSW[which(dealerSW$Year == yr),]

# Using function and plot to make density plots
plotfunc <- function(x, y, color, img) {
  nc <- prettyNum(c(sum(x$NUM_SOLD)),big.mark=',')
  n17 <- prettyNum(c(sum(y$NUM_SOLD)),big.mark=',')
  dclim <- hist(x$WholeWeight, breaks=ifelse(nc<1000000,20,40), plot=F)
  dyr <- hist(y$WholeWeight, breaks=ifelse(nc<1000000,20,40), plot=F)
  ymx <- max(c(dclim$counts/sum(dclim$counts), dyr$counts/sum(dyr$counts)))*1.05
  plot(dclim$mids, dclim$counts/sum(dclim$counts), xlim=c(0,212), ylim=c(0,ymx), type='l', las=1, xlab=NA, ylab=NA, main=NA, yaxs='i')
  polygon(c(-5, dclim$mids,205), c(0,dclim$counts/sum(dclim$counts),0), col='gray70', border='gray50')
  mtext(side=1, line=2.5, text='Weight (kg)', cex=0.7)
  mtext(side=2, line=3.5, text='Proportion of Fish', cex=0.7)
  par(new=T)
  plot(dyr$mids, dyr$counts/sum(dyr$counts), xlim=c(0,212), ylim=c(0,ymx), type='l', xlab=NA, ylab=NA, main=NA, axes=F, yaxs='i', col=color)
  polygon(c(-5, dyr$mids, 205), c(0,dyr$counts/sum(dyr$counts),0), col=color, border=color)
  text(160,ymx*0.92, pos=4, labels=bquote(paste('n'['climatology']*' = ', .(prettyNum(nc,big.mark=','), sep=' '))), cex=0.9)
  text(160, ymx*0.85, pos=4, labels=bquote(paste('n'[.(yr)]*' = ', .(prettyNum(n17,big.mark=','), sep=' '))), cex=0.9)
  grid.raster(beye, 0.85,0.6, width=.11)
  grid.raster(sword, 0.85,0.26, width=.14)
  
}
# Make three part Density Kernel plot
setwd(file.path(mainDir, 'Output'))
png(paste0('SAFEdealerWeightDensDist_', yr, '.png'), width=5, height=7, res=300, units='in')
par(mfrow=c(3,1))
par(mar=c(4,5,0.25,0.5))
plotfunc(dealerClim, dealerYr, '#63A298A0') #alpha=180 of 250 for the violin version in plot above
plotfunc(dealerBEClim, dealerBEYr, '#5B84A1A0') 
plotfunc(dealerSWClim, dealerSWYr, '#FFBA37A0')  
dev.off()
