# This script calculates spatial statistics of SPDs on a 4°x4° grid.
# It checks for input arguments that control partial computation in a parallel mode.
# Examples: "Rscript grid_growth.r $SLURM_ARRAY_TASK_ID" within a slurm.sh
#           "Rscript grid_growth.r 1 3" to process only from patch 1 to 3
#           "for (( i=1; i<9; i++ )); do (Rscript grid_growth.r $(expr $i \* 8 - 7) $(expr $i \* 8) &) ; done"

# SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
# SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
# SPDX-License-Identifier: GPL-3.0-or-later

# Input files:
# - 'data/seamask_norm_0.05.mat': Matlab seamask for entire domain
# - 'c14mat/C14_europe_<ri>.mat': Matlab C14 data for single grid patch

# Output files:
# - 'out/plots/neo_Gr_<ri>.png': Graphics output for each patch
# - 'out/bin/eurogrid_<ri>.bin': Binary data for each patch

# Variables:
# - ri0, ri1: Range of patches to process
# - eurospatial: Spatial permutation test results
# - pclust, pclust0: Indices to dates
# - ktot, nclustsum: Counters for total and sum of clusters
# - v: Seamask data
# - dat: C14 data for a single grid patch
# - loi, lai: Indices for longitude and latitude boundaries
# - ng, nl: Number of grid points and land points
# - rland: Ratio of land to sea
# - SDc: Standard deviation cutoff for filtering dates
# - ii: Indices of filtered dates
# - edates: Calibrated dates
# - lat, lon: Latitude and longitude of dates
# - sites0, usite, usitei, usitei0: Site information
# - sites: Vector with site index
# - ebins: Binned data
# - sites2, sites1: Data frames of site locations
# - timeRange, breaks, nb: Time window settings
# - eurospd: SPD for site distribution
# - rgra, rgr: RGR for time sequence
# - base: Basemap for plotting
# - xrange, yrange: Bounding coordinates of site distribution
# - pfile: File path for graphics output

library(rworldmap)
library(rcarbon)
library('R.matlab')
library("RColorBrewer")
library(sf)

scdir='out/'

myargs = commandArgs(trailingOnly=TRUE)
if (length(myargs)>0) {
  ri0 = as.numeric(myargs[1])
  if (length(myargs)>1) {
    ri1 = as.numeric(myargs[2])
  } else {
    ri1 <- ri0
  }
  if (ri1 < ri0) {
    ri1  <-  ri0
  }
} else {
  ri0 <- 1
  ri1 <- 64
}

# clean variables
print(c(ri0,ri1))
eurospatial <- c()
pclust <- c()
ktot <- 0
nclustsum <- 0

# read matlab seamask for entire domain
file <- 'data/seamask_norm_0.05.mat'
v <- readMat(file)
dim(v$value)

# loop over patch / data segments (entire set exceeds memory)
for (ri in ri0:ri1) {

  # read matlab C14 data for single grid patch
  fname <- paste0('c14mat/C14_europe_',ri,'.mat')
  print(paste('read C14 data from',fname))
  dat <- readMat(fname)
  n0  <- length(dat$C14age) # number of dates

  # retrieves boundary long-lat for patch
  loi <- which(v$long>=dat$bo[1] & v$long<dat$bo[3])
  lai <- which(v$latg>=dat$bo[2] & v$latg<dat$bo[4])

  # estimates fraction of sea/land
  ng <- length(v$value[loi,lai])
  nl <- length(which(v$value[loi,lai]>=0))
  rland  <-  ng/(nl+1)
  rland  <-  min(10,rland)
  print(paste('sea grid',ng,nl,round(rland,digits=3),'lo',loi[1],loi[length(loi)],'la',lai[1],lai[length(lai)]))

  # filters for high quality=low SD C14-dates depending on density: total number and land ratio
  SDc  <- 220 - sqrt(n0*rland)*2.5
  SDc  <- max(SDc,40)
  ii   <- which(dat$C14SD < SDc)

  n <- length(ii) # number of filtered dates
  if (n>0) {

    # calibration using intcal20
    edates <- calibrate(dat$C14age[ii],dat$C14SD[ii],normalised=FALSE,verbose=FALSE,calCurves='intcal20')
    print(paste('dates done ... ',n,'from',n0))

    # check for same site position
    lat <- round(dat$lat[ii],digits=3)
    lon <- round(dat$lon[ii],digits=3)
    sites0 <- data.frame(lat=lat,lon=lon)
    usite  <- unique(sites0[,c('lat','lon')])
    usitei <- as.numeric(rownames(usite))
    usitei0<- usitei
    usitei <- c(usitei,length(ii)+1)

    # creates vector with site index
    sites  <- NULL
    for (ui in seq(length(usitei)-1)) {
      sites=c(sites,rep(ui,usitei[ui+1]-usitei[ui]))
    }
    print(paste('sites',length(sites),sites[1],sites[length(sites)],ui))

    # bins the data
    ebins <- binPrep(sites=sites,ages=dat$C14age[ii],h=100)
    print('ebins ready ...')

    sites2  <- data.frame(id=sites[usitei0],lat=lat[usitei0],lon=lon[usitei0])
    sites1  <- st_as_sf(sites2,coords=c('lon','lat'),crs=4326)

    # vector of indices to dates
    pclust  <- dat$SiteID[ii[usitei0]]
    pclust0 <- dat$SiteID[ii]

    # Create a data.frame of site locations extracting spatial coordinates
    rownames(sites2) <- sites2$id
    sites2 <- sites2[,-1]

    # Convert to a SpatialPoints class object:
    coordinates(sites2) <- c("lon","lat")
    proj4string(sites2) <- CRS("+proj=longlat +datum=WGS84")

    ## outcomment for version(rcarbon) < 1.5
    # Compute distance matrix
    ##d <- spDists(sites2,sites2,longlat=TRUE)
    # Compute spatial weights
    ##w <- spweights(d,h=100)
    ##print('d and w ready ...')

    # define time windows
    dt <- 400
    timeRange <- c(10200,2600)
    breaks    <- seq(timeRange[1],timeRange[2],-dt)
    nb <- length(breaks)-2

    # calc SPD for prepared site distribution
    eurospd = spd(x = edates,bins=ebins,timeRange = timeRange)

    # calc RGR for given time seqence
    rgra <- spd2rc(eurospd,breaks = breaks)
    rgr  <- rgra$roca

    # spatial permutation tests of sample sites to detect local deviations in rates of change in the SPD
    ## outcomment for version(rcarbon) < 1.5
    ##eurospatial <- sptest(calDates=edates, bins=ebins,timeRange=timeRange,locations=sites2,permute="locations",nsim=1000, breaks=breaks,spatialweights=w, ncores=4) #

    # sptest after rcarbon-1.5 based on sf
    eurospatial <- sptest(calDates=edates, bins=ebins,timeRange=timeRange,locations=sites1,locations.id.col='id',h=100,kernel='gaussian',permute="locations",nsim=100,breaks=breaks,ncores=4,verbose=FALSE)
    print(dim(eurospatial$qval))

    # output to a map
    base <- getMap(resolution="low")

    # extract bounding coordinates of the site distribution
    xrange <- bbox(sites2)[1,]
    yrange <- bbox(sites2)[2,]

    # prepare graphics output
    pfile <- paste0(scdir,'plots/neo_Gr_', ri,'.png')
    print(paste("plot to ",pfile))
    png(filename=pfile,width = 12, height = 12, units = "cm", res=300 )

    # tiles of subplots for all time segments
    nx <- round(sqrt(nb))
    par(oma = c(1, 1, 2, 1), mar = c(0.08, 0.08, 0.9, 0.3),mfrow=c(nx,ceiling(nb/nx)), cex.lab=0.5,cex.sub=0.5, cex.main=1.,cex.axis=0.7)

    ncolor <- 2
    cmap <- palette(rainbow(ncolor))
    pr   <-  sample(ncolor)

    # loop over time segments
    for (ti in seq(nb)) {
      transtis<- as.character(breaks[1+ti])
      file    <- paste0(scdir,"plots/EuroGr_",ri,'_',transtis,".pdf")

      # background map
      plot(base,col=rgb(0.95, 0.95, 0.95),border="antiquewhite3",xlim=xrange,ylim=yrange,main=paste(transtis,"kBP ",round(rgr[ti]*1E3)))

      plot(eurospatial,index=ti,xlim=xrange,ylim=yrange,option="test",add=TRUE,legend=(ti==1),legSize=0.7,location="left")
    }
    dev.off()
  }
  print(length(pclust))

  # write binary data
  file <- paste0(scdir,'bin/eurogrid_', ri,'.bin')
  print(paste('data ready to ',file))
  save(file=file,eurospatial,pclust,pclust0,usitei0,sites,breaks,rgr)
}

# run postprocessing only if this script has been executed for all patches
# otherwise "cluster_growth.r" has to be called manually
# writeMat(paste0(scdir,'clusti.mat'),clusti=clusti)
