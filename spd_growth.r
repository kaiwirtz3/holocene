# This script calculates Summed Probability Density (SPD) and related growth (RGR) for each region and time slice using an array of methods.
# It processes C14 dates, calibrates the dates, creates time bins, calculates SPD, computes relative change rates (RGR),
# smooths the data, and generates plots for visual analysis. The script also saves the processed data and plots to files.

# SPDX-FileCopyrightText: 2024 2023-2024 Helmholtz-Zentrum hereon GmbH
# SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
# SPDX-License-Identifier: GPL-3.0-or-later

# Input files:
# - c14mat/C14_europe0.mat
# - out/bin/clusti_<cc[1]><ti>_120.bin

# Output files:
# - out/plots/spd_clu_<tag>_<ti>_<nplot>.png
# - out/mat/PrePop_<tag>_<tii>_<i>.mat
# - out/mat/AllPop_<tag>_<tii>.mat

# Variables:
# - cc: Continent name
# - breaks: Time windows with given region outlay
# - binv: Binning size in SPD calculation
# - ntag: Number of SPD methods
# - nbrk: Number of time breaks
# - myargs: Command line arguments for method and time slice
# - an: Argument number
# - tagi0, tagi1: Start and end indices for SPD methods
# - ti0, ti1: Start and end indices for time slices
# - scdir: Directory for input/output files
# - tagv: List of SPD methods
# - nmaxclu: Maximum number of clusters
# - dtl: Time window length
# - dt: Time step
# - norms: Normalization options
# - dat: Data read from input file
# - ncol, nrow: Number of columns and rows for figure layout
# - colv: Color vector for plots
# - nplot: Plot counter
# - ymsa: Matrix to store smoothed density
# - ti: Current time slice
# - file: File name for region/cluster info
# - clustn: Cluster numbers
# - nregions: Number of regions
# - ymv, rgrv: Matrices to store density and RGR
# - n0, ndates: Vectors to store number of sites and dates
# - lat, lon: Latitude and longitude of sites
# - ii: Indices of sites in current region
# - lati, loni: Latitude and longitude of sites in current region
# - ci, sites: Vectors to store site indices
# - i2: Indices of valid C14 dates
# - X: String to print number of valid dates
# - eurodates: Calibrated C14 dates
# - eurobins: Time bins for SPD calculation
# - timeRange: Time range for SPD calculation
# - steps: Time steps for SPD calculation
# - tirgr: Time vector for RGR calculation
# - clu_spd: SPD calculation result
# - clu_rgra: Relative change rate (RGR) calculation result
# - rgr: Relative change rates (RGR)
# - y, ym: Density and smoothed density
# - mi: Indices for storing data
# - tm, tm2: Time vectors for storing data
# - fac: Scaling factor for RGR plot

rm(list = ls())
library(rcarbon)  # RCARBON by Crema2017
library('R.matlab')
source("movavg.r")

# settings
cc = 'Europe'  # continent name
breaks = seq(3000, 9800, 400) # time windows with given region outlay
binv = 100    # only non-normalized
ntag = length(binv)
nbrk = length(breaks)

# method and time slice as index from input argument
myargs = commandArgs(trailingOnly = TRUE)
print(paste('args:', length(myargs), ' nbrk=', nbrk))

if (length(myargs) > 0) {
  an = as.numeric(myargs[1])
  tagi0 = 1 + floor(an / nbrk)
  tagi1 = tagi0
  ti0 = 1 + (an %% nbrk)
  ti1 = ti0
} else {
  ti0 = 1
  ti1 = nbrk
  tagi0 = 1
  tagi1 = ntag
}

# input/output directory
scdir = 'out/'
# list of SPD methods
tagv = c('_NoNorm_Bin100')
nmaxclu = 45
dtl = 325
dt = 25  # time-step

# read matlab C14 data
norms = c('', 'No')
dat <- readMat(paste0('c14mat/C14_europe0.mat'))

# figure settings
ncol = 2
nrow = 3 # figure outlay
colv = rainbow(5)

# loop over SPD methods
for (tagi in tagi0:tagi1) {
  nplot = 0 # reset of variables
  tag = paste0('_', norms[1 + (tagi %% 2)], 'Norm_Bin', binv[tagi])
  print(paste(tagi, tag, ti0, ti1))
  ymsa = NULL

  # loop over time slices
  for (tii in seq(ti0, ti1)) {
    ti = breaks[tii]

    # read region/cluster info
    file <- paste0(scdir, 'bin/clusti_', cc[1], ti, '_120.bin')
    load(file)    # variables: clusti, k, wi, cluc
    clustn = seq(k)
    nregions = k    # number of regions
    nplot = 0

    # prepares fields
    ymv = array(NaN, c(nregions, ceiling((2 * dtl) / dt)))
    rgrv = array(NaN, c(nregions, ceiling((2 * dtl) / dt) - 1))
    n0 = array(0, nregions)
    ndates = n0

    # locations of regionalized C14 dates
    lat = round(dat$lat, digits = 3)
    lon = round(dat$lat, digits = 3)

    # loop over regions
    for (i in 1:nregions) {
      if (i %% (ncol * nrow) == 1) {
        if (nplot > 0) {
          dev.off()
        }
        print(paste('new figure', tii, i, nplot))
        file = paste0(scdir, "plots/spd_clu", tag, "_", ti, "_", nplot, ".png")
        png(file, width = 980, height = 940, units = 'px')
        par(oma = c(1, 0, 1, 1), mar = c(0.1, 0.1, 0.14, 0.5), mfrow = c(nrow, ncol), cex.lab = 1.5, cex.sub = 0.5, cex.main = 1., cex.axis = 2)
        nplot = nplot + 1
      }

      # pointer to region/cluster and positions
      ii <- which(clusti == i)
      lati = round(dat$lat[ii], digits = 3)
      loni = round(dat$lon[ii], digits = 3)
      n0[i] = length(ii)

      # clean vectors
      ci = c()
      sites = ci

      # loop over dates to fill site vector
      for (ij in 1:n0[i]) {
        ni = which(lat == lati[ij] & lon == lati[ij])
        ci = c(ci, ni)
        sites = c(sites, rep(ij, length(ni)))
      }

      # filter NaNs
      i2 <- which(!is.na(dat$C14agesn[ci]))
      ii = ci[i2]
      sites = sites[i2]
      ndates[i] = length(i2)

      # print number of valid C14 dates per region
      X <- sprintf('%2d %d sites valid dates: %d ', i, n0[i], length(ii))
      print(X)

      # calibration using intcal20
      eurodates <- calibrate(dat$C14agesn[ii], dat$C14SDsn[ii], calCurves = 'intcal20', ncores = 4, normalised = (tagi %% 2 == 0))
      print(paste('eurodates with normalised=', (tagi %% 2 == 0), ' ready ...'))

      # create time bins
      eurobins <- binPrep(sites = sites, ages = dat$C14agesn[ii], h = binv[tagi])

      # set time vector and calculate SPD
      timeRange <- ti + c(dtl, -dtl)
      steps <- seq(timeRange[1], timeRange[2], -dt)
      tirgr <- steps[2:(length(steps) - 1)]
      clu_spd <- spd(x = eurodates, bins = eurobins, timeRange = timeRange)

      # calculate relative change rate (RGR)
      clu_rgra <- spd2rc(clu_spd, breaks = steps)
      rgr <- clu_rgra$roca

      # retrieve and smooth density
      y <- clu_spd$grid$PrDens
      ym <- movavg(y, dt)

      # store into matrix
      mi <- seq(round(0.5 * dt), floor(length(clu_spd$grid$calBP) / dt) * dt, dt)
      ym <- ym[mi]
      ymsa <- cbind(ymsa, ym)
      tm <- clu_spd$grid$calBP
      tm2 <- tm[mi]

      ymv[i, 1:length(ym)] <- ym
      rgrv[i, 1:length(rgr)] <- rgr

      # plot SPD with own method
      plot(tm, y, col = "indianred", lwd = 1, lty = 1)

      # plot smoothed SPD, also of other methods
      lines(tm2, ymv[i, ], col = colv[1], lwd = 2, lty = 1)
      lines(tirgr, mean(ymv[i, ]) + 5 * rgrv[i, ], col = colv[1], lwd = 3, lty = 2)

      text(ti - dtl * 0.5, mean(y) * 0.15, labels = paste(i, ti), cex = 3, col = NULL)
      text(ti, mean(y) * 0.05, labels = length(ii), cex = 2, col = NULL)

      if ((i) %% (ncol * nrow) == -1) {
        legend(x = "topleft", legend = tagv)
      }
      if (i %% 8 == 0 | i == nregions) {
        file = paste0(scdir, "mat/PrePop", tag, '_', tii, '_', i, ".mat")
        print(paste0("write data to ", file))
        writeMat(file, poptime = tm, ymv = ymv, trgr = tirgr, rgr = rgrv, nreg = nregions, nsites = n0, ndates = ndates)
      }
    }

    # save population data as Matlab binary
    file = paste0(scdir, "mat/AllPop", tag, '_', tii, ".mat")
    print(paste0("write data to ", file))
    writeMat(file, poptime = tm2, ymv = ymv, trgr = tirgr, rgr = rgrv, nreg = nregions, nsites = n0, ndates = ndates)

    print(paste('close figure', tii, i, nplot))
    dev.off()

    tii = tii + 1
  }
}
