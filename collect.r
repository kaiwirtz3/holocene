# This script processes grid cell data to prepare clustering using an extended k-means algorithm.
# It merges geographical and statistical information from multiple grid patches, calculating
# growth rate classifications based on p-values and q-values, and generates growth matrices.
# Finally, it saves the processed data in a binary format.

# SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
# SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
# SPDX-License-Identifier: GPL-3.0-or-later

# Input files:
# - 'out/bin/eurogrid_*.bin': Contains grid cell data for multiple patches.

# Output files:
# - 'out/bin/eurodat_all.bin': Processed data including locations, growth matrices, and clustering information.

# Variables:
# - scdir: Directory for input and output files.
# - ri0, ri1: Range of grid patches to process.
# - plim: Statistical significance threshold for p-values and q-values.
# - lons, lats: Longitude and latitude of sites.
# - pps, rmm: Growth matrices for each patch.
# - pcl: Cluster data for each patch.
# - sit: Site information for each patch.
# - s0: Counter for site indices.

library(sf)

scdir <- 'out/'  # main data IO directory

# range of grid patches
ri0 <- 1
ri1 <- 64

# statistical rigor
plim  <- 0.05

# clean variables
lons <- c()
lats <- lons
pps <- lons
rmm <- lons
pcl <- lons
sit <- lons
s0 <- 0

# loop over grid patches
for (ri in ri0:ri1) {

  file <- paste0(scdir, 'bin/eurogrid_', ri, '.bin')
  print(file)

  if (file.exists(file)) {

    load(file)
    print(paste(ri, length(pcl), length(pclust), 'pclust=', pclust[1], pclust[length(pclust)]))

    # checks if patch contains sites/dates
    if (length(pclust) > 0) {

      pcl <- c(pcl, pclust)  # merges patch data

      locs <- st_coordinates(eurospatial$locations)  # geolocation data

      # position of sites
      lon <- as.vector(locs[, 1])
      lat <- as.vector(locs[, 2])

      lons <- c(lons, lon)  # merges sites' georefs
      lats <- c(lats, lat)
      sit <- c(sit, s0 + sites)
      s0  <- s0 + sites[length(sites)]

      nn <- length(eurospatial$pvalHi[, 1])
      dt <- length(eurospatial$pvalHi[1, ])

      # create growth matrix
      pp <- matrix(nrow = nn, ncol = dt)
      pp0 <- pp
      j <- 1

      # loop over time segments
      for (ti in seq(1, dt)) {

        pvalHi <- eurospatial$pvalHi[, ti]
        pvalLo <- eurospatial$pvalLo[, ti]
        qvalHi <- eurospatial$qvalHi[, ti]
        qvalLo <- eurospatial$qvalLo[, ti]

        # classify according to p-value
        iHi <- which(pvalHi <= plim)
        iLo <- which(pvalLo <= plim)

        # classify also according to q-value
        iHi2 <- which(qvalHi[iHi] <= plim)
        iLo2 <- which(qvalLo[iLo] <= plim)

        # mean rate of change
        rt <- 0.5 * asinh(rgr[ti] / 2.E-3)
        if (is.na(rt) | is.nan(rt)) { rt <- 0 }
        fac <- 1
        p <- rep(rt, nn)
        p0 <- p
        p[iHi] <- p[iHi] + 1 * fac
        p[iHi[iHi2]] <- p[iHi[iHi2]] + 2 * fac
        p[iLo] <- p[iLo] - 1 * fac
        p[iLo[iLo2]] <- p[iLo[iLo2]] - 2 * fac
        pp[, j] <- p  # discrete growth levels into matrix
        pp0[, j] <- p0  # discrete growth levels into matrix

        j <- j + 1
      }  # end for ti

      # append growth matrices
      pps <- rbind(pps, pp)
      rmm <- rbind(rmm, pp0)
    }
  } else {
    print(paste(file, 'not found!'))
  }
}  # end for ri

# write binary data
datc <- data.frame(lons, lats, pps)
file <- paste0(scdir, 'bin/eurodat_all', '.bin')
print(paste('data ready to ', file))
save(file = file, datc, rmm, pcl, breaks, sit)
