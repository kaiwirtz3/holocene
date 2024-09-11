# This script runs a Generalized Linear Model (GLM) for different sets of input variables,
# calculates autoregressive lag corrections for RGR (Relative Growth Rate) based on Holocene periods,
# and generates probability outputs for various variable combinations. The script produces plots comparing RGR
# and the difference between boom and bust probabilities. Model results and coefficients are saved in MATLAB files.

# SPDX-FileCopyrightText: 2023-2024 Helmholtz-Zentrum hereon GmbH
# SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
# SPDX-License-Identifier: GPL-3.0-or-later

# Input files:
# - 'target_ts_0.mat' (contains time series data and input variables)

# Output files:
# - GLM results: 'out/glmres_<variable_combination>.mat'
# - Plot comparing RGR with probability: 'out/plots/comp_GLM_3V.png'

# Variables used:
# - scdir: directory for output files
# - tmax: maximum time for analysis
# - irgr: index for RGR (Relative Growth Rate) in the input data
# - vari: indices for the input variables (e.g., TSI, climate stability)
# - nvar: number of input variables
# - rgr0: original RGR time series
# - lag: time lag for autoregressive adjustments
# - time: time vector for analysis
# - var: input variable matrix after squeezing time domain
# - stdc: standard deviation correction for RGR
# - varname, varshort, shortname: arrays storing full, short, and variable names
# - nc: number of variable combinations
# - m: matrix defining combinations of input variables
# - prob: probability matrix
# - cores: model coefficients
# - tprob: difference between boom and bust probabilities

library('R.matlab')

scdir <- 'out/'
tmax <- 8.8  # maximum time

# indices of variables in input matrix
irgr <- 2  # 4:area weighted +smoothed +detrended
# index to input variables
vari <- c(6, 9, irgr)  # 6:- TSI 9: climate stability
nvar <- length(vari)
bobu <- c('booms', 'busts', 'bobu', 'single')
corn <- c('hit.bin', 'Pearson', 'sum', 'single')

# load data
fname <- paste0(scdir, 'target_ts_0', '.mat')
print(paste('reading...', fname))
d <- readMat(fname)

time <- d$dat[,1]
dt <- time[2] - time[1]

rgr0 <- d$dat[, irgr]  # RGR

# autoregressive part by lagged RGR

# lag for early Holocene
lag <- floor(0.36 / dt)  # 360a for t>7kaBP
ii <- which(time >= 7 & time < max(time) - 0.35)
d$dat[ii, irgr] <- rgr0[lag + ii]

# lag for early Holocene
ii <- which(time < 7)
lag <- floor(0.21 / dt)  # 210a for t>7kaBP
d$dat[ii, irgr] <- rgr0[lag + ii]
d$dat[is.na(d$dat[, irgr]), irgr] <- 0

# squeeze on time domain of interest
ii <- which(time >= 3 & time <= tmax)
time <- time[ii]
var <- d$dat[ii, vari]
rgr0 <- rgr0[ii]

stdc <- 0.5 * sd(rgr0)
var[is.nan(var)] <- 0
var[is.na(var)] <- 0

datnames <- paste(unlist(d$legdat))
gsub("[\r\n]", "", datnames)  # clear names from special chars

# create array of variable names in different versions
varname <- c()
varshort <- c()
shortname <- c()

for (cti in 1:nvar) {
  fullname <- datnames[vari[cti]]
  varname <- c(varname, fullname)
  fullname <- gsub(" ", "", fullname)
  str <- substring(fullname, 1, 4)
  varshort <- c(varshort, str)
  shortname <- c(shortname, paste0('V', as.character(cti)))
  print(varname[cti])
}

# create all possible combinations of input variables
nc <- 3  # number of input combinations
m <- matrix(c(c(0, 0, 1), c(1, 1, 0), c(1, 1, 1)), ncol = nc, nrow = 3)

# loop over variable combinations
for (loi in 1:nc) {
  vi <- which(m[, loi] == 1)
  source('do_lgm.r')
  print(cores)
  tprob <- prob[1, ] - prob[2, ]

  # save model coefficients and results to matlab file
  fname <- paste0(scdir, 'glmres_', paste0(varshort[vi], collapse = '_'), '.mat')
  writeMat(fname, timres = time, prob = prob, tprob = tprob, statcoeff = coresb, vind = vi)
}

# plot comparison RGR with probability
file <- paste0(scdir, 'plots/comp_GLM_3V', '.png')
png(file)
plot(time, rgr0 / sd(rgr0), 'l', col = "blue", main = paste('RGR vs. prob(boom)-prob(bust)'))
lines(time, tprob / sd(tprob), 'l', col = "red")
dev.off()
