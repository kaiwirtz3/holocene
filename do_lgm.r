
# This script runs logistic regression models for different sets of input variables,
# predicting the probability of a boom or bust (binary outcomes) based on a set of selected variables.
# It computes regression coefficients and stores results for further analysis.

# SPDX-FileCopyrightText: 2024 2023-2024 Helmholtz-Zentrum hereon GmbH
# SPDX-FileContributor: Kai W. Wirtz <kai.wirtz@hereon.de>
# SPDX-License-Identifier: GPL-3.0-or-later

# Input files:
# - Data input variables are expected to be loaded or processed earlier in the script

# Output files:
# - None explicitly produced, but regression coefficients and probabilities are calculated and stored in variables.

# Variables used:
# - vi: selected input variable indices
# - nvi: number of input variables selected
# - signbb: sign array for boom and bust distinction
# - prob: array storing predicted probabilities
# - coresb: array storing regression coefficients
# - rgr0: original RGR (Relative Growth Rate) time series
# - stdc: standard deviation correction for RGR
# - vart: selected variable matrix used in the logistic regression
# - mydata: data frame containing binary response and selected variables
# - mylogit: logistic regression model
# - predicted: predicted probabilities from the logistic regression
# - statres: model summary containing coefficients

nvi <- length(vi)

if (nvi == 1) {
  vi <- c(vi, vi)
}

signbb <- c(1, -1)
prob <- array(rep(0, 2 * length(ii)), c(2, length(ii)))
coresb <- array(rep(0, (1 + nvi) * 4), c(2, 1 + nvi, 4))

# Loop over booms and busts
for (bb in 1:2) {

  # Create (binary) response variable input matrix
  rgr1 <- ifelse(rgr0 * signbb[bb] > stdc, 1, 0)
  vart <- var[, vi]
  colnames(vart) <- shortname[1:max(2, nvi)]
  mydata <- data.frame(rgr = rgr1, vart)

  # Create Logistic Regression Model
  switch(nvi,
    mylogit <- glm(rgr ~ V1, data = mydata, family = "binomial"),
    mylogit <- glm(rgr ~ V1 + V2, data = mydata, family = "binomial"),
    mylogit <- glm(rgr ~ V1 + V2 + V3, data = mydata, family = "binomial")
  )

  # Predict probability
  predicted <- predict(mylogit, mydata, type = "response")
  prob[bb, ] <- predicted

  statres <- summary(mylogit)
  cores <- statres$coefficients

  # Store coefficients
  coresb[bb, , ] <- cores

  ## CIs using profiled log-likelihood
  ## confidence intervals for the coefficient estimates
  # confint(mylogit)
  # plot(predicted, rgr, 'p', col = "blue")

  # First skill measure
  # i2 <- which(rgr0[ii] * signbb[bb] > 0.5 * sd(rgr0[ii]))
  # Pearson <- cor(predicted[i2], rgr0[ii[i2]])

  # prediscr <- ifelse(predicted > 0.5, 1, 0)
  # isid <- abs(prediscr - rgr1)
} # for bb
