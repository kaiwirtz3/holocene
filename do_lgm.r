    nvi=length(vi)
    if (nvi==1) {
       vi <- c(vi,vi) }
    signbb <- c(1,-1)
    prob <- array(rep(0,2*length(ii)), c(2,length(ii)))
    coresb <- array(rep(0,(1+nvi)*4), c(2,1+nvi,4))
    # loop over booms and busts
    for (bb in 1:2) {
      # create (binary) response var input matrix
      rgr1 <- ifelse(rgr0*signbb[bb]> stdc, 1, 0)
      vart <- var[,vi]
      colnames(vart) <- shortname[1:max(2,nvi)]
      mydata <- data.frame(rgr=rgr1, vart)
      # create Logistic Regression Model
      switch(nvi,
        mylogit <- glm(rgr ~ V1, data = mydata,family = "binomial"),
        mylogit <- glm(rgr ~ V1+V2, data = mydata,family = "binomial"),
        mylogit <- glm(rgr ~ V1+V2+V3, data = mydata,family = "binomial"))
      #
      #predict probability of defaulting
      predicted <- predict(mylogit, mydata, type="response")
      prob[bb,] <- predicted
      statres <- summary(mylogit)
      cores <- (statres$coefficients)
      # For every one unit change in V1, the log odds of rgr (versus non-rgr) increases by ..
      coresb[bb,,] <- cores
      ## CIs using profiled log-likelihood
      ## confidence intervals for the coefficient estimates
      #confint(mylogit)
      # plot(predicted,rgr,'p',col="blue")

      # first skill measure
      #i2=which(rgr0[ii]*signbb[bb]>0.5*sd(rgr0[ii]))
      #Pearson=cor(predicted[i2],rgr0[ii[i2]])

      # prediscr  <- ifelse(predicted > 0.5, 1, 0)
      # isid= abs(prediscr - rgr1)
    } # for bb
