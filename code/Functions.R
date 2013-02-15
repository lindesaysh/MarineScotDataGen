
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~ Simulate, refit, cis, preds ~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wrapperfunc<-function(model, grid, seed){
  set.seed(i)
  # generate one set of simulated data
  simData<-getSimData(fitted(model))
  
  # fit model to simData
  newfit<-updateResp(model, simData)
  
  # make confidence interval on betas
  betaCIs<-makeCIs(coef=newfit$coefficients, df=newfit$df.residual, se=summary(newfit)$coef[,2])
  
  # check for true coef inclusion in CIs
  betaCheck<- checkBetaCIs(model$coefficients, betaCIs)
  
  # make prediction to a grid
  preds<- makePredictions(newfit, grid)
  
  return(list(simData=simData, fit=newfit, betaCIs=betaCIs, betaCheck = betaCheck, preds=preds))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~ getSimData - poisson sim ~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate data from a baseModel, nsim times
# INPUT:
# baseModel:  model fitted object such as glm
# OUTPUT:
# simfits:  n x nsim matrix of values.  each column is a simulated dataset (from a poisson)
getSimData<- function(data){
  # simulate data. make simfits which is 100 x nrow(data)
  n<- length(data)
  simdata<-rpois(n, data)
  return(simdata)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~ updateResp - own resp update func ~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# updating a model with a new response vector.  
# prevents the need for update() which needs certain parameters to exist in the workspace rather than the function space.
updateResp<-function(model, newresp){
  teval<-paste('fit<-glm(newresp ~', paste(attr(model$terms, "term.labels"), sep='', collapse='+'), ', family=poisson, data=model$data)')
  eval(parse(text=teval))
  return(fit)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~ makeCIs - get cis for each coeff ~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# uses usual ci formula to get 95% cis for each parameter in model
# INPUT: 
# coef: vector of model coefficients (length b)
# df:   degrees of freedom
# se:   standard error of model coefficients
# OUTPUT:
# cis:  ncoeffx2 matrix of upper and lower 95% cis for each coeff

makeCIs<- function(coef, df, se){
  cis<- matrix(NA, length(coef), 2)
  for(i in 1:nrow(cis)){
    cis[i,1]<-coef[i] + qt(0.025, df)*se[i]
    cis[i,2]<-coef[i] - qt(0.025, df)*se[i]
  }
  return(cis)
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~ sim beta CI true beta check ~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# work out if interval contain truth
# INPUT:
  # vector of true coefficients, length b
  # matrix of cis for a set of simulated coeffs (bxn)
# OUTPUT:
  # vector of 0/1 (1 for truth within CI) of length b
checkBetaCIs<-function(truecoef, cis){
  marker<- rep(0, nrow(cis))
  for(i in 1:nrow(cis)){
    marker[i]<-ifelse(truecoef[i]< cis[i,2] & truecoef[i]>cis[i,1], 1, 0)
  }
  return(marker)
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~ make predictions for a model ~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
makePredictions<- function(model, grid){
  #preds <- (coef[-1] %*% grid) + coef[1]
  preds<- predict(model, grid)
  return(preds)
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~ plotCoeffIntervals - beta check ~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to plot histogram of simulated coefficients with main blue line as the true coefficient from the baseModel and dashed red lines are the 95% quantiles of the simulated coefficients with mean (red line)
# data:   vector of beta coefficients
# truth:  true beta coefficient from baseModel
# title:  title for the plot
plotCoeffIntervals<-function(data, truth, title){
  hist(data, xlab='', main=title)
  abline(v=truth[1], col='blue', lwd=3)
  abline(v=mean(data), col='red', lwd=3)
  abline(v=quantile(data, probs=c(0.025, 0.975)), col='red', lwd=3, lty=2)
  abline(v=truth[2:3], col='blue', lwd=3, lty=2)
}