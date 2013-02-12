# Functions for marine scotland tender

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# simulate data from a baseModel, nsim times
# INPUT:
    # baseModel:  model fitted object such as glm
    # nsim:       number of sets of simulated data to generate
    # seed:       set the seed so that data generation can be replicated
# OUTPUT:
    # simfits:  n x nsim matrix of values.  each column is a simulated dataset (from a poisson)
getSimData<- function(baseModel, nsim = 100, seed=1234){
  set.seed(seed)
  # simulate data. make simfits which is 100 x nrow(data)
  n<- length(baseModel$y)
  simfits<- matrix(NA, nsim, n)
  for(i in 1:n){
    simfits[,i]<- rpois(n=nsim, fitted(baseModel)[i])
  }
  # transpose simfits so each column is a dataset
  simfits<- t(simfits)
  return(simfits)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# function to check the true coeff lie within 95% of sim coefficients
# INPUT:
    # baseModel:  model fitted object such as glm
    # simfits:    n x nsim matrix of values.  each column is a simulated dataset (from a poisson)
# OUTPUT:
    # simcoeff: nsim x number of coeff matrix of coefficients for each of the simulated datasets in simfits.
checkSimCoeff<- function(baseModel, simfits){
  # fit models to each simfits data
  simcoeff<- matrix(NA, ncol(simfits), length(baseModel$coeff))
  for(i in 1:ncol(simfits)){
    if((i/10)%%1 == 0 ){cat(i, '\n')}else{cat('.')}
    resp<<- simfits[,i] # need to have new variable in workspace
    fit<- update(baseModel, resp~.)
    # store coeffs
    simcoeff[i,]<- fit$coefficients
  }
  rm(resp, pos=.GlobalEnv) # removes randomly made resp
  return(simcoeff)
} # end func

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Function to plot histogram of simulated coefficients with main blue line as the true coefficient from the baseModel and dashed red lines are the 95% quantiles of the simulated coefficients with mean (red line)
# data:   vector of beta coefficients
# truth:  true beta coefficient from baseModel
# title:  title for the plot
plotCoeffIntervals<-function(data, truth, title){
  hist(data, xlab='', main=title)
  abline(v=truth, col='blue', lwd=3)
  abline(v=mean(data), col='red', lwd=3)
  abline(v=quantile(data, probs=c(0.025, 0.975)), col='red', lwd=3, lty=2)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~