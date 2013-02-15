# Functions for marine scotland tender

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~ getSimData - poisson sim ~~~~~~~~~~~~~~~
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
# ~~~~~~~~~~~ checkSimCoeff - using update ~~~~~~~~~~~
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
# ~~~~~~~~~ checkSimCoeff2 - own update func ~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to check the true coeff lie within 95% of sim coefficients
# same as checkSimCoeff func above but uses my own function for updating the glm formulae so that the new responses dont have to be in the workspace. 
# NB ~ 4% slower than func above
# INPUT:
# baseModel:  model fitted object such as glm
# simfits:    n x nsim matrix of values.  each column is a simulated dataset (from a poisson)
# OUTPUT:
# simcoeff: nsim x number of coeff matrix of coefficients for each of the simulated datasets in simfits.
checkSimCoeff2<- function(baseModel, simfits){
  # fit models to each simfits data
  simcoeff<- matrix(NA, ncol(simfits), length(baseModel$coeff))
  e1<- new.env(parent=environment(update))
  for(i in 1:ncol(simfits)){
    if((i/10)%%1 == 0 ){cat(i, '\n')}else{cat('.')}
    fit<-updateResp(baseModel, simfits[,i])
    # store coeffs
    simcoeff[i,]<- fit$coefficients
  }
  return(simcoeff)
} # end func

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
# ~~~ checkSimCoeff_cis - own update + CI calc ~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to check the true coeff lie within 95% of sim coefficients and get cis of simulated models
# same as checkSimCoeff func above but uses my own function for updating the glm formulae so that the new responses dont have to be in the workspace. 
# INPUT:
# baseModel:  model fitted object such as glm
# simfits:    n x nsim matrix of values.  each column is a simulated dataset (from a poisson)
# OUTPUT (list):
# simcoeff: nsim x number of coeff matrix of coefficients for each of the simulated datasets in simfits.
# CIs:      array of confidence intervals for each coeff for each simulated model (nbetas x 2 x nsims)
checkSimCoeff_cis<- function(baseModel, simfits){
  # fit models to each simfits data
  simcoeff<- matrix(NA, ncol(simfits), length(baseModel$coeff))
  CIs<- array(NA, c(length(baseModel$coefficients), 2, ncol(simfits)))
  
  for(i in 1:ncol(simfits)){
    if((i/10)%%1 == 0 ){cat(i, '\n')}else{cat('.')}
    fit<-updateResp(baseModel, simfits[,i])
    # store coeffs
    simcoeff[i,]<- fit$coefficients
    CIs[,,i]<- makeCIs(fit)    
  }
  return(list(simcoeff=simcoeff, CIs=CIs))
} # end func

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~ makeCIs - get cis for each coeff ~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# uses usual ci formula to get 95% cis for each parameter in model
makeCIs<- function(fit){
  cis<- matrix(NA, length(fit$coefficients), 2)
  for(i in 1:nrow(cis)){
    cis[i,1]<-fit$coefficients[i] + qnorm(0.025)*summary(fit)$coef[i,2]
    cis[i,2]<-fit$coefficients[i] - qnorm(0.025)*summary(fit)$coef[i,2]  
  }
  return(cis)
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# work out how many intervals contain truth
checkBetaCIs<-function(coef, cis){
  marker<- rep(0, ncol(cis))
  for(i in 1:ncol(cis)){
    marker[i]<-ifelse(coef< cis[2, i] & coef>cis[1, i], 1, 0)
  }
  return(marker)
}