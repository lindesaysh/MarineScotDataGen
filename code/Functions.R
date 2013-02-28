



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~ Simulate, refit, cis, preds ~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wrapperfunc<-function(model, grid, seed, preds_imp, d=NULL){
    
  # generate one set of simulated data
  simData<-addNoise(preds_imp, d=d, seed)
  
  # fit model to simData
  newfit<-updateResp(model, simData[1:nrow(grid)], grid)
  
  # fitted values to be saved
  #fits_simdata<- fitted(newfit)
  # calculate 95% CIs for each gridcell
  pred_cis<- makeCIs(newfit, grid)
  # make a marker to say which CIs contain preds_imp
  CIcheck<-checkCIs(exp(preds_imp)[1:nrow(grid)], pred_cis[,2:3])
  
  return(list(simData=simData, fit=newfit, pred_cis=pred_cis, CIcheck = CIcheck))
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~ addNoise - poisson sim ~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate data from a baseModel, nsim times
# INPUT:
# baseModel:  model fitted object such as glm
# OUTPUT:
# simfits:  n x nsim matrix of values.  each column is a simulated dataset (from a poisson)
addNoise<- function(data, d=NULL, seed){
  # simulate data. make simfits which is nsim x nrow(predgrid)
  n<- length(data)
  if(is.null(d)){
    set.seed(i)
    simdata<-rpois(n, data)
  }else{
    set.seed(i)
    simdata<-rpois.od(n, lambda=exp(data), d=d)  
  }
  return(simdata)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ rpois.od - generate overdispersed poisson ~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INPUTS:
  # n:      how many points to generated (same length as lambda)
  # lambda: vector of lambdas to generate from
  # d:      scaling parameter to increase/decrease overdispersion
# OUTPUTS:
  # data length n of generations from an overdisperesed poisson distribution
rpois.od<- function(n, lambda, d=1){
  if(d==1){
    rpois(n, lambda)
  }else{
    rnbinom(n, size=(lambda/(d-1)), mu=lambda)
  }
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~ updateResp - own resp update func ~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# updating a model with a new response vector.  
# prevents the need for update() which needs certain parameters to exist in the workspace rather than the function space.
updateResp<-function(model, newresp, data){
    
#   if(class(model)[1]=='glm'){
#     teval<-paste('fit<-glm(newresp ~', paste(attr(model$terms, "term.labels"), sep='', collapse='+'), ', family=',model$family[1],', data=data)')  
#   }
#   if(class(model)[1]=='gam'){
#     teval<-paste('fit<-gam(newresp ~', strsplit(as.character(model$formula), '~')[3], ', family=',model$family[1],', data=data)')  
#   }
  form<- strsplit(as.character(model$call), '~')
  teval<-paste("fit<-", form[[1]], "(newresp~", form[[2]][2], ", family=", form[[3]], ", data=data)")
  
  
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

# makeCIs<- function(coef, df, se){
#   cis<- matrix(NA, length(coef), 2)
#   for(i in 1:nrow(cis)){
#     cis[i,1]<-coef[i] + qt(0.025, df)*se[i]
#     cis[i,2]<-coef[i] - qt(0.025, df)*se[i]
#   }
#   return(cis)
# }

makeCIs<- function(newfit, grid){
  preds<- predict(newfit, grid, type='link', se.fit=T)
  uppers<- exp(preds$fit + qt(0.975, newfit$df.residual) * preds$se.fit)
  lowers<- exp(preds$fit - qt(0.975, newfit$df.residual) * preds$se.fit)
  return(data.frame(preds = exp(preds$fit), lowers=lowers, uppers=uppers))  
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
checkCIs<-function(truth, cis){
  marker<- rep(0, nrow(cis))
  for(i in 1:nrow(cis)){
    marker[i]<-ifelse(truth[i]< cis[i,2] & truth[i]>cis[i,1], 1, 0)
  }
  return(marker)
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~ make predictions for a model ~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
makePredictions<- function(model, grid){
  #preds <- (coef[-1] %*% grid) + coef[1]
  preds<- predict(model, grid, type='response')
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



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~ adding hotspot to data ~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# can be used to add an impact effect over the surface
# INPUT:
# inputdat: data for hotspot to be added to (vector)
# x, y :    single numbers representing the coordinates 
# altitude: height of hotspot to add (eg. numbers of animals)
# sigma:    spread of hotspot (related to range of distances in the grid area)
# OUTPUT:
# density:  vector of inputdat with the hotspot added
#
addHotspot<-function (inputdat, x, y, altitude, sigma, add=F){ 
  if (!is.numeric(x) | !is.numeric(y)) 
    stop("\nThe x/y-position must be numeric.\n")
  if (!is.numeric(altitude)) 
    stop("\nThe altitude must be numeric.\n")
  if (!is.numeric(sigma)) 
    stop("\n<sigma> must be numeric.\n")
  if (sigma <= 0) 
    stop("\nsigma must be greater then zero.\n")
  
  coords <- matrix(0, nrow = nrow(grid), ncol = 2)
  coords[, 1] <- grid[,2]
  coords[, 2] <- grid[,3]
  distance <- sqrt((x - coords[, 1])^2 + (y - coords[, 2])^2)
  hs <- dnorm(distance, mean = 0, sd = sigma)
  hs.max <- dnorm(0, mean = 0, sd = sigma)
  hs <- hs/hs.max
  hs <- hs * altitude
  if(add==T){
    density <- inputdat + hs 
  }else{
    density <- inputdat - hs
  }
  
  #if(min(density)<0) density[which(density<0)]<-0
  
  return(list(density=density, hotspot = hs))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~ adding impact effect to data ~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

addImpact<- function(data, Impact){
  if(Impact$impact=='decrease'){
    impEffect<- log(Impact$effect)
      simgen_imp<-c(data, data - impEffect)
  }
  if(Impact$impact=='redistribution'){
    impEffect<- addHotspot(data, Impact$x, Impact$y, Impact$altitude, Impact$sigma)
    simgen_imp<- c(data, impEffect$density)
  }
  return(simgen_imp)  
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~ adding correlation to data ~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getCorrelation<-function(blockid, rho, sdTuning){
  # define blocks - use survey day
  # get length for each block
  # nblocks<- Autocorr$nblocks
  # blockLength<- ceiling(nrow(grid)/nblocks)
  # blockid<- rep(1:nblocks, each=blockLength)[1:nrow(grid)]
  blocks<-table(blockid)
  
  dataCor<- NULL
  for(i in 1:length(blocks)){
    # calc arim.sim() for each block
    blockCor<-arima.sim(list(ar=rho), n=blocks[i], sd = sdTuning)
    dataCor<- c(dataCor, blockCor)
  }
  return(dataCor)
}





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~ making local radial for CReSS ~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LocalRadialFunction<- function(radiusIndices, dists, radii,aR){
  for (i in 1:length(aR)){
    zhold<- dists[,aR[i]]
    zhold<- exp(-zhold/(radii[radiusIndices[i]]**2)) 
    if (i==1) {B<- zhold} else {B<- cbind(B,zhold)}
  }
  B <- data.frame(B)
  for(i in 1:ncol(B)){
    names(B)[i] <- paste('b', i, sep='')
  }
  B=as.matrix(B)
  return(B)
}