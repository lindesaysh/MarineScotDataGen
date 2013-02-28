
# load data and source functions
load("~/Dropbox/MarineScotland_LSH/data/three-phases-of-rodsand.RData")
source('~/Dropbox/MarineScotland_LSH/MarineScotDataGen/code/Functions.R')

load("~/../Dropbox/MarineScotland_LSH/data/three-phases-of-rodsand.RData")
source('~/../Dropbox/MarineScotland_LSH/MarineScotDataGen/code/Functions.R')


data<-danish.abc[danish.abc$phase=='A',]
plot(data$X, data$y, pch=20, cex=0.2, asp=1)

# load existing pred grid as it has depth in it?
grid<- read.csv('../data/PredictionGrid.csv')

#grid2<- data.frame(Depth=grid$Depth, phase = rep('A', nrow(grid)))
#grid2<- rbind(grid2, data.frame(Depth=grid$Depth, phase = rep('B', nrow(grid))))
#grid2<- rbind(grid2, data.frame(Depth=grid$Depth, phase = rep('C', nrow(grid))))
                                 

# fit model
# simple glm
baseModel<- glm(ducks ~ Depth , family=poisson, data=data)

# simple gam
require(mgcv)
baseModel<- gam(ducks ~ Depth  , family=poisson, data=data)

# gam with a smooth on depth
baseModel<- gam(ducks ~ s(Depth, k=4) , family=poisson, data=data)



# simulation setup
# Impact list of parameters
Impact<- list(impact = "none", x = mean(grid[,2]), y = mean(grid[,3]), altitude= mean(data$ducks) + 0.5, sigma=5000)
# Autocorr list of parameters
Autocorr <- list( rho = 0.6, sdTuning = 0.02, nblocks = length(table(data$Dato)), d=20)


nsim=5

# generate data 
preds<- matrix(NA, nrow(grid), nsim) 
simData <- matrix(NA, nrow(data), nsim)
markers = coeff = matrix(NA, nsim, length(baseModel$coefficients))
#tim<-proc.time()
for(i in 1:nsim){

  if((i/10)%%1 == 0 ){cat(i, '\n')}else{cat('.')}
  
  dataObj<-wrapperfunc(baseModel, grid, data, seed=i)

  simData[,i] <- dataObj$simData
  preds[,i] <-  dataObj$preds
  markers[i,] <- dataObj$betaCheck
  coeff[i,] <- dataObj$fit$coeff
  rm(dataObj)
}

#tim - proc.time()

# plot out betas to check for biases
betacicheck<-apply(markers,2,sum)/nsim
if(class(baseModel)[1]=='glm'){se<-summary(baseModel)$coef[,2]}
if(class(baseModel)[1]=='gam'){se<-summary(baseModel)$se}
baseMcoefCIs<- makeCIs(baseModel$coef, df=baseModel$df.residual, se=se)
png('results/gam_sDepthbetacheck.png')
par(mfrow=c(3,2))
for(i in 1:ncol(coeff)){
  plotCoeffIntervals(data=coeff[,i], c(baseModel$coefficients[i], baseMcoefCIs[i,]), paste('beta', i-1, ': ', betacicheck[i],sep=''))
}
dev.off()

# check predicted data is sensible
require(fields)
i=2
predrange<-range(preds[,i])
par(mfrow=c(2,2))
quilt.plot(grid$longitude, grid$lat, preds[1:nrow(grid),i], asp=1, nrow=100, ncol=50, zlim=predrange)
quilt.plot(grid$longitude, grid$lat, preds[(nrow(grid)+1):(nrow(grid)*2),i], asp=1, nrow=100, ncol=50, zlim=predrange)
quilt.plot(grid$longitude, grid$lat, preds[((nrow(grid)*2)+1):nrow(grid2),i], asp=1, nrow=100, ncol=50, zlim=predrange)
  
  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~ adding to generated data ~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# take away a hotspot
predHot<-addHotspot(preds[1:nrow(grid),1], x = mean(grid[,2]), y = mean(grid[,3]), altitude= mean(data$ducks) + 0.5, sigma=5000, add=F)

png('results/addhotspot.png')
par(mfrow=c(2,2))
quilt.plot(grid[,2], grid[,3], preds[1:nrow(grid),1],asp=1, nrow=100, ncol=50)
quilt.plot(grid[,2], grid[,3], predHot$hotspot, asp=1, nrow=100, ncol=50)
quilt.plot(grid[,2], grid[,3], predHot$density, asp=1, nrow=100, ncol=50)
dev.off()                                                  
                                                    
                                                