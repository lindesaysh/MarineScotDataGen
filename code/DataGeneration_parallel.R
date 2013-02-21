
# load data and source functions
load("~/Dropbox/MarineScotland_LSH/data/three-phases-of-rodsand.RData")
source('~/Dropbox/MarineScotland_LSH/MarineScotDataGen/code/Functions.R')

load("~/../Dropbox/MarineScotland_LSH/data/three-phases-of-rodsand.RData")
source('~/../Dropbox/MarineScotland_LSH/MarineScotDataGen/code/Functions.R')

require(doParallel)
require(parallel)

data<-danish.abc
plot(data$X, data$y, pch=20, cex=0.2, asp=1)

# load existing pred grid as it has depth in it?
grid<- read.csv('../data/PredictionGrid.csv')
grid2<- data.frame(rbind(cbind(grid$Depth, rep('A', nrow(grid))), cbind(grid$Depth, rep('B', nrow(grid))), cbind(grid$Depth, rep('C', nrow(grid)))))
names(grid2)<- c('Depth', 'phase')
class(grid2$Depth)<- 'numeric'
                                  

# fit model
baseModel<- glm(ducks ~ Depth + as.factor(phase), family=poisson, data=data)

nsim=100

# generate data 
preds<- matrix(NA, nrow(grid2), nsim) 
simData <- matrix(NA, nrow(data), nsim)
markers = coeff = matrix(NA, nsim, length(baseModel$coefficients))

n.cores<-detectCores()
cl<- makeCluster(n.cores)
registerDoParallel(cl)
system.time(simData<-foreach(i = 1:nsim) %dopar% wrapperfunc(baseModel, grid2, seed=i))

  #if((i/10)%%1 == 0 ){cat(i, '\n')}else{cat('.')}
#  simData <- dataObj$simData
  #preds[,i] <-  dataObj$preds
  #markers[i,] <- dataObj$betaCheck
  #coeff[i,] <- dataObj$fit$coeff
 # rm(dataObj)
#}


betacicheck<-apply(markers,2,sum)/nsim
baseMcoefCIs<- makeCIs(baseModel$coef, df=baseModel$df.residual, se=summary(baseModel)$coef[,2])
par(mfrow=c(2,2))
for(i in 1:ncol(coeff)){
  plotCoeffIntervals(data=coeff[,i], c(baseModel$coefficients[i], baseMcoefCIs[i,]), paste('beta', i-1, ': ', betacicheck[i],sep=''))
}
 
require(fields)
par(mfrow=c(1,1))
quilt.plot(grid$longitude, grid$lat, preds[1:nrow(grid),1], asp=1, nrow=100, ncol=50)
quilt.plot(grid$longitude, grid$lat, preds[(nrow(grid)+1):(nrow(grid)*2),1], asp=1, nrow=100, ncol=50)
quilt.plot(grid$longitude, grid$lat, preds[((nrow(grid)*2)+1):nrow(grid2),1], asp=1, nrow=100, ncol=50)
  
  