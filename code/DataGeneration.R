
# load data and source functions
load("~/Dropbox/MarineScotland_LSH/data/three-phases-of-rodsand.RData")
source('~/Dropbox/MarineScotland_LSH/MarineScotDataGen/code/Functions.R')

load("~/../Dropbox/MarineScotland_LSH/data/three-phases-of-rodsand.RData")
source('~/../Dropbox/MarineScotland_LSH/MarineScotDataGen/code/Functions.R')


data<-danish.abc[danish.abc$phase=='A',]
plot(data$X, data$y, pch=20, cex=0.2, asp=1)
quilt.plot(data$X, data$Y, data$ducks,asp=1, nrow=104, ncol=55, main='raw data')

# load existing pred grid as it has depth in it?
grid<- read.csv('../data/PredictionGrid.csv')
names(grid)[2:3]<- c('X', 'Y')
#grid2<- data.frame(Depth=grid$Depth, phase = rep('A', nrow(grid)))
#grid2<- rbind(grid2, data.frame(Depth=grid$Depth, phase = rep('B', nrow(grid))))
#grid2<- rbind(grid2, data.frame(Depth=grid$Depth, phase = rep('C', nrow(grid))))
                                 

# fit model
# simple glm
#baseModel<- glm(ducks ~ Depth , family=poisson, data=data)

# simple gam
require(mgcv)
#baseModel<- gam(ducks ~ Depth  , family=poisson, data=data)

# gam with a smooth on depth
baseModel<- gam(ducks ~ s(Depth,k=4) + s(X,Y, k=10), family=quasipoisson, data=data)


# simulation setup
# Impact list of parameters
Impact<- list(impact = "decrease", effect=1.5)
#Impact<- list(impact = "redistribution", x = mean(grid[,2]), y = mean(grid[,3]), altitude= log(mean(data$ducks) + 0.5), sigma=5000)
# ~~~~~~~~~~~~~~~~~~

# make predictions and add impact
# make predictions
preds<- predict(baseModel, grid, type='link', se.fit=T)
# add impact effect
preds_imp<-addImpact(preds$fit, Impact)
idpre<- 1:nrow(grid); idpost<- (nrow(grid)+1): length(preds_imp)
par(mfrow=c(1,2))
quilt.plot(grid$X, grid$Y, exp(preds_imp)[idpre],asp=1, nrow=104, ncol=55, main='pre', zlim=range(exp(preds_imp)))
quilt.plot(grid$X, grid$Y, exp(preds_imp)[idpost],asp=1, nrow=104, ncol=55, main='post', zlim=range(exp(preds_imp)))

#~~~~~~~~~~~~~~~
nsim=50
#newgrid<-data.frame(rbind(grid, grid))
# generate data 
preds<- array(NA, c(nrow(grid), 3, nsim))
simData <- matrix(NA, nrow(grid)*2, nsim)
markers <- matrix(NA, nrow(grid), nsim)
#markers = matrix(NA, nsim, nrow(grid))
#tim<-proc.time()
for(i in 1:nsim){

  if((i/10)%%1 == 0 ){cat(i, '\n')}else{cat('.')}
  
  dataObj<-wrapperfunc(baseModel, grid, seed=i, preds_imp, d=30)

  simData[,i] <- dataObj$simData
  preds[,,i] <-  matrix(as.matrix(dataObj$pred_cis, ncol=2), ncol=3)
  markers[,i] <- dataObj$CIcheck
  #coeff[i,] <- dataObj$fit$coeff
  # rm(dataObj)
}

plot(baseModel)
plot(dataObj$fit)
summary(baseModel)$dispersion
summary(dataObj$fit)$dispersion


mean(data$ducks)
var(data$ducks)
mean(exp(preds_imp[idpre]))
var(exp(preds_imp[idpre]))

mean(apply(simData[idpre,], 2, mean))
#apply(simData[idpost,], 2, mean)
mean(apply(simData[idpre,], 2, var))
#apply(simData[idpost,], 2, var)
#apply(simData[idpre,], 2, sum)
#apply(simData[idpost,], 2, sum)
mean(apply(preds[,1,], 2, mean))
mean(apply(preds[,1,], 2, var))

par(mfrow=c(2,2))
j=1
quilt.plot(grid$X, grid$Y, simData[idpre,j], asp=1, nrow=100, ncol=50, main='noisy data: pre')
#quilt.plot(grid$X, grid$Y, simData[idpost,2], asp=1, nrow=100, ncol=50, main='noisy data: post')
quilt.plot(grid$X, grid$Y, preds[,1,j], asp=1, nrow=100, ncol=50, main='simdatapreds: pre')
points(knotgrid[newknots,], pch=20, col='red')
#quilt.plot(grid$X, grid$Y, preds[,2,j], asp=1, nrow=100, ncol=50, main='simdatapreds: pre lower')
#quilt.plot(grid$X, grid$Y, preds[,3,j], asp=1, nrow=100, ncol=50, main='simdatapreds: pre upper')

# no impact effect in model so cant pick up difference...
# ranges are so different because of noise so confidence intervals dont match up at all

quilt.plot(grid$X, grid$Y, apply(markers,1,mean), asp=1, nrow=104, ncol=55, main='CIcheck', zlim=c(0,1))
#quilt.plot(grid$X, grid$Y, markers[,3], asp=1, nrow=100, ncol=50, main='CIcheck')

# relative % bias
quilt.plot(grid$X, grid$Y, (apply(preds[,1,], 1, mean)-exp(preds_imp[idpre]))/exp(preds_imp[idpre]), asp=1, nrow=104, ncol=55, main='% relative bias')

#tim - proc.time()

# ~~~~~~~~~~~~~~~~~~~~~~~
# check for whether CI inclusion of truth is positively or negatively biased
temp<- vector(length=nrow(grid))
j=100
for(i in 1:nrow(grid)){
  if(exp(preds_imp)[idpre][i]< preds[i,3,j] & exp(preds_imp)[idpre][i]>preds[i,2,j]) temp[i] <-0
  if(exp(preds_imp)[idpre][i]< preds[i,2,j])  temp[i] <- -1
  if(exp(preds_imp)[idpre][i]> preds[i,3,j])  temp[i] <- 1
}

hist(temp)
quilt.plot(grid$X, grid$Y, temp, asp=1, nrow=104, ncol=55, main='CIcheck')

# ~~~~~~~~~~~~~~~

# sample data
# use unique x values from Dato==200000404, round to be better
xid<-which(data$Dato==20000404)
Tid<-unique(signif(data$X[xid], digits=2.5))
Tid
plot(data$X[xid], data$Y[xid], pch=20)
points(Tid, rep(6045000, length(Tid)), pch=20, col='red')

#tidy up
Tid<-Tid[-c(1,21)]
points(Tid, rep(6045000, length(Tid)), pch=20, col='green')

# ask which rounded X's == Tid
t1<-which(signif(grid$X, digits=2.5)==Tid[1])
points(grid$X[t1], grid$Y[t1], pch=20, col='blue')

transectid = whichtransect = NULL
for(i in 1:length(Tid)){  
  t<-which(signif(grid$X, digits=2.5)==Tid[i])
  transectid<- c(transectid,t)
  whichtransect<- c(whichtransect, rep(i, length(t)))  
}

points(grid$X[transectid], grid$Y[transectid], pch=20, col='red')

# or interpolate to get points??


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add autocorrelation to sample

dataCorr<-getCorrelation(blockid, rho=0.6, sdTuning = 0.2)
                                                