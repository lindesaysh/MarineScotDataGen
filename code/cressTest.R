# fitting cress models


require(splines)
require(akima)
bd = c(min(data$Depth, grid$Depth), max(data$Depth, grid$Depth))
knots = quantile(grid$Depth, c(0.25,0.5,0.75))
bob<-bs(data$Depth, knots=knots, Boundary.knots=bd)

distGrid<-matrix(as.matrix(dist(cbind(grid$X, grid$Y)), nrow=nrow(grid)), nrow=nrow(grid))

# select knots - 20 of them (roughly what gam was picking)
set.seed(5)
sampleid<- sample(1:nrow(data), 1000)
temp<-cover.design(R=data[sampleid, 5:6], nd=50, nruns=5)$best.id  
knotid<-(1:nrow(data))[sampleid[temp]]
knotgrid<-cbind(data$X[knotid], data$Y[knotid])

plot(data$X, data$Y)
points(data$X[sampleid], data$Y[sampleid], pch=20, col='blue')
points(knotgrid, pch=20, col='red')

# ~~~~~~~~~~~~~~~~~~~~~~~~~
#for fitting
# use interpp function to get distances for data to knots
# first cut interpolate grid to knots
g2k<- matrix(NA, nrow=length(knotid), ncol=nrow(grid))
for(i in 1:ncol(distGrid)){
  cat(paste('Calculating iteration:', i, '\n', sep=' '))
  g2k[,i]<-interpp(grid$X, grid$Y, distGrid[,i], data$X[knotid], data$Y[knotid], duplicate='mean')$z
}

# check distances
quilt.plot(grid$X, grid$Y, g2k[1,], asp=1, nrow=104, ncol=55)
quilt.plot(grid$X, grid$Y, g2k[10,], asp=1, nrow=104, ncol=55)
quilt.plot(grid$X, grid$Y, g2k[20,], asp=1, nrow=104, ncol=55)

# second interp from grid2knots to data2knots
d2k<- matrix(NA, nrow=nrow(data), ncol=length(knotid))
for(i in 1:nrow(g2k)){
  cat(paste('Calculating iteration:', i, '\n', sep=' '))
  d2k[,i]<-interpp(grid$X, grid$Y, g2k[i,], data$X, data$Y, duplicate='mean')$z
}

# check distances
quilt.plot(data$X, data$Y, d2k[,1], asp=1, nrow=104, ncol=55)
quilt.plot(data$X, data$Y, d2k[,10], asp=1, nrow=104, ncol=55)
quilt.plot(data$X, data$Y, d2k[,20], asp=1, nrow=104, ncol=55)
# ~~~~~~~~~~~~~~~~~~~~~~
# calculate radii
NumberOfRadii = 50
#establish smallest knot-knot and observation-knot distances
rmin<- sqrt(max(distGrid)/21)
rmax<- sqrt(max(distGrid)/3e-7)
r_seq <- exp(seq(log(rmin), log(rmax), length=NumberOfRadii))
r_seq
r_seq[30]
# ~~~~~~~~~~~~~~~~~~~~~~
# 1d setup  
bd = c(min(data$Depth, grid$Depth), max(data$Depth, grid$Depth))
knots = quantile(grid$Depth, c(0.25,0.5,0.75))
newknots<- cover.design(R=knotgrid, nd=10, nruns=5)$best.id  
dist<-d2k
baseModel<-glm(ducks ~ bs(Depth, knots=knots, Boundary.knots=bd) + LocalRadialFunction(rep(1, 10), dist, 10619, newknots), family=quasipoisson, data=data)

dist<-t(g2k)
preds_temp<- predict(baseModel, cbind(grid, t(g2k)), type='link', se.fit=T)
quilt.plot(grid$X, grid$Y, exp(preds_temp$fit), asp=1, nrow=104, ncol=55)

preds_imp<-addImpact(preds_temp$fit, Impact)
idpre<- 1:nrow(grid); idpost<- (nrow(grid)+1): length(preds_imp)
par(mfrow=c(1,2))
quilt.plot(grid$X, grid$Y, exp(preds_imp)[idpre],asp=1, nrow=104, ncol=55, main='pre', zlim=range(exp(preds_imp)))
quilt.plot(grid$X, grid$Y, exp(preds_imp)[idpost],asp=1, nrow=104, ncol=55, main='post', zlim=range(exp(preds_imp)))


#~~~~~~~~~~~~~~~
nsim=500
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

round(cbind(baseModel$coeff, dataObj$fit$coeff))

summary(baseModel)$dispersion
summary(dataObj$fit)$dispersion

mean(data$ducks)
var(data$ducks)
mean(exp(preds_imp[idpre]))
var(exp(preds_imp[idpre]))
mean(apply(simData[idpre,], 2, mean))
mean(apply(simData[idpre,], 2, var))
mean(apply(preds[,1,], 2, mean))
mean(apply(preds[,1,], 2, var))


par(mfrow=c(2,2))
j=1
quilt.plot(grid$X, grid$Y, simData[idpre,j], asp=1, nrow=104, ncol=55, main='noisy data: pre')
quilt.plot(grid$X, grid$Y, preds[,1,j], asp=1, nrow=104, ncol=55, main='simdatapreds: pre, i=2')
points(knotgrid[newknots,], pch=20, col='red')
quilt.plot(grid$X, grid$Y, apply(markers,1,mean), asp=1, nrow=104, ncol=55, main='CIcheck', zlim=c(0,1))
# relative % bias
relbias<-matrix(NA, nrow(preds), dim(preds)[3])
for(i in 1:nsim){
  relbias[,i]<-(preds[,1,i]-exp(preds_imp[idpre]))/exp(preds_imp[idpre])
}
quilt.plot(grid$X, grid$Y, apply(relbias,1,median), asp=1, nrow=104, ncol=55, main='% relative bias')

# plot again with no limits on z
par(mfrow=c(1,1))
quilt.plot(grid$X, grid$Y, apply(markers,1,mean), asp=1, nrow=104, ncol=55, main='CIcheck')
ciid<- which(apply(markers,1,mean)<.9)
points(grid$X[ciid], grid$Y[ciid], pch=20, cex=0.8, col='white')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# partial plots

# get model matrix
bob<-model.matrix(dataObj$fit)

# get coeffs for particular covariate and subset model matrix accordingly
betas<- matrix(as.matrix(dataObj$fit$coeff, ncol=1), ncol=1)

# model matrix times coeffs + intercept? on scale of link
d_id<- c(2:7)
space_id <- c(1, 8:length(betas))

d_partial<- resid(dataObj$fit, type='pearson') + exp(bob[,d_id]%*%(betas[d_id]))
space_partial<- resid(dataObj$fit, type='pearson') + exp(bob[,space_id]%*%betas[space_id])

par(mfrow=c(1,2))
plot(grid$Depth, d_partial, pch=20)
#lines(lowess(grid$Depth, d_partial), col='red', lwd=2)
quilt.plot(grid$X, grid$Y, space_partial, asp=1, col=rev.terrain.colors(100), nrow=104, ncol=55)
