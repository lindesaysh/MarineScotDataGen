# Data generation for marine scotland tender
# 11/2/13
# change 3

# load data
load("~/Dropbox/MarineScotland_LSH/data/three-phases-of-rodsand.RData")
source('~/Dropbox/MarineScotland_LSH/MarineScotDataGen/code/Functions.R')

load("~/../Dropbox/MarineScotland_LSH/data/three-phases-of-rodsand.RData")
source('~/../Dropbox/MarineScotland_LSH/MarineScotDataGen/code/Functions.R')

data<-danish.abc
plot(data$X, data$y, pch=20, cex=0.2, asp=1)

# make grid for prediction
# load existing pred grid as it has depth in it?
grid<-read.csv('../data/PredictionGrid.csv')

# fit model
baseModel<- glm(ducks ~ Depth + phase, family=poisson, data=data)
coeffs<- baseModel$coefficients
coefCIs<- makeCIs(baseModel)

# make sim data and get sim coeffs
nsim=500
test.simdata<- getSimData(baseModel, nsim)
#test.coef<-checkSimCoeff(baseModel, test.simdata)
#test.coef2<-checkSimCoeff2(baseModel, test.simdata)

test.coef.cis<-checkSimCoeff_cis(baseModel, test.simdata)
#benchmark(checkSimCoeff(baseModel, test.simdata), checkSimCoeff2(baseModel, test.simdata), replications=0)

cicheck<- matrix(NA, nsim, length(baseModel$coefficients))
for(i in 1:length(baseModel$coefficients)){
  cicheck[,i]<-checkBetaCIs(baseModel$coefficients[i], test.coef.cis$CIs[i,,])  
}
betacicheck<-apply(cicheck,2,sum)/nsim

# make plots of intervals for each coeff
png(file='results/simplePoiss_betas.png')
par(mfrow=c(2,2))
for(i in 1:length(baseModel$coeff)){
  plotCoeffIntervals(test.coef.cis$simcoeff[,i],c(coeffs[i], coefCIs[i,]), paste('beta', i-1, ': ', betacicheck[i],sep=''))
}
dev.off()





