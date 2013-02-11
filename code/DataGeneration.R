# Data generation for marine scotland tender
# 11/2/13
# change 2

# load data
load("~/Dropbox/MarineScotland/data/three-phases-of-rodsand.RData")

# fit model
baseModel<- glm(ducks ~ Depth + phase, family=poisson, data=danish.abc)

# make sim data and get sim coeffs
test.simdata<- getSimData(baseModel, 10)
test.coef<-checkSimCoeff(baseModel, test.simdata)

# make plots of intervals for each coeff
par(mfrow=c(2,2))
for(i in 1:length(baseModel$coeff)){
  plotCoeffIntervals(test.coef[,i], baseModel$coeff[i], paste('beta', i-1, sep=''))
}
