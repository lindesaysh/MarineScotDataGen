arima(baseModel$residuals, order=c(1,1,0))

Call:
  arima(x = baseModel$residuals, order = c(1, 0, 0))

  Coefficients:
         ar1  intercept
      0.0283    -0.0293
s.e.  0.0056     0.0605

sigma^2 estimated as 108.9:  log likelihood = -118582.8,  aic = 237171.6

sim<-arima.sim(n=nrow(data), list(order=c(1,0,0), ar=0.0283), sd=sqrt(108.9))
ts.plot(sim[1:100])
acf(sim[1:1000], lag=10)
acf(baseModel$residuals, lag=10)



N = 11478
rho = 0.6
log.lambda<- 1+ arima.sim(list(ar=rho), n=N)
y<- rpois(N, lambda = exp(log.lambda))
mean(y)
var(y)
hist(y)
acf(y)


modelPreds<- predict(baseModel, grid, type='link')



#ts.plot(dataCor)
range(dataCor)

# data to simulate from
# also add at this point an impact effect.
simgen<-modelPreds+dataCor




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



require(fields)
par(mfrow=c(2,2))
quilt.plot(data$X, data$Y, data$ducks, nrow=104, ncol=55, xlim=range(grid$longitude), ylim=range(grid$latitude), asp=1, main='raw')
quilt.plot(grid$longitude, grid$latitude, exp(modelPreds), asp=1, nrow=104, ncol=55, main='preds')
quilt.plot(grid$longitude, grid$latitude, exp(simgen), asp=1, nrow=104, ncol=55, main='preds+corr')

n<- length(simgen)
set.seed(1)
#simdata<-rpois(n, exp(simgen))
#quilt.plot(grid$longitude, grid$latitude, simdata, asp=1, nrow=104, ncol=55, main='rpois')
d<- Autocorr$d
simdata<-rpois.od(n, lambda=exp(simgen), d=d)
quilt.plot(grid$longitude, grid$latitude, simdata, asp=1, nrow=104, ncol=55, main=paste('rpois.od, d=', d, sep=''))

mean(data$ducks)
mean(exp(simgen))
mean(simdata)
#mean(simdata2)

var(data$ducks)
var(exp(simgen))
var(simdata)
#var(simdata2)


acf(residuals(baseModel))
tempModel<-updateResp(baseModel, simdata, data=grid)
acf(residuals(tempModel))