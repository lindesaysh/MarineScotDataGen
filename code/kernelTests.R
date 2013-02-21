# kernSmooth
require(KernSmooth)
data(geyser, package='MASS')
x<- cbind(geyser$duration, geyser$waiting)
est<-bkde2D(x, bandwidth=c(0.7,7))
contour(est$x1, est$x2, est$fhat)
persp(est$fhat)


# ks package using kde

require(ks)

x<- data[data$phase=='A',c(5, 6, 9)]
H.pi<- Hpi(x, binned=T)
fhat <- kde(x=x, H=H.pi, binned = T)  
plot(fhat)


fhat<-kde(x, H=Hpi(x, binned=T), binned=T, compute.cont=T)
plot(fhat, display='density', col=heat.colors(100)zlim=)

require(fields)
quilt.plot(fhat$eval.points[[1]], fhat$eval.points[[2]], fhat$eval.points[[3]])


bob<-fhat$estimate
bob2<-matrix(NA, 51,51)
for(i in 1:51){
  for(j in 1:51){
    bob2[i, j]<-sum(bob[i,j,])
  }
}

image.plot(fhat$eval.points[[1]], fhat$eval.points[[2]],bob2)


#smooth.2d from the fields package
bob<-smooth.2d(Y=data$ducks, x=cbind(data$X, data$Y))

#np package
require(np)
xdat<- data[data$phase=='A',c(5, 6, 9)]
#attach(xdat)
#bw <- npregbw(ducks~ X + Y + factor(phase) + Depth, regtype="ll", bwmethod="cv.aic")

bw <- npregbw(ducks~ X+Y, regtype="ll", bwmethod="cv.aic", data=xdat)
plot(bw, plot.errors.method="asymptotic")
points(age, logwage, cex=.2, col="red")


data("cps71")
attach(cps71)

bw <- npregbw(logwage~age, regtype="ll", bwmethod="cv.aic")

# Next, plot the regression function...

plot(bw, plot.errors.method="asymptotic")
points(age, logwage, cex=.2, col="red")

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)

# Next, plot the derivative...

plot(bw, plot.errors.method="asymptotic", gradient=TRUE)

detach(cps71)

# Sleep for 5 seconds so that we can examine the output...

Sys.sleep(5)