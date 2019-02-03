library(plotly)
d = 6
# No. observations
n= 5000
# s set size
s = 2500
X = matrix(runif(n*d),ncol = d,nrow=n)
X[,6] <- rbinom(5000,1,.5)
xi <- function(x){
   return (1+1/(1+exp(-20*(x-1/3))))
 }
tau = xi(X[,1])*xi(X[,2])
tauxy <- function(X_0,X_1){
  xi(X_0)*xi(X_1)
}
# Quick plot of tau as function of X_1, X_2 assuming continuous support
plotFunc <- function(func){
  X_0 <- seq(0,1,.01)
  X_1 <- seq(0,1,.01)

  Z <- outer(X_0, X_1, func)
  plot_ly(x=~X_0, y=~X_1, z=~Z, type='contour',
          contours = list(start = 1,end = 4,showlines = FALSE,coloring = 'heatmap'))
}

plotFunc(tauxy)


subSampleMask = sample(seq(0, n), s)
setIMask = sample(subSampleMask, as.integer(s/2))
setI = data.frame(X[setIMask,])
setJMask = subSampleMask[!subSampleMask%in%setIMask]
setJ = data.frame(X[setJMask,])

require(party)
fit.data <- data.frame(tau[setJMask],setJ)
fit <- ctree(fit.data$tau.setJMask.~., data=fit.data, controls=ctree_control(minsplit=2,minbucket=2,testtype="Univariate"))
tau_hat <- predict(fit,setI)


plotResults <- function(SetI,tau_hat){
  require(mgcv)
  # Interpolate
  fit.data <- data.frame(xi=SetI[,1],yi=SetI[,2],tau_hat)
  g.spline  <- gam(tau_hat ~ s(xi, bs = "cr")+s(yi, bs = "cr"),data=fit.data)
  pred.func <- function(xi,yi){
    g.pred <- predict(g.spline, newdata = data.frame(xi,yi), type = "link",se.fit = TRUE)
    g.pred$fit
  }
  #sAEnbk1 <- g.pred$se.fit
  # Set up a regular grid of interpolation points
  xi <-  seq(min(SetI[,1]),max(SetI[,1]), .01)
  yi <-  seq(min(SetI[,2]),max(SetI[,2]), .01)
  zi <- outer(xi, yi, pred.func)

  plot_ly(x=~xi, y=~yi, z=~zi, type='contour',
          contours = list(start = 1,end = 4,showlines = FALSE,coloring = 'heatmap'))
}
plotResults(setI,tau_hat)




#grf paper
rm(list = ls())
n = 10000
ntree = 4000
n.test = 10000

sigma = 1
d = 20

effect = function(x) {
  4/((1 + exp(-12 * (x[1] - 0.5))) * (1 + exp(-12 * (x[2] - 0.5))))
}

baseline = function(x) { 0 }

X = matrix(runif(n * d, 0, 1), n, d) # features
W = rbinom(n, 1, 0.5) #treatment condition
Y = apply(X, 1, baseline) +  (W - 0.5) * apply(X, 1, effect) + sigma * rnorm(n)

X.test = matrix(runif(n.test * d, 0, 1), n.test, d)
true.eff = apply(X.test, 1, effect)

forest = grf::causal_forest(X, Y, W, num.trees = ntree)
tauhat.rf <- predict(forest, X.test)


print(mean((true.eff - tauhat.rf$predictions)^2))
minp = min(true.eff, tauhat.rf$predictions)
maxp = max(true.eff, tauhat.rf$predictions)
rngp = maxp - minp

ncol=100
true.scl = pmax(ceiling(ncol * (true.eff - minp) / rngp), 1)
rf.scl = pmax(ceiling(ncol * (tauhat.rf$predictions - minp) / rngp), 1)
hc = heat.colors(ncol)

plot(X.test[,1], X.test[,2], pch = 16, col = hc[true.scl], xlab = "", ylab = "")
dev.off()

plot(X.test[,1], X.test[,2], pch = 16, col = hc[rf.scl], xlab = "", ylab = "")
dev.off()




