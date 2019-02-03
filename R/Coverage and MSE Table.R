library(causalForest)
library(mgcv)
library(randomForestCI)
library(FNN)
library(Hmisc)
library(xtable)

rm(list = ls())

n = 5000
s = 2500
ntree = 2000
sigma = .1

n.test = 1000

dvals = c(2, 3, 4, 5, 6, 8)
simu.reps = 40

effect = function(x) {
  x[1] + x[2]
}
baseline <- function(x){10}

simu.fun = function(d) {
  X = matrix(runif(n * d, 0, 1), n, d) # features
  W = rbinom(n, 1, 0.5) #treatment condition
  odds = apply(X, 1,baseline) +  (W - 0.5) * apply(X, 1, effect) + sigma * rnorm(n)
  Y = rbinom(n, 1, odds/(1+odds))

  Y.fit.1 <- SuperLearner(Y=Y[W==1], X=as.data.frame(X)[W==1,], SL.library=c("SL.glm", "SL.earth"), family="binomial")
  Y.fit.0 <- SuperLearner(Y=Y[W==0], X=as.data.frame(X)[W==0,], SL.library=c("SL.glm", "SL.earth"), family="binomial")
  Y.hat.1 <- predict(Y.fit.1,X[W==1,])$pred
  Y.hat.0 <- predict(Y.fit.0,X[W==0,])$pred
  Y.hat <-  Y
  Y.hat[W==1] <- Y.hat.1
  Y.hat[W==0] <- Y.hat.0
  Y.hat <- Y.hat/(1-Y.hat)


  X.test = matrix(runif(n.test * d, 0, 1), n.test, d)
  true.eff = apply(X.test, 1, effect)


  tauhat.rf <- predict(forest, X.test)

  #
  # random forest
  #
  forest = grf::causal_forest(X, Y.hat, W,  num.trees = ntree)
  tauhat.rf = predict(forest, X.test,estimate.variance=TRUE)
  predictions = tauhat.rf$predictions

  se.hat = sqrt(tauhat.rf$variance.estimates)
  rf.cov = abs(predictions - true.eff) <= 1.96 * se.hat
  rf.covered = mean(rf.cov)
  rf.mse = mean((predictions - true.eff)^2)

  c(rf.covered = rf.covered,
    rf.mse = rf.mse)
}

results.raw = lapply(dvals, function(d) {
  print(paste("NOW RUNNING:", d))
  res.d = sapply(1:simu.reps, function(iter) simu.fun(d))
  res.fixed = data.frame(t(res.d))
  print(paste("RESULT AT", d, "IS", colMeans(res.fixed)))
  res.fixed
})

results.condensed = lapply(results.raw, function(RR) {
  RR.mu = colMeans(RR)
  RR.var = sapply(RR, var) / (nrow(RR) - 1)
  rbind("mu"=RR.mu, "se"=sqrt(RR.var))
})

results.condensed

save.image("table2_easy_sigmoid.RData")

results.parsed = lapply(results.condensed, function(RR) {
  apply(RR, 2, function(arg) {
    paste0(round(arg[1], 2), " (", round(100 * arg[2], 0), ")")
  })
})

results.table = data.frame(cbind(d=dvals, Reduce(rbind, results.parsed)))
results.table

results.table = results.table[,c(1, 3, 5, 7, 2, 4, 6)]
xtab = xtable(results.table)
print(xtab, include.rownames = FALSE)
