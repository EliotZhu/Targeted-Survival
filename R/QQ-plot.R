library(causalForest)

rm(list = ls())


n = 2000
ntree = 3000
n.test = 1000
sigma = .1
d = 10
simu.reps=20

effect = function(x) {
  4/((1 + exp(-12 * (x[1] - 0.5))) * (1 + exp(-12 * (x[2] - 0.5))))
}
baseline = function(x) { 0 }

single.run = function(n) {
  print("estimating")
  X = as.data.frame(matrix(runif(n * d, 0, 1), n, d)) # features
  X.test = matrix(runif(n.test * d, 0, 1), n.test, d)
  W = rbinom(n, 1, 1/(1+exp(-X[,1]))) #treatment condition
  Z = apply(X, 1, baseline) +  (W - 0.5) * apply(X, 1, effect) + sigma * rnorm(n)
  pr = 1/(1+exp(-Z)) #inv-logit function
  Y = rbinom(n, 1, pr)
  pihat <- SuperLearner(W, X=X, SL.library=c("SL.glm", "SL.gam"), family="binomial")
  pihat <- pihat$SL.predict
  #causal forest
  Y.fit <- SuperLearner(Y=Y, X=X, SL.library=c("SL.glm", "SL.gam","SL.rpart"), family="binomial")
  Y.hat <- Y.fit$SL.predict
  forest = grf::causal_forest(X, Y, W,W.hat = pihat,Y.hat=Y.hat,  num.trees = 3000)
  predictions <- predict(forest, X.test,estimate.variance = TRUE)
  predictions$predictions
}


results = lapply(1:simu.reps, function(rr) single.run(n))

preds = Reduce(cbind, results)
preds.std = t(apply(preds, 1, function(xx) (xx - mean(xx)) / sd(xx)))
ymax = max(abs(preds.std))

gauss.quantiles = qnorm((1:simu.reps) / (simu.reps + 1))


all.data = data.frame(Reduce(rbind,
                             lapply(1:nrow(preds.std),
                                    function(iter) cbind(Q= gauss.quantiles, T=sort(preds.std[iter,]))
                                    ))
                      )


pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(NA, NA, xlim = range(c(gauss.quantiles, -2, 2)), ylim = c(-ymax, ymax), xlab = "theoretical quantiles", ylab = "sample quantiles")
abline(0, 1, lwd = 2, lty = 3)
boxplot(T~Q, all.data, xaxt="n", at = gauss.quantiles, boxwex=0.2, add=TRUE)
xaxes = seq(-2, 2, by = 1)
axis(1, at=xaxes, label=xaxes)
par=pardef
