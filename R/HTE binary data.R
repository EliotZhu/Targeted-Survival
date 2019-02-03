#Simple binary data
rm(list = ls())
n = 10000
ntree = 3000
n.test = 5000

sigma = .1
d = 10

#generated data, please substitue X,W and Y with data your data
effect = function(x) {
  4/((1 + exp(-12 * (x[1] - 0.5))) * (1 + exp(-12 * (x[2] - 0.5))))
}
baseline = function(x) { 0 }
X = as.data.frame(matrix(runif(n * d, 0, 1), n, d)) # features
W = rbinom(n, 1, 1/(1+exp(-X[,1]))) #treatment condition
Z = apply(X, 1, baseline) +  (W - 0.5) * apply(X, 1, effect) + sigma * rnorm(n)
pr = 1/(1+exp(-Z)) #inv-logit function
Y = rbinom(n, 1, pr)

#fit the causal forest using test X
Z1 = apply(X, 1, baseline) + 0.5 * apply(X, 1, effect)
Z0 = apply(X, 1, baseline) - 0.5 * apply(X, 1, effect)
Y1 = rbinom(n, 1, 1/(1+exp(-Z1)))
Y0 = rbinom(n, 1, 1/(1+exp(-Z0)))
true.eff <- 1/(1+exp(-Z1))-1/(1+exp(-Z0))
true.ate <- mean(Y1-Y0)


pihat <- SuperLearner(W, X=X, SL.library=c("SL.glm", "SL.gam"), family="binomial")
pihat <- pihat$SL.predict

#causal forest
Y.fit <- SuperLearner(Y=Y, X=X, SL.library=c("SL.glm", "SL.gam","SL.rpart"), family="binomial")
Y.hat <- Y.fit$SL.predict
forest = grf::causal_forest(X, Y, W,W.hat = pihat,Y.hat=Y.hat,  num.trees = 3000)
tauhat.rf <- predict(forest, X)
tauhat.rf.adj <- tauhat.rf$predictions
hist(tauhat.rf.adj, col=rgb(0,0,1,1/4))
hist(true.eff,add=T, col=rgb(1,0,0,1/4) )

grf::average_treatment_effect(forest,target.sample="all",method="TMLE")


#bcf
bcf_fit = bcf::bcf(Y, W,as.matrix(X), as.matrix(X), pihat, nburn=2000, nsim=2000)
tau.post <- bcf_fit$tau
tauhat.rf.adj <- colMeans(tau.post)


#Plot the heat map
#set true.eff = tauhat.rf.adj if the true effect is unkonw.
print(mean((true.eff - tauhat.rf.adj)^2))
minp = min(true.eff, tauhat.rf.adj)
maxp = max(true.eff, tauhat.rf.adj)
rngp = maxp - minp

ncol=100
true.scl = pmax(ceiling(ncol * (true.eff - minp) / rngp), 1)
rf.scl = pmax(ceiling(ncol * (tauhat.rf.adj - minp) / rngp), 1) #Tmle adjusted
rf.scl.2 = pmax(ceiling(ncol * (tauhat.rf$predictions - minp) / rngp), 1) #unadjusted

hc = heat.colors(ncol)

plot(X[,1], X[,2], pch = 16, col = hc[true.scl], xlab = "", ylab = "",main = "TRUE")
plot(X[,1], X[,2], pch = 16, col = hc[rf.scl], xlab = "", ylab = "", main = "EST")
#plot(X.test[,1], X.test[,2], pch = 16, col = hc[rf.scl.2], xlab = "", ylab = "",main = "Unadjusted")














Y.fit.1 <- SuperLearner(Y=Y[W==1], X=X[W==1,], SL.library=c("SL.glm", "SL.gam","SL.rpart"), family="binomial")
Y.fit.0 <- SuperLearner(Y=Y[W==0], X=X[W==0,], SL.library=c("SL.glm", "SL.gam","SL.rpart"), family="binomial")
Y.hat.1 <- predict(Y.fit.1,as.data.frame(X))$pred
Y.hat.0 <- predict(Y.fit.0,as.data.frame(X))$pred
Y.hat <- Y.hat.1
Y.hat[W==1] <- Y.hat.1[W==1]
Y.hat[W==0] <- Y.hat.0[W==0]


#TMLE adjustment
tau.hat.pointwise <- tauhat.rf.adj
Y.hat.0 <- forest$Y.hat - forest$W.hat * tau.hat.pointwise
Y.hat.1 <- forest$Y.hat + (1 - forest$W.hat) * tau.hat.pointwise
eps.tmle.robust.0 <- lm(B ~ A + 0, data = data.frame(A = 1/(1 - forest$W.hat[forest$W.orig == 0]),
                                                     B = forest$Y.orig[forest$W.orig == 0] - Y.hat.0[forest$W.orig == 0]))
eps.tmle.robust.1 <- lm(B ~ A + 0, data = data.frame(A = 1/forest$W.hat[forest$W.orig == 1],
                                                     B = forest$Y.orig[forest$W.orig == 1] - Y.hat.1[forest$W.orig == 1]))
delta.tmle.robust.0 <- predict(eps.tmle.robust.0,
                               newdata = data.frame(A = (1/(1 - forest$W.hat))))
delta.tmle.robust.1 <- predict(eps.tmle.robust.1,
                               newdata = data.frame(A = (1/forest$W.hat)))
dr.correction <- delta.tmle.robust.1 - delta.tmle.robust.0

sigma2.hat <- sandwich::vcovHC(eps.tmle.robust.0) * mean(1/(1 - forest$W.hat))^2 +
  sandwich::vcovHC(eps.tmle.robust.1) * mean(1/forest$W.hat)^2
tauhat.rf.adj <- tauhat.rf$predictions + dr.correction
mean(tauhat.rf.adj)
tau.se <- sqrt(sigma2.hat)
