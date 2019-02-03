#simu data
dat <-  get.data(54321,2000,"scenario 3",endtime =12,ratDiv=500)
data_out <- dat$datt_out
data_out <- data_out[data_out$T.tilde<100 & data_out$T.tilde>0,]
data_out <- data_out[complete.cases(data_out),]
table(data_out$T.tilde)
table(data_out$Delta[data_out$T.tilde<12])/nrow(data_out)
table(data_out$Delta[data_out$T.tilde>=12])/nrow(data_out)


testdat <-  get.data(24321,1100,"moderate",endtime =12,ratDiv=70)
testdat <- testdat$datt_out
testdat <- testdat[testdat$T.tilde<100 & testdat$T.tilde>0,]
testdat <- testdat[complete.cases(testdat),]
testdat <- testdat[1:1000,]
# ---------------------------------------------------------------------------------------
# KM
library(survival)
n.data <- nrow(data_out)
km.fit <- survfit(Surv(time = T.tilde, event = Delta) ~ A, data = data_out[data_out$T.tilde<=24,])
ggsurvplot(km.fit, risk.table = TRUE, pval = TRUE, conf.int = TRUE, xlim = c(0,24), ylim = c(0,1),
           break.time.by =1, ggtheme =theme_hc(),risk.table.y.text.col = T, risk.table.y.text = FALSE)




# ---------------------------------------------------------------------------------------
SL.gam_.5 <- function (... , deg.gam = .5) {
  SL.gam (... , deg.gam = deg.gam)
}
SL.gam_2 <- function (... , deg.gam = 2) {
  SL.gam (... , deg.gam = deg.gam)
}
library(MOSSATE)
detach('package:mgcv',unload=TRUE)
options(mc.cores = 8)
set.seed(25431, "L'Ecuyer-CMRG")

onestepfit <- MOSSATE::MOSS$new(data_out, dW = 1,
                        verbose = TRUE, epsilon.step = 1e-10, max.iter = 1)
onestepfit$onestep_curve(g.SL.Lib = c("SL.gam","SL.glmnet","SL.mean","SL.earth"),
                           Delta.SL.Lib = c("SL.mean","SL.earth","SL.glmnet","SL.gam"),
                           ht.SL.Lib = c("SL.gam","SL.mean","SL.earth","SL.gam_.5","SL.glmnet"),
                           env = parent.frame())

onestepfit1 <- MOSS::MOSS$new(data_out, dW = 1,
                                verbose = TRUE, epsilon.step = 1e-1, max.iter = 50)
onestepfit1$onestep_curve(g.SL.Lib = c("SL.gam","SL.glmnet","SL.mean","SL.earth"),
                         Delta.SL.Lib = c("SL.mean","SL.earth","SL.glmnet","SL.gam"),
                         ht.SL.Lib = c("SL.gam","SL.mean","SL.earth","SL.gam_.5","SL.glmnet"))

onestepfit0 <- MOSS::MOSS$new(data_out, dW = 0,
                             verbose = TRUE, epsilon.step = 1e-1, max.iter = 50)
onestepfit0$onestep_curve(g.SL.Lib = c("SL.gam","SL.glmnet","SL.mean","SL.earth"),
                         Delta.SL.Lib = c("SL.mean","SL.earth","SL.glmnet","SL.gam"),
                         ht.SL.Lib = c("SL.gam","SL.mean","SL.earth","SL.gam_.5","SL.glmnet"))



onestepdiff <- MOSS_difference$new(data_out, verbose = TRUE, epsilon.step = 1e-3,
                                   max.iter = 10,
                                   ftimeMod = onestepfit$ftimeMod,
                                   ctimeMod = onestepfit$ctimeMod,
                                   trtMod = onestepfit$trtMod)
onestepdiff$initial_fit(env = parent.frame())
onestepdiff$onestep_diff_curve( )


tpts <- seq(1:12)
ptps <- endtime <- 12
A <- onestepdiff$dat$A
controls <- onestepdiff$dat[,grep("W",names(onestepdiff$dat),value = T)]

#Get the true effect by treatment
true.func.1 <- function(x,points,ratDiv=500){
  x <- as.matrix(x,nrow=1)
  rate <- as.numeric(x[5]+x[6]+x[7]+x[8]+x[1]+x[2]+x[3]+x[4])
  s_diff_true <-  exp(-rate/ratDiv*seq(0,points,1))
  return(s_diff_true)
}
true.1 <- t(apply(controls, 1, function(x) true.func.1(x,24)))

true.func.0 <- function(x,points,ratDiv=500){
  x <- as.matrix(x,nrow=1)
  rate <- as.numeric(x[5]+x[6]+x[7]+x[8])
  s_diff_true <-  exp(-rate/ratDiv*seq(0,points,1))
  return(s_diff_true)
}
true.0 <- t(apply(controls, 1, function(x) true.func.0(x,24)))

plot(km.fit, col=c('palegreen3','darkgreen'), lty = 2,main = paste('n=', n.data, '\n # of nuisance covariate = ?'),xlab = 'Time')

lines(colMeans(true.1), type="l", cex=0.2, col = 'red')
lines(onestepfit1$Psi.hat, type="l", cex=.2, col = 'red',lty=2) #This show the bias
lines(colMeans(onestepdiff$MOSS_A1$Qn.A1.t_full), type="l", cex=.2, col = 'red',lty=3) #This show the bias
lines(colMeans(onestepdiff$MOSS_A1$Qn.A1.t_initial), type="l", cex=.2, col = 'red',lty=4) #This show the bias

lines(colMeans(true.0), type="l", cex=0.2, col = 'blue')
lines(onestepfit0$Psi.hat, type="l", cex=.2, col = 'blue',lty=2) #This show the bias
lines(colMeans(onestepdiff$MOSS_A0$Qn.A1.t_full), type="l", cex=.2, col = 'blue',lty=3) #This show the bias
lines(colMeans(onestepdiff$MOSS_A0$Qn.A1.t_initial), type="l", cex=.2, col = 'blue',lty=4) #This show the bias

plot(colMeans(true.1-true.0), type="l", cex=0.2, col = 'green')
lines(onestepdiff$Psi.hat, type="l", cex=.2, col = 'blue',lty=4) #This show the bias
lines(colMeans(onestepdiff$MOSS_A1$Qn.A1.t_initial-onestepdiff$MOSS_A0$Qn.A1.t_initial), type="l", cex=.2, col = 'blue',lty=3) #This show the bias



#Plug in the estimator
Y.hat.1 <- onestepfit.1$Qn.A1.t_initial[,1:endtime]
Y.hat.0 <- onestepfit.0$Qn.A1.t_initial[,1:endtime]
prediction <- Y.hat.1-Y.hat.0
prediction_tmle <- onestepfit.1$Qn.A1.t_full[,1:endtime]-onestepfit.0$Qn.A1.t_full[,1:endtime]

#Y.hat.0.adj <- as.matrix(Y.hat.0+as.tibble(t(colMeans(prediction))) %>% dplyr::slice(rep(1:n(), each = nrow(Y.hat.0))))
#ks.test(Y.hat.1[,ptps],Y.hat.0.adj[,ptps])


#Get the real and average effect
true.eff <- (true.1[,1:endtime])-(true.0[,1:endtime])
plot(colMeans(true.eff), type="l", cex=0.2, col = 'green')
lines(colMeans(prediction), type="l", cex=0.2, col = 'blue',lty=4) #This show the bias
lines(colMeans(prediction_tmle), type="l", cex=0.2, col = 'red',lty=2) #This show the bias



out.ite <- function(filename,tau.hat.tmle,tau.hat,tau=NULL,endtime){
  point.rmse <- function(ptps,tau.hat,tau=NULL){
    if(!is.null(tau)){
      RMSE <- sqrt(mean((tau.hat[,ptps] - tau[,ptps])^2))
    }else{
      RMSE <- NA
    }
    return(RMSE)
  }
  se <- rmse <- matrix(1:endtime)
  for (i in 1:endtime){
    rmse[i] <-  point.rmse(ptps=i, tau.hat=tau.hat, tau=tau)
    se[i] <-  sd(tau.hat)
  }
  out.dat <- data.frame(rmse=rmse,se=se)
  write.csv(out.dat, file = paste0(getwd(),'/result/summary_',filename,".csv") , row.names = FALSE)
  write.csv(tau.hat, file = paste0(getwd(),'/result/raw_',filename,".csv") , row.names = FALSE)
  write.csv(tau.hat.tmle, file = paste0(getwd(),'/result/tmle_',filename,".csv") , row.names = FALSE)
}
out.ite("CH_3_By_Size_AUG_20000",tau.hat.tmle=prediction_tmle,tau.hat=prediction,tau=true.eff,12)



plot(true.eff[,ptps] , prediction[,ptps])




n <- 500
p <- 50
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X) <- paste("X", 1:p, sep="")
X <- data.frame(X)
Y <- X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + rnorm(n)
SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.polymars", "SL.mean")
test <- SuperLearner(Y = Y, X = X, SL.library = SL.library,
                     verbose = TRUE, method = "method.NNLS")
test$SL.predict
