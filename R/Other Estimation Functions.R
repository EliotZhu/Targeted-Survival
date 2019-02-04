options(repr.plot.width = 5, repr.plot.height = 4)  ## resizing plots

library(np)
library(extrafont)
library(ggthemes,ggplot2)
library(survtmle2)
library(SuperLearner)

finite.differences <- function(x, y) {
  if (length(x) != length(y)) {
    stop('x and y vectors must have equal length')
  }
  
  n <- length(x)
  
  # Initialize a vector of length n to enter the derivative approximations
  fdx <- vector(length = n)
  
  # Iterate through the values using the forward differencing method
  for (i in 2:n) {
    fdx[i-1] <- (y[i-1] - y[i]) / (x[i-1] - x[i])
  }
  
  # For the last value, since we are unable to perform the forward differencing method
  # as only the first n values are known, we use the backward differencing approach
  # instead. Note this will essentially give the same value as the last iteration
  # in the forward differencing method, but it is used as an approximation as we
  # don't have any more information
  fdx[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])
  
  return(fdx)
}
initial_ <- function(data_out,end_time,s_diff_true,s_diff_true0){
  # Plot KM estimators and construct the difference
  #data_out$delta <- ifelse(data_out$T.tilde>end_time,0,data_out$delta)
  #data_out$T.tilde <- ifelse(data_out$T.tilde>end_time,end_time,data_out$T.tilde)
  full.pre1 <- survfit(Surv(T.tilde,delta)~A, data = data_out[data_out$A==1,])
  full.pre0 <- survfit(Surv(T.tilde,delta)~A, data = data_out[data_out$A==0,])
  
  interval1 <- c(diff(full.pre1$time),1)
  interval0 <- c(diff(full.pre0$time),1)
  s_diff_simu <- data.frame(treatment = NA, difference = NA,  time = c(0:(end_time-1)), group = "Simulated(KM)", mse = NA, mse2 =NA)
  s_diff_simu$treatment <- c(rep(full.pre1$surv,interval1)[1:end_time])
  s_diff_simu$difference <- c(rep(full.pre1$surv,interval1)[1:end_time]) - c(rep(full.pre0$surv,interval0)[1:end_time])
  s_diff_simu$mse <- (s_diff_simu$treatment-s_diff_true)^2
  s_diff_simu$mse2 <- (s_diff_simu$difference-(s_diff_true-s_diff_true0))^2
  
  #Inverse probability
  # nm <- plogis(predict(glm(A ~ 1, family="binomial"(link="logit"), data = data_out)))
  # Ws <- paste(grep("W",names(data_out),value = T),collapse = "+")
  # dn <- plogis(predict(glm(eval(paste("A~",Ws,sep="")), family="binomial"(link="logit"),data = data_out)))
  # data_out$weights <- ifelse(data_out$A==1, nm/dn, (1-nm)/(1-dn))
  #
  
  # quint = c(0,end_time,as.numeric(quantile(data_out[data_out$delta==0,]$T.tilde,probs=c(0.20, 0.40, 0.60, 0.80))))
  # data_out$start = -1
  # data_out$drop = ifelse(data_out$delta==0, 1, 0)
  # df.split = survSplit(data=data_out,
  #                      cut=as.numeric(quint),
  #                      end="T.tilde",
  #                      start="start",
  #                      event="delta")
  # df.split$start = 0
  # colnames(df.split)[which(colnames(df.split) %in% c("start","T.tilde"))] = c("in.t", "out.t")
  # library(data.table)
  # df.split.dt = data.table(df.split)
  # setkey(df.split.dt, "ID", "in.t")
  # first <- with(df.split.dt, c(TRUE, diff(ID) !=0)) #first id for each subject
  # last <- c(first[-1], TRUE) #last id
  # df.split.dt$drop = ifelse(last & df.split.dt$drop==1, 1, 0)
  # df.split.dt[, `:=` (nm2 = 1 - plogis(predict(glm(drop ~ 1, family="binomial", data = df.split.dt))),
  #                     dn2 = 1 - plogis(predict(glm(eval(paste("drop~",Ws,sep="")), family="binomial",data = df.split.dt))))]
  # df.split.dt[, num.do:= cumprod(nm2), by=list(ID)][, den.do:=cumprod(dn2), by=list(ID)]
  # first <- with(df.split.dt, c(TRUE, diff(ID) !=0)) #first id for each subject
  # last <- c(first[-1], TRUE) #last id
  # df.split.dt = within(df.split.dt, {
  #   w2 = ifelse(drop, (1-num.do)/(1-den.do), num.do/den.do)
  #   w3 = weights*w2
  # })
  # var.list.2 = c("ID", "in.t", "out.t", "delta", "drop", "num.do", "den.do", "w2", "w3")
  # df.split.dt[1:10, var.list.2, with=F]
  # df.split.dt=as.data.frame(df.split.dt)
  
  
  # full.pre1_ipw <- survfit(Surv(out.t,delta)~1, weights = weights, data = df.split.dt[df.split.dt$A==1,])
  # full.pre0_ipw <- survfit(Surv(out.t,delta)~1, weights = weights, data = df.split.dt[df.split.dt$A==0,])
  # interval1_ipw <- c(diff(full.pre1_ipw$time),0)
  # interval0_ipw <- c(diff(full.pre0_ipw$time),0)
  # s_diff_simu_ipw <- data.frame(treatment = NA, difference = NA,  time = c(0:end_time), group = "Simulated(KM-IPW)", mse = NA, mse2 =NA)
  # s_diff_simu_ipw$treatment <- c(1, rep(full.pre1_ipw$surv,interval1_ipw)[1:end_time])
  # s_diff_simu_ipw$difference <- c(1, rep(full.pre1_ipw$surv,interval1_ipw)[1:end_time]) -
  #   c(1, rep(full.pre0_ipw$surv,interval0_ipw)[1:end_time])
  # s_diff_simu_ipw$mse <- (s_diff_simu_ipw$treatment-s_diff_true)^2
  # s_diff_simu_ipw$mse2 <- (s_diff_simu_ipw$difference-(s_diff_true-s_diff_true0))^2
  
  return(list(s_diff_simu))
}
RSF_ <- function(data_out,end_time,s_diff_true,s_diff_true0){
  
  true_diff <- s_diff_true-s_diff_true0
  data_out <- data_out[data_out$T.tilde<=end_time,]
  W_names <- grep('W',names(data_out),value = T)
  
  
  fit.rsf <- rfsrc(as.formula( paste("Surv(T.tilde,delta)~", paste(c(W_names,'A'),collapse = "+")) ),data_out ,ntree = 1000)
  
  gHatSL_1 <- SuperLearner(Y=data_out$A, X=data_out[,W_names], SL.library=c("SL.rpart","SL.glmnet"), family="binomial")
  pscore <-  gHatSL_1$SL.predict
  require(MatchIt)
  match <- matchit(as.formula( paste("A~", paste(W_names,collapse = "+")) ), data = data_out, method = "subclass", distance = pscore)
  matched.data <- match.data(match)
  tx.indices <- (matched.data$A==1)
  control.indices <- (matched.data$A==0)
  
  
  
  data1 <- data.frame(matched.data[,!names(matched.data)=="A"],A=1)
  fit.rsf_1 <- predict(fit.rsf,data1)
  data0 <- data.frame(matched.data[,!names(matched.data)=="A"],A=0)
  fit.rsf_0 <- predict(fit.rsf,data0)
  
  interval1 <- c(diff(fit.rsf_1$time.interest),1)
  if(sum(interval1)<end_time) interval1[length(interval1)] <- end_time-sum(interval1)+1
  interval0 <- c(diff(fit.rsf_0$time.interest),1)
  if(sum(interval0)<end_time) interval0[length(interval0)] <- end_time-sum(interval0)+1
  
  Qn.A1.t_1 <- fit.rsf_1$survival
  Qn.A1.t_full_1 <- Qn.A1.t_1[,c(rep(seq(1,ncol(Qn.A1.t_1)),interval1))]
  Qn.A1.t_0 <- fit.rsf_0$survival
  Qn.A1.t_full_0 <- Qn.A1.t_0[,c(rep(seq(1,ncol(Qn.A1.t_0)),interval0))]
  
  s_diff_rsf <- data.frame(treatment = NA,difference = NA, time = seq(0,end_time-1), group = "RSF", mse =NA, mse2 =NA, varest = NA)
  s_diff_rsf$treatment <- colMeans(Qn.A1.t_full_1)
  psi <- sweep((Qn.A1.t_full_1 - Qn.A1.t_full_0),1, matched.data$weights,"*")
  s_diff_rsf$difference <- colSums(psi/sum(matched.data$weights))
  s_diff_rsf$varest <- colSums(((psi-s_diff_rsf$difference)^2*matched.data$weights)/(sum(matched.data$weights) * sum(matched.data$weights)))
  
  s_diff_rsf$mse<- (s_diff_rsf$treatment-s_diff_true)^2
  s_diff_rsf$mse2 <- (s_diff_rsf$difference-true_diff)^2
  
  s_diff_rsf_his <- s_diff_rsf
  
  return(list(s_diff_rsf_his))
}
RSFc_ <- function(data_out,end_time,s_diff_true,s_diff_true0,strata = 'W1'){
  
  true_diff <- s_diff_true-s_diff_true0
  data_out <- data_out[data_out$T.tilde<=end_time,]
  W_names <- grep('W',names(data_out),value = T)
  
  
  fit.rsf <- rfsrc(as.formula( paste("Surv(T.tilde,delta)~", paste(c(W_names,'A'),collapse = "+")) ),data_out ,ntree = 1000,forest.wt=T)
  
  gHatSL_1 <- SuperLearner(Y=data_out$A, X=data_out[,W_names], SL.library=c("SL.rpart","SL.glmnet"), family="binomial")
  pscore <-  gHatSL_1$SL.predict
  require(MatchIt)
  match <- matchit(as.formula( paste("A~", paste(W_names,collapse = "+")) ), data = data_out, method = "subclass", distance = pscore)
  matched.data <- match.data(match)
  tx.indices <- (matched.data$A==1)
  control.indices <- (matched.data$A==0)
  
  assign(paste(strata),5)
  rhs <- paste(strata, 1:10, "<-", lhs, sep="")
  eval(parse(text=rhs))
  data1 <- data.frame(matched.data[,c(W_names[!W_names %in% strata],'delta','T.tilde')],A=1, ),
  fit.rsf_1 <- predict(fit.rsf,data1,forest.wt = T)
  data0 <- data.frame(matched.data[,W_names],A=0)
  fit.rsf_0 <- predict(fit.rsf,data0)
  
  interval1 <- c(diff(fit.rsf_1$time.interest),1)
  if(sum(interval1)<end_time) interval1[length(interval1)] <- end_time-sum(interval1)+1
  interval0 <- c(diff(fit.rsf_0$time.interest),1)
  if(sum(interval0)<end_time) interval0[length(interval0)] <- end_time-sum(interval0)+1
  
  Qn.A1.t_1 <- fit.rsf_1$survival
  Qn.A1.t_full_1 <- Qn.A1.t_1[,c(rep(seq(1,ncol(Qn.A1.t_1)),interval1))]
  Qn.A1.t_0 <- fit.rsf_0$survival
  Qn.A1.t_full_0 <- Qn.A1.t_0[,c(rep(seq(1,ncol(Qn.A1.t_0)),interval0))]
  
  s_diff_rsf <- data.frame(treatment = NA,difference = NA, time = seq(0,end_time-1), group = "RSF", mse =NA, mse2 =NA, varest = NA)
  s_diff_rsf$treatment <- colMeans(Qn.A1.t_full_1)
  psi <- sweep((Qn.A1.t_full_1 - Qn.A1.t_full_0),1, matched.data$weights,"*")
  s_diff_rsf$difference <- colSums(psi/sum(matched.data$weights))
  s_diff_rsf$varest <- colSums(((psi-s_diff_rsf$difference)^2*matched.data$weights)/(sum(matched.data$weights) * sum(matched.data$weights)))
  
  s_diff_rsf$mse<- (s_diff_rsf$treatment-s_diff_true)^2
  s_diff_rsf$mse2 <- (s_diff_rsf$difference-true_diff)^2
  
  s_diff_rsf_his <- s_diff_rsf
  
  return(list(s_diff_rsf_his))
}
BNN_ <- function(data_out,end_time,s_diff_true,s_diff_true0){
  
  true_diff <- s_diff_true-s_diff_true0
  wnames <- colnames(data_out[,!names(data_out)%in%c("ID")])
  
  bnn_model <- bnnSurvival(Surv(T.tilde,delta)~., data = data_out[,wnames],k = 10, num_base_learners = length(wnames)-2, num_features_per_base_learner = length(wnames)-2)
  fit.bnn_1 <- predict(bnn_model, data_out[data_out$A==1,!names(data_out)%in%"ID"])
  fit.bnn_0 <- predict(bnn_model,data_out[data_out$A==0,!names(data_out)%in%"ID"])
  interval1 <- c(diff(timepoints(fit.bnn_1)),1)
  interval0 <- c(diff(timepoints(fit.bnn_0)),1)
  
  s_diff_bnn <- data.frame(treatment = NA,difference = NA, time = seq(0,end_time-1), group = "BNN", mse =NA, mse2 =NA )
  
  s_diff_bnn$treatment <- c(rep(colMeans(predictions(fit.bnn_1)),interval1))[1:end_time]
  s_diff_bnn$difference <-s_diff_bnn$treatment - c(rep(colMeans(predictions(fit.bnn_0)),interval0))[1:end_time]
  s_diff_bnn$mse<- (s_diff_bnn$treatment-s_diff_true[1:end_time])^2
  s_diff_bnn$mse2 <- (s_diff_bnn$difference-true_diff[1:end_time])^2
  
  s_diff_bnn_his <- s_diff_bnn
  
  return(list(s_diff_bnn_his))
}
COX_ <- function(data_out,end_time,s_diff_true,s_diff_true0){
  y <- Surv(data_out$T.tilde,data_out$delta)
  Wnames <- grep("W",names(data_out),value=T)
  full.ph <-  coxph(y~. ,  data =data_out[,c(Wnames,'A')], method="breslow")
  coxfit <- survfit(full.ph)
  #mean(as.matrix(data_out[data_out$A==0,c(Wnames,'A')]) %*%  as.matrix(full.ph$coefficients))
  survConcordance(full.ph$y~predict(full.ph))
  
  bh=basehaz(full.ph)
  x_1 <- as.matrix(data.frame(data_out[,c(Wnames)],A=1)) %*% full.ph$coef
  x_0 <- as.matrix(data.frame(data_out[,c(Wnames)],A=0)) %*% full.ph$coef
  s_diff_cox1 <- sapply(exp(-bh[,1]), function(s) s*(exp(x_1)))
  s_diff_cox0 <- sapply(exp(-bh[,1]), function(s) s*(exp(x_0)))
  
  interval <- c(bh[-1,2] - bh[-length(bh[,2]),2],max(end_time-bh[length(bh[,2]),2],1))
  
  
  s_diff_cox <- data.frame(treatment = NA,difference = NA, time = seq(0,end_time-1), group = "COX", mse =NA, mse2 =NA )
  s_diff_cox$treatment <- c(1, rep(colMeans(s_diff_cox1),interval)[1:end_time-1])
  s_diff_cox$difference <- c(1, rep(colMeans(s_diff_cox1)-colMeans(s_diff_cox0),interval)[1:end_time-1])
  
  s_diff_cox$mse <- (s_diff_cox$treatment-s_diff_true)^2
  s_diff_cox$mse2 <- (s_diff_cox$difference-(s_diff_true-s_diff_true0))^2
  
  s_diff_cox_his <- s_diff_cox
  return(list(s_diff_cox_his))
}
survplot_0 <- function(result,s_diff_true,s_diff_true0,end_time,size_eval=1){
  s_diff <- data.frame(treatment = s_diff_true[1:end_time],difference = s_diff_true[1:end_time]-s_diff_true0[1:end_time],
                       time = seq(0,end_time-1), group = "True")
  s_diff$mse <- 0
  s_diff$mse2 <- 0
  s_diff$SE <- 0
  s_diff$varest <- 0
  s_diff_his <- s_diff
  
  result <-  TMLE_result
  
  s_diff_tmle <- list()
  for (i in seq(1,size_eval*2,2)){
    s_diff_tmle <- append(s_diff_tmle,list(TMLE_result[[i]]))
  }
  if(length(result)>1){
    s_diff_SL <- list()
    for (i in seq(2,size_eval*2,2)){
      s_diff_SL <- append(s_diff_SL,list(TMLE_result[[i]]))
    }
    s_diff_SL_his <-  Reduce("+", lapply(s_diff_SL, function(x) replace(x, is.na(x), NA)))
    s_diff_SL_his$time <- seq(0,end_time-1)
    s_diff_SL_his$group <- "Second Method"
    s_diff_SL_his[,c('mse','mse2','treatment','difference','varest')] <- s_diff_SL_his[,c('mse','mse2','treatment','difference','varest')]/size_eval
    s_diff_SL_his$SE <- plyr::ldply(s_diff_SL, function(x) x$difference) %>% sapply(sd) %>% as.vector()/sqrt(size_eval)
    s_diff_his <- s_diff_his%>%rbind(s_diff_SL_his)
  }
  
  s_diff_tmle_his <-  Reduce("+", lapply(s_diff_tmle, function(x) replace(x, is.na(x), NA)))
  # s_diff_tmle_his <-  Reduce("+", lapply(map(TMLE_result, 1), function(x) replace(x, is.na(x), NA)))
  # s_diff_SL_his <-  Reduce("+", lapply(map(TMLE_result, 2), function(x) replace(x, is.na(x), NA)))
  
  s_diff_tmle_his$time <- seq(0,end_time-1)
  s_diff_tmle_his$group <- "First Method (TMLE)"
  s_diff_tmle_his[,c('mse','mse2','treatment','difference','varest')] <- s_diff_tmle_his[,c('mse','mse2','treatment','difference','varest')]/size_eval
  s_diff_tmle_his$SE <- plyr::ldply(s_diff_tmle, function(x) x$difference) %>% sapply(sd) %>% as.vector()/sqrt(size_eval)
  s_diff_his <- s_diff_his%>%rbind(s_diff_tmle_his)
  
  
  g <-  ggplot(data=s_diff_his, aes(x=time, y=s_diff_his$difference, colour = group, lty = group )) +
    geom_line(size=0.5)+
    labs(x = "Time",y = "Difference in survival probablity")+
    scale_color_manual(values=c("#85B8CB","#5A0651","#283B42","#E81E25","#2A3132","#FE7773"))+
    scale_linetype_manual(values=c("solid", "dotted", "twodash","solid","solid","solid"))+
    theme_hc()+
    theme(legend.position="right")
  return(g)
  
}


#write_csv(dat[[1]],path='isimuwei.csv')

#exp(-seq(0,10,.1)/exp(0.767864))
#s_diff_his[,c(1,2)] <- 0


#Simulate and check the data.
endtime <- 20
dat <-  get.data(4321,5e3,"moderate",tpts=endtime)
table(dat$datt_out$T.tilde)


curve <- surv_curve (data_out,dW = rep(1, nrow(data_out)),T.cutoff = 10,method =  'SL')
colMeans(curve$Psi.hat_1)-colMeans(curve$Psi.hat_0)


size =2000
nsim = 1
simulation <- get_data(dat,size,nsim,curve,curve_c,mode = "true_censor",intervaltype='type1')
g = ggplot()+theme_hc()+
  xlim(0, 100)+labs(x = "Survival Time",y = "Density")
for (i in 1:nsim){
  temp <- data.frame(x=simulation$T.tildes[simulation$deltas[,i],i])
  g <-  g + geom_line(data = temp,aes(x, colour="Simulated Time"), alpha=0.3, stat="density", show.legend = TRUE)
}
g <- g + geom_line(data = dat,aes(T.tilde, colour="Original Time"), stat="density", show.legend = TRUE)
g <- g+ scale_colour_manual( values= c( "#F62A00","#828282"),
                             labels = c("Originael Data", "Simulated Data (50 sets)"))+
  theme(legend.title=element_blank())+
  theme(legend.position=c(.8, .9))

dat_eval <-  simulation$data_out

plot(survfit(Surv(T.tilde,delta)~A, data = dat_eval[[1]]))


end_time <- curve$T.max
interval <- c(diff(curve$T.uniq),1)
curve_1 <- c(rep(colMeans(curve$Psi.hat_1),interval)[1:curve$T.max])
curve_0 <- c(rep(colMeans(curve$Psi.hat_0),interval)[1:curve$T.max])
s_diff_true <-  curve_1
s_diff_true0 <-  curve_1 + pexp(seq(0,curve$T.max-1,1)/100, 2)/10





#evluation
#s_diff <- data.frame(treatment = s_diff_true[1:101],difference = s_diff_true[1:101]-s_diff_true0[1:101], time = c(0:end_time), group = "True")
s_diff <- data.frame(treatment = s_diff_true[1:end_time],difference = s_diff_true[1:end_time]-s_diff_true0[1:end_time],
                     time = seq(0,end_time-1), group = "True")
s_diff$mse <- 0
s_diff$mse2 <- 0
s_diff$SE <- 0
s_diff$varest <- 0
s_diff_his <- s_diff

size_eval <-  nsim
dat_eval <-  simulation$data_out
pb <- txtProgressBar(min = 0, max = size_eval, style = 3)

#s_diff_his <- s_diff_his%>%rbind(s_diff_simu_his_ipw)

RSF_result <- list()
for (i in 1:size_eval){
  Sys.sleep(0.1)
  # update progress bar
  setTxtProgressBar(pb, i)
  RSF_result <-  append(RSF_result, RSF_(dat_eval[[i]],end_time,s_diff_true,s_diff_true0))
}
close(pb)
s_diff_rsf_his <-  Reduce("+", lapply(RSF_result, function(x) replace(x, is.na(x), NA)))
s_diff_rsf_his$time <- seq(0,end_time-1)
s_diff_rsf_his$group <- "RSF"
s_diff_rsf_his[,c('mse','treatment','mse2','difference','varest')] <- s_diff_rsf_his[,c('mse','treatment','mse2','difference','varest')]/size_eval
s_diff_rsf_his$SE <- plyr::ldply(RSF_result, function(x) x$difference) %>% sapply(sd) %>% as.vector()/sqrt(size_eval)
s_diff_his <- s_diff_his%>%rbind(s_diff_rsf_his)
s_diff_his$difference <- ifelse(s_diff_his$difference==0,NA,s_diff_his$difference)


# SL+TMLE fitting
#TMLE_result <- parLapply(cl, dat, function(x) TMLE_(data_out=x,end_time,s_diff_true,s_diff_true0))

TMLE_result <- list()
for (i in 1:size_eval){
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
  TMLE_result <- append(TMLE_result,TMLE_(dat_eval[[i]],end_time,s_diff_true,s_diff_true0))
}

s_diff_tmle <- list()
for (i in seq(1,size_eval*2,2)){
  s_diff_tmle <- append(s_diff_tmle,list(TMLE_result[[i]]))
}
s_diff_SL <- list()
for (i in seq(2,size_eval*2,2)){
  s_diff_SL <- append(s_diff_SL,list(TMLE_result[[i]]))
}

s_diff_tmle_his <-  Reduce("+", lapply(s_diff_tmle, function(x) replace(x, is.na(x), NA)))
s_diff_SL_his <-  Reduce("+", lapply(s_diff_SL, function(x) replace(x, is.na(x), NA)))
# s_diff_tmle_his <-  Reduce("+", lapply(map(TMLE_result, 1), function(x) replace(x, is.na(x), NA)))
# s_diff_SL_his <-  Reduce("+", lapply(map(TMLE_result, 2), function(x) replace(x, is.na(x), NA)))

s_diff_tmle_his$time <- seq(0,end_time-1)
s_diff_SL_his$time <- seq(0,end_time-1)
s_diff_tmle_his$group <- "TMLE"
s_diff_SL_his$group <- "SL"
s_diff_tmle_his[,c('mse','mse2','treatment','difference','varest')] <- s_diff_tmle_his[,c('mse','mse2','treatment','difference','varest')]/size_eval
s_diff_SL_his[,c('mse','mse2','treatment','difference','varest')] <- s_diff_SL_his[,c('mse','mse2','treatment','difference','varest')]/size_eval
s_diff_tmle_his$SE <- plyr::ldply(s_diff_tmle, function(x) x$difference) %>% sapply(sd) %>% as.vector()/sqrt(size_eval)
s_diff_SL_his$SE <- plyr::ldply(s_diff_SL, function(x) x$difference) %>% sapply(sd) %>% as.vector()/sqrt(size_eval)


s_diff_his <- s_diff_his%>%rbind(s_diff_SL_his)
s_diff_his <- s_diff_his%>%rbind(s_diff_tmle_his)


ggplot(data=s_diff_his, aes(x=time, y=s_diff_his$difference, colour = group, lty = group )) +
  geom_line(size=0.5)+
  labs(x = "Time",y = "Difference in survival probablity")+
  scale_color_manual(values=c("#85B8CB","#5A0651","#283B42","#E81E25","#2A3132","#FE7773"))+
  scale_linetype_manual(values=c("solid", "dotted", "twodash","solid","solid","solid"))+
  theme_hc()+
  theme(legend.position="right")


ggsave("One-estimation (Moderate-10var-5000sim).png", plot = last_plot(), device = 'png', path = '~/planets/GFX1',
       scale = 1, width = 150, height = 100, units = c("mm"), dpi = 300)
ggsave("MSE (1000 Sample).png", plot = last_plot(), device = 'png', path = '~/planets/GFX',
       scale = 1, width = 150, height = 100, units = c("mm"), dpi = 300)
ggsave("Density (50 sets).png", plot = last_plot(), device = 'png', path = '~/planets/GFX',
       scale = 1, width = 150, height = 100, units = c("mm"), dpi = 300)
ggsave("SE (50 sets).png", plot = last_plot(), device = 'png', path = '~/planets/GFX',
       scale = 1, width = 150, height = 100, units = c("mm"), dpi = 300)







BNN_result <- list()
for (i in 1:size_eval){
  Sys.sleep(0.1)
  # update progress bar
  setTxtProgressBar(pb, i)
  BNN_result <-  append(BNN_result, BNN_(dat_eval[[i]],end_time,s_diff_true,s_diff_true0))
}
close(pb)
s_diff_bnn_his <-  Reduce("+", lapply(BNN_result, function(x) replace(x, is.na(x), NA)))
s_diff_bnn_his$time <- seq(0,end_time-1)
s_diff_bnn_his$group <- "BNN"
s_diff_bnn_his[,c('mse','treatment','mse2','difference')] <- s_diff_bnn_his[,c('mse','treatment','mse2','difference')]/size_eval
s_diff_bnn_his$SE <- plyr::ldply(BNN_result, function(x) x$difference) %>% sapply(sd) %>% as.vector()/sqrt(size_eval)
s_diff_his <- s_diff_his%>%rbind(s_diff_bnn_his)
s_diff_his$difference <- ifelse(s_diff_his$difference==0,NA,s_diff_his$difference)





#KM
#KM_result <- parLapply(cl, dat_eval, function(x) initial_(data_out = x, end_time,s_diff_true,s_diff_true0))
KM_result <- list()
for (i in 1:size_eval){
  KM_result <- append(KM_result,initial_(dat_eval[[i]], end_time,s_diff_true,s_diff_true0))
}


#s_diff_simu_his <-  Reduce("+", lapply(map(KM_result, 1), function(x) replace(x, is.na(x), NA)))
s_diff_simu_his <-  Reduce("+", lapply(KM_result, function(x) replace(x, is.na(x), NA)))
s_diff_simu_his$time <- seq(0,end_time-1)
s_diff_simu_his$group <- "KM"
#s_diff_simu_his_ipw$time <- c(0:end_time)
#s_diff_simu_his_ipw$group <- "Simulated(KM-IPW)"
s_diff_simu_his[,c('mse','treatment','mse2','difference')] <- s_diff_simu_his[,c('mse','treatment','mse2','difference')]/size_eval
s_diff_simu_his$SE <- plyr::ldply(KM_result, function(x) x$difference) %>% sapply(sd) %>% as.vector()/sqrt(size_eval)
#s_diff_simu_his_ipw[,c('mse','treatment','mse2','difference')] <- s_diff_simu_his_ipw[,c('mse','treatment','mse2','difference')]/iti
s_diff_his <- s_diff_his%>%rbind(s_diff_simu_his)
s_diff_his$difference <- ifelse(s_diff_his$difference==0,NA,s_diff_his$difference)











grey_black <- c("#2A3132","#505160","#0E0301")
green_gradient <- c("#003B46","#07575B","#66A5AD")
red_gradient <- c("#E81E25","#FE7773")
blue_gradient<- c("#85B8CB","#1D6A96","#283B42")
contrast_color <- c("#C5001A","#113743","#E4E3DB","C5BEBA")











