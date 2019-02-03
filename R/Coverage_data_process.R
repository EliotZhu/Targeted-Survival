get_data <- function(one_sample,size,nsim,T.uniq,curve1,curve0,mode=NULL,idxs=NULL){
  one_sample$ID <- seq(1,nrow(one_sample),1)

  interval <- c(diff(T.uniq),1)
  #curve_individual_1 <- apply(curve1, 1, function(x) c(1, rep(x,interval)[1:max(T.uniq)])) %>% t()
  #curve_individual_0 <- apply(curve0, 1, function(x) x +pexp(seq(0,max(T.uniq)-1,1)/100, 2)/10%>% t()) %>% t()
  curve_individual_1 <- curve1
  curve_individual_0 <- curve0

  if(mode == "true_censor") curve_c <- curve
  curve_individual_c_1 <- curve1
  curve_individual_c_0 <- curve0

  n <- nrow(curve_individual_0)
  ids <- tnew <- ynew <- data.frame(matrix(nrow = size, ncol = nsim))
  data_out <- list()

  for(sim in 1:nsim) {
    if(is.null(idxs)){
      idxs <- sample(n, size, replace = TRUE)
    }
    ids[,sim] <- one_sample$ID[idxs]
    # event time
    u <- runif(size, 0, 1)
    # the first time survival drops below u
    stime <-  ifelse(one_sample$A[idxs]==1,
                     apply(curve_individual_1[idxs,] < u, 1, function(x) which(x)[1]),
                     apply(curve_individual_0[idxs,] < u, 1, function(x) which(x)[1]))
    w <- ifelse(one_sample$A[idxs]==1,
                curve_individual_1[idxs,length(T.uniq)] > u,
                curve_individual_0[idxs,length(T.uniq)] > u)
    stime <- T.uniq[stime]
    stime[w] <- max(T.uniq) + 1


    if(mode == "true_censor") ctime <-ifelse(one_sample$Delta[idxs]==0, one_sample$T.tilde[idxs],max(T.uniq))

    # put it together
    tnew[,sim] <- pmin(stime, ctime)
    names(tnew) <- paste("T.tilde", 1:nsim, sep = "")
    ynew[,sim] <- stime == tnew[,sim]
    names(ynew) <- paste("Delta", 1:nsim, sep = "")
    data_out[[sim]] <- data.frame(one_sample[ids[,sim],!names(one_sample) %in% c("T.tilde","Delta")],
                                  Delta =as.numeric(ynew[,sim]), T.tilde = tnew[,sim])
    data_out[[sim]]$ID <- seq(1:nrow(data_out[[sim]]))
    rownames(data_out[[sim]]) <- NULL

  }

  to.return <- list(data_out = data_out,
                    deltas = ynew,
                    T.tildes = tnew,
                    ids =ids
  )
  return(to.return)
}
surv_simulate <- function(data_out,size,nsim,g.SL.Lib,Delta.SL.Lib,ht.SL.Lib,idxs){
  truefit <- MOSSATE::MOSS$new(data_out, dW = 1,
                               verbose = TRUE, epsilon.step = 1e-1, max.iter = 0)
  truefit$onestep_curve(g.SL.Lib, Delta.SL.Lib, ht.SL.Lib,env = parent.frame())
  truefit.1 <- MOSSATE::MOSS$new(data_out, dW = 0, verbose = TRUE, epsilon.step = 1e-1,
                                 max.iter = 0,pred_data = T,
                                 ftimeMod = truefit$ftimeMod,
                                 ctimeMod = truefit$ctimeMod,
                                 trtMod = truefit$trtMod)
  truefit.1$onestep_curve( env = parent.frame())
  truefit.0 <- MOSSATE::MOSS$new(data_out, dW = 1, verbose = TRUE, epsilon.step = 1e-1,
                                 max.iter = 0,pred_data = T,
                                 ftimeMod = truefit$ftimeMod,
                                 ctimeMod = truefit$ctimeMod,
                                 trtMod = truefit$trtMod)
  truefit.0$onestep_curve( env = parent.frame())
  curve0 <- truefit.0$Qn.A1.t_initial
  curve1 <- truefit.1$Qn.A1.t_initial
  T.uniq <- truefit.1$T.uniq
  simulation <- get_data(data_out,size,nsim,T.uniq,curve1,curve0,mode="true_censor",idxs=idxs)
  return(simulation)
}


SL.gam_.5 <- function (... , deg.gam = .5) {
  SL.gam (... , deg.gam = deg.gam)
}
SL.gam_2 <- function (... , deg.gam = 2) {
  SL.gam (... , deg.gam = deg.gam)
}


#simu data
dat <-  get.data(54321,1000,"scenario 3",endtime =12,ratDiv=115)
data_out <- dat$datt_out
data_out <- data_out[data_out$T.tilde<100 & data_out$T.tilde>0,]
data_out <- data_out[complete.cases(data_out),]
table(data_out$T.tilde)
table(data_out$Delta[data_out$T.tilde<12])/nrow(data_out)
table(data_out$Delta[data_out$T.tilde>=12])/nrow(data_out)

a <- surv_simulate(data_out,size = 800,nsim = 100,
              g.SL.Lib = c("SL.gam","SL.glmnet"),
              Delta.SL.Lib = c("SL.glmnet","SL.gam"),
              ht.SL.Lib = c("SL.glmnet","SL.gam"),
              idxs = c(1:800))
simulation <-a
simulation$data_out
simulation$T.tildes
g = ggplot()+theme_hc()+
  xlim(0, 100)+labs(x = "Survival Time",y = "Density")
for (i in 1:100){
  temp <- data.frame(x=simulation$T.tildes[simulation$deltas[,i],i])
  g <-  g + geom_line(data = temp,aes(x, colour="Simulated Time"), alpha=0.3, stat="density", show.legend = TRUE)
}
g <- g + geom_line(data = data_out,aes(T.tilde, colour="Original Time"), stat="density", show.legend = TRUE)
g <- g+ scale_colour_manual( values= c( "#F62A00","#828282"),
                             labels = c("Original Data", "Simulated Data (50 sets)"))+
  theme(legend.title=element_blank())+
  theme(legend.position=c(.8, .9))
g


