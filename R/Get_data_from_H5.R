df <- read.csv('af_noac.csv', header <- T)
library(rhdf5)

estimation <- matrix(0,nrow=500)
error <- matrix(0,nrow=500)
trueeffect <- matrix(0,nrow=500)
for (i in 1:499){
  simudata <- h5read("linear_survival_data.h5",name ="train")
  #simudata <- h5read("simudata.h5",name =paste("group",i,sep=""))
  x <- t(as.matrix(simudata$x))
  t <- sapply(simudata$t,as.numeric)
  e <- sapply(simudata$e,as.numeric)
  hr <- sapply(simudata$hr,as.numeric)
  hr1 <- sapply(simudata$hr1,as.numeric)
  s_diff_true <- c(1,sapply(simudata$sp1,as.numeric))
  s_diff_true0 <- c(1,sapply(simudata$sp0,as.numeric))

  timeline1 <- sapply(simudata$t1,as.numeric)
  timeline1 <- timeline1[timeline1<100][1:min(101,length(timeline1[timeline1<100]))]
  interval1 <- c(diff(timeline1),max(0,101-max(timeline1)))
  timeline1 <- rep(timeline1,interval1)
  s_diff_true[timeline1]

  timeline0 <- sapply(simudata$t0,as.numeric)
  timeline0 <- timeline0[timeline0<100][1:min(101,length(timeline0[timeline0<100]))]
  interval0 <- c(diff(timeline0),max(0,101-max(timeline0)))
  timeline0 <- rep(timeline0,interval0)
  s_diff_true0[timeline0]

  df =as.tibble(as.matrix(cbind (t,e,x)))
  names(df)[1]='T.tilde'
  names(df)[2]='delta'
  names(df)[6]='A'
  names(df) <- sub("V","W",names(df))
  df$T.tilde <- round(df$T.tilde*100,0)
  df$ID <- as.numeric(rownames(df))
  Wname <- grep('W', colnames(df), value = TRUE)
  train <- as_tibble(df[,c('ID', Wname, 'A', "T.tilde", "delta")])
  checkbin <- function(x){
    vSet <- unique(x)
    return(length(vSet)==2)
  }
  for (it in names(train)){
    if(!checkbin(train[,it]) & it!="ID" & it!="T.tilde"){
      train[,it] <- ifelse(train[,it]<=0, 0, 1)
    }
  }
  train <- data.frame(train)
}
cat("\n Final estimated effect:", mean(estimation), mean(error))
#cat("\n Final simulated effect:", mean(trueeffect))




get_data <- function(one_sample,size,nsim,curve,curve_c,mode=NULL,intervaltype='type1'){
one_sample$ID <- seq(1,nrow(one_sample),1)

if (intervaltype == 'type1'){
interval <- c(diff(curve$T.uniq),1)
curve_individual_1 <- apply(curve$Psi.hat_1, 1, function(x) c(1, rep(x,interval)[1:curve$T.max])) %>% t()
curve_individual_0 <- apply(curve_individual_1, 1, function(x) x +pexp(seq(0,curve$T.max-1,1)/100, 2)/10%>% t()) %>% t()
#curve_individual_0 <- apply(curve$Psi.hat_0, 1, function(x) c(1, rep(x,interval)[1:curve$T.max])) %>% t()

if(mode == "true_censor") curve_c <- curve
interval <- c(diff(curve_c$T.uniq),1)
curve_individual_c_1 <- apply(curve_c$Psi.hat_1, 1, function(x) c(1, rep(x,interval)[1:curve_c$T.max])) %>% t()
curve_individual_c_0 <- apply(curve_c$Psi.hat_0, 1, function(x) c(1, rep(x,interval)[1:curve_c$T.max])) %>% t()
}else{
  curve_individual_1 <- curve$Psi.hat_1
  curve_individual_0 <- curve$Psi.hat_0
  curve_individual_c_1 <- curve_c$Psi.hat_1
  curve_individual_c_0 <- curve_c$Psi.hat_0
}


n <- nrow(curve_individual_0)
ids <- tnew <- ynew <- data.frame(matrix(nrow = size, ncol = nsim))
data_out <- list()

for(sim in 1:nsim) {
  idxs <- sample(n, size, replace = TRUE)
  ids[,sim] <- one_sample$ID[idxs]
  # event time
  u <- runif(size, 0, 1)
  # the first time survival drops below u
  stime <-  ifelse(one_sample$A[idxs]==1,
                   apply(curve_individual_1[idxs,] < u, 1, function(x) which(x)[1]),
                   apply(curve_individual_0[idxs,] < u, 1, function(x) which(x)[1]))
  w <- ifelse(one_sample$A[idxs]==1,
              curve_individual_1[idxs,length(curve$T.uniq)] > u,
              curve_individual_0[idxs,length(curve$T.uniq)] > u)
  stime <- curve$T.uniq[stime]
  stime[w] <- max(curve$T.uniq) + 1

  # censoring time
  u <- runif(size, 0, 1)
  ctime <-  ifelse(one_sample$A[idxs]==1,
                   apply(curve_individual_c_1[idxs,] < u, 1, function(x) which(x)[1]),
                   apply(curve_individual_c_0[idxs,] < u, 1, function(x) which(x)[1]))
  w <- ifelse(one_sample$A[idxs]==1,
              curve_individual_c_1[idxs,length(curve$T.uniq)] > u,
              curve_individual_c_0[idxs,length(curve$T.uniq)] > u)
  ctime <- curve$T.uniq[ctime]
  ctime[w] <- max(curve_c$T.uniq)
  if(mode == "true_censor") ctime <-ifelse(one_sample$delta[idxs]==0, one_sample$T.tilde[idxs],curve_c$T.max )

  # put it together
  tnew[,sim] <- pmin(stime, ctime)
  names(tnew) <- paste("T.tilde", 1:nsim, sep = "")
  ynew[,sim] <- stime == tnew[,sim]
  names(ynew) <- paste("delta", 1:nsim, sep = "")
  data_out[[sim]] <- data.frame(one_sample[ids[,sim],!names(one_sample) %in% c("T.tilde","delta")],
                                delta =as.numeric(ynew[,sim]), T.tilde = tnew[,sim])
  data_out[[sim]]$ID <- seq(1:nrow(data_out[[sim]]))
}

to.return <- list(data_out = data_out,
                  deltas = ynew,
                  T.tildes = tnew,
                  ids =ids
)
return(to.return)
}


