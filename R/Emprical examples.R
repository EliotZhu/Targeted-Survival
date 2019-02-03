
data_out <- Blanca[,c('female',"age","imd2015","prev_vka","af_hospital","mi_count",
                      "chf_count","pvd_count","cvd_count","dem_count","copd_count","rhe_count","dtm_u_count","dtm_c_count",
                      "ckd_count","can_count","charlson_index_2y","vas_count","hyp_count",
                      "chadvasc_index_2y","bleed_count","alch_count","nsaid_asp_count","l_inr_count","hasbled_index_2y")]
data_out[is.na(data_out)] <- 0
data_out <- as.data.frame(apply(data_out, 2, as.numeric))
apply(data_out, 2,function(x) round(length(which(x==0))/length(x)*100,4))

normalise <-  function(x){
  if(length(unique(x))<=2){
    x <- x
  }else{
    x <- (x-min(x))/(max(x)-min(x))
  }
}

data_out <- as.data.frame(apply(data_out,2,function(x) normalise(x)))

names(data_out) <- paste("W",names(data_out),sep = "_")
Blanca$combodate <- as.Date(Blanca$combodate, "%Y-%m-%d")
Blanca$af_start <- as.Date(Blanca$af_start, "%Y-%m-%d")
Blanca$af_end <- as.Date(Blanca$af_end, "%Y-%m-%d")
follow_up <- Blanca$af_end-Blanca$af_start
event_time <- Blanca$combodate-Blanca$af_start
data_out$T.tilde <- ifelse(is.na(event_time),follow_up,ifelse(event_time>follow_up,follow_up,event_time))
data_out$Delta <- ifelse(is.na(event_time),0,ifelse(event_time>follow_up,0,1))
data_out$A <- ifelse(Blanca$anticoagulant=="NOAC",1,0)
data_out$T.tilde <- round(data_out$T.tilde/7)
data_out <- data_out[data_out$T.tilde<=200 & data_out$T.tilde>0,]



table(data_out$T.tilde)
hist(data_out$T.tilde)


table(data_out$Delta[data_out$T.tilde<100])/nrow(data_out)
table(data_out$A)/nrow(data_out)

idx <- sample(c(1:20000),1000)
data_out <- data_out[idx,]
data_out <- data_out[complete.cases(data_out),]


library(dplyr,abind)
library(tidyverse)
library(survival,simsurv)
library(survminer)
library(ggthemes,ggplot2)
library(survtmle2)
library(SuperLearner)
library(np)
library(extrafont)
library(MOSSATE)
detach('package:mgcv',unload=TRUE)
options(mc.cores = 8)
set.seed(25431, "L'Ecuyer-CMRG")
onestepfit <- MOSSATE::MOSS$new(data_out, dW = 1,verbose = TRUE, epsilon.step = 1e-10, max.iter = 1)

onestepfit$onestep_curve(g.SL.Lib = c("SL.gam"),
                         Delta.SL.Lib = c("SL.gam"),
                         ht.SL.Lib = c("SL.gam"),
                         env = parent.frame())


onestepdiff <- MOSSATE::MOSS_difference$new(data_out, verbose = TRUE, epsilon.step = 1e-2,
                                   max.iter = 50,
                                   ftimeMod = onestepfit$ftimeMod,
                                   ctimeMod = onestepfit$ctimeMod,
                                   trtMod = onestepfit$trtMod)
onestepdiff$initial_fit(env = parent.frame())
onestepdiff$onestep_diff_curve( )

onestepdiff$Psi.hat


onestepfit1 <- MOSS::MOSS$new(data_out, dW = 1,
                              verbose = TRUE, epsilon.step = 1e-2, max.iter = 50)
onestepfit1$onestep_curve(g.SL.Lib = c("SL.gam"),
                          Delta.SL.Lib = c("SL.gam"),
                          ht.SL.Lib = c("SL.gam"))

onestepfit0 <- MOSS::MOSS$new(data_out, dW = 0,
                              verbose = TRUE, epsilon.step = 1e-2, max.iter = 50)
onestepfit0$onestep_curve(g.SL.Lib = c("SL.gam","SL.glmnet"),
                          Delta.SL.Lib = c("SL.glmnet","SL.gam"),
                          ht.SL.Lib = c("SL.gam","SL.glmnet"))

fit_survtmle <- function(T.tilde, Delta, A, W_df) {
  t_0 <- max(T.tilde)
  fit <- survtmle(
    ftime = T.tilde,
    ftype = Delta,
    trt = A,
    adjustVars = W_df,
    SL.trt = c("SL.mean", "SL.glm", "SL.gam"),
    SL.ftime = c("SL.mean", "SL.glm", "SL.gam"),
    SL.ctime = c("SL.mean", "SL.glm", "SL.gam"),
    method = "hazard",
    returnIC = TRUE,
    verbose = FALSE
  )
  # extract cumulative incidence at each timepoint
  tpfit <- timepoints(fit, times = seq_len(t_0))
  len_groups <- as.numeric(
    unique(lapply(lapply(tpfit, FUN = `[[`, "est"), FUN = length))
  )
  names_groups <- unique(
    lapply(lapply(tpfit, FUN = `[[`, "est"), FUN = rownames)
  )[[1]]
  est_only <- t(matrix(unlist(
    lapply(tpfit, FUN = `[[`, "est")
  ), ncol = len_groups, byrow = TRUE))
  est_only <- as.data.frame(est_only)
  rownames(est_only) <- names_groups
  colnames(est_only) <- paste0("t", seq_len(ncol(est_only)))

  s_0 <- 1 - as.numeric(est_only[1, ])
  s_1 <- 1 - as.numeric(est_only[2, ])
  return(data.frame(time = 1:t_0, s_0 = s_0, s_1 = s_1))
}

SL <- fit_survtmle(data_out$T.tilde,data_out$Delta,data_out$A,data_out[,grep("W",names(data_out),value = T)])




do_once <- function(n_sim = 2e2) {
  simulated <- simulate_data(n_sim = n_sim)
  df <- simulated$dat
  true_surv <- simulated$true_surv1

  sl_lib_g <- c("SL.mean", "SL.glm", "SL.gam")
  sl_lib_censor <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
  sl_lib_failure <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
  range(df$T.tilde)
  df$T.tilde <- df$T.tilde + 1
  k_grid <- 1:max(df$T.tilde)

  message("KM")
  n_sample <- nrow(df)
  km_fit <- survfit(Surv(time = T.tilde, event = Delta) ~ A, data = df)
  surv1_km <- tail(km_fit$surv, km_fit$strata["A=1"])
  time1_km <- tail(km_fit$time, km_fit$strata["A=1"])
  surv0_km <- tail(km_fit$surv, km_fit$strata["A=0"])
  time0_km <- tail(km_fit$time, km_fit$strata["A=0"])
  library(zoo)
  impute_KM <- function(time, km) {
    surv1_km_final <- rep(NA, max(df$T.tilde))
    surv1_km_final[time] <- km
    surv1_km_final <- na.locf(surv1_km_final, na.rm = FALSE)
    surv1_km_final[is.na(surv1_km_final)] <- 1
    surv1_km_final <- c(1, surv1_km_final)
    surv1_km_final <- surv1_km_final[-length(surv1_km_final)]
    return(surv1_km_final)
  }
  surv1_km_final <- impute_KM(time = time1_km, km = surv1_km)
  surv0_km_final <- impute_KM(time = time0_km, km = surv0_km)
  km_fit_1 <- survival_curve$new(t = k_grid, survival = surv1_km_final)
  km_fit_0 <- survival_curve$new(t = k_grid, survival = surv0_km_final)

  message("SL")
  sl_fit <- initial_sl_fit(
    ftime = df$T.tilde,
    ftype = df$Delta,
    trt = df$A,
    adjustVars = data.frame(df[, c("W", "W1")]),
    t_0 = max(df$T.tilde),
    SL.trt = sl_lib_g,
    SL.ctime = sl_lib_censor,
    SL.ftime = sl_lib_failure
  )
  sl_fit$density_failure_1$hazard_to_survival()
  sl_fit$density_failure_0$hazard_to_survival()
  # WILSON hack no data is t_tilde = 2
  sl_fit$density_failure_1$t <- k_grid
  sl_fit$density_failure_0$t <- k_grid

  # ipcw
  message("ipcw + ee")
  ipcw_fit_1_all <- repeat_t_grid$new(
    method = ipcw,
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_1,
    density_censor = sl_fit$density_censor_1,
    g1W = sl_fit$g1W,
    A_intervene = 1
  )$fit(k_grid = k_grid)
  ipcw_fit_0_all <- repeat_t_grid$new(
    method = ipcw,
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_0,
    density_censor = sl_fit$density_censor_0,
    g1W = sl_fit$g1W,
    A_intervene = 0
  )$fit(k_grid = k_grid)
  ee_fit_1_all <- repeat_t_grid$new(
    method = ee,
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_1,
    density_censor = sl_fit$density_censor_1,
    g1W = sl_fit$g1W,
    A_intervene = 1
  )$fit(k_grid = k_grid)
  ee_fit_0_all <- repeat_t_grid$new(
    method = ee,
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_0,
    density_censor = sl_fit$density_censor_0,
    g1W = sl_fit$g1W,
    A_intervene = 0
  )$fit(k_grid = k_grid)
  ipcw_fit_1 <- survival_curve$new(t = k_grid, survival = ipcw_fit_1_all)
  ipcw_fit_0 <- survival_curve$new(t = k_grid, survival = ipcw_fit_0_all)
  ee_fit_1 <- survival_curve$new(t = k_grid, survival = ee_fit_1_all)
  ee_fit_0 <- survival_curve$new(t = k_grid, survival = ee_fit_0_all)

  sl_density_failure_1_marginal <- sl_fit$density_failure_1$clone(deep = TRUE)
  sl_density_failure_0_marginal <- sl_fit$density_failure_0$clone(deep = TRUE)
  sl_density_failure_1_marginal$survival <- matrix(colMeans(sl_density_failure_1_marginal$survival), nrow = 1)
  sl_density_failure_0_marginal$survival <- matrix(colMeans(sl_density_failure_0_marginal$survival), nrow = 1)
  message("moss")
  moss_fit <- MOSS$new(
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_1,
    density_censor = sl_fit$density_censor_1,
    g1W = sl_fit$g1W,
    A_intervene = 1,
    k_grid = k_grid
  )
  psi_moss_1 <- moss_fit$onestep_curve(
    epsilon = 1e-1 / n_sim,
    # epsilon = 1e-5,
    max_num_interation = 1e2,
    verbose = F
  )
  moss_fit <- MOSS$new(
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_0,
    density_censor = sl_fit$density_censor_0,
    g1W = sl_fit$g1W,
    A_intervene = 0,
    k_grid = k_grid
  )
  psi_moss_0 <- moss_fit$onestep_curve(
    epsilon = 1e-1 / n_sim,
    # epsilon = 1e-5,
    max_num_interation = 1e2,
    verbose = F
  )
  moss_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_1)
  moss_fit_0 <- survival_curve$new(t = k_grid, survival = psi_moss_0)

  # tmle
  message("tmle")
  tmle_fit <- tryCatch({
    tmle_fit <- fit_survtmle(
      T.tilde = df$T.tilde,
      Delta = df$Delta,
      A = df$A,
      W_df = data.frame(df[, c("W", "W1")])
    )
  },
  error = function(cond) {
    message("tmle error")
    NULL
  }
  )
  if (is.null(tmle_fit)) {
    tmle_fit_1 <- sl_density_failure_1_marginal$clone(deep = TRUE)
    tmle_fit_0 <- sl_density_failure_0_marginal$clone(deep = TRUE)
  } else {
    s_1 <- c(1, tmle_fit$s_1)
    s_1 <- s_1[-length(s_1)]
    s_0 <- c(1, tmle_fit$s_0)
    s_0 <- s_0[-length(s_0)]
    tmle_fit_1 <- survival_curve$new(t = k_grid, survival = s_1)
    tmle_fit_0 <- survival_curve$new(t = k_grid, survival = s_0)
  }
  message("moss with hazard submodel")
  moss_hazard_fit <- MOSS_hazard$new(
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_1,
    density_censor = sl_fit$density_censor_1,
    g1W = sl_fit$g1W,
    A_intervene = 1,
    k_grid = k_grid
  )
  psi_moss_hazard_1 <- moss_hazard_fit$iterate_onestep(epsilon = 1e-2, verbose = FALSE)
  moss_hazard_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_1)

  # evaluate against truth
  survival_truth_1 <- survival_curve$new(t = k_grid, survival = simulated$true_surv1(k_grid - 1))
  survival_truth_0 <- survival_curve$new(t = k_grid, survival = simulated$true_surv0(k_grid - 1))

  evaluate_moss <- evaluate_metric$new(
    survival = moss_fit_1, survival_truth = survival_truth_1
  )
  df_entropy_moss_1 <- evaluate_moss$evaluate_cross_entropy()
  df_entropy_moss_1$metric_name <- "cross_entropy"
  df_mse_moss_1 <- evaluate_moss$evaluate_mse()
  df_mse_moss_1$metric_name <- "mse"
  evaluate_moss_hazard <- evaluate_metric$new(
    survival = moss_hazard_fit_1, survival_truth = survival_truth_1
  )
  df_entropy_moss_hazard_1 <- evaluate_moss_hazard$evaluate_cross_entropy()
  df_entropy_moss_hazard_1$metric_name <- "cross_entropy"
  df_mse_moss_hazard_1 <- evaluate_moss_hazard$evaluate_mse()
  df_mse_moss_hazard_1$metric_name <- "mse"
  evaluate_sl <- evaluate_metric$new(
    survival = sl_density_failure_1_marginal, survival_truth = survival_truth_1
  )
  df_entropy_sl_1 <- evaluate_sl$evaluate_cross_entropy()
  df_entropy_sl_1$metric_name <- "cross_entropy"
  df_mse_sl_1 <- evaluate_sl$evaluate_mse()
  df_mse_sl_1$metric_name <- "mse"
  evaluate_ipcw <- evaluate_metric$new(
    survival = ipcw_fit_1, survival_truth = survival_truth_1
  )
  df_entropy_ipcw_1 <- evaluate_ipcw$evaluate_cross_entropy()
  df_entropy_ipcw_1$metric_name <- "cross_entropy"
  df_mse_ipcw_1 <- evaluate_ipcw$evaluate_mse()
  df_mse_ipcw_1$metric_name <- "mse"
  evaluate_ee <- evaluate_metric$new(
    survival = ee_fit_1, survival_truth = survival_truth_1
  )
  df_entropy_ee_1 <- evaluate_ee$evaluate_cross_entropy()
  df_entropy_ee_1$metric_name <- "cross_entropy"
  df_mse_ee_1 <- evaluate_ee$evaluate_mse()
  df_mse_ee_1$metric_name <- "mse"
  evaluate_tmle <- evaluate_metric$new(
    survival = tmle_fit_1, survival_truth = survival_truth_1
  )
  df_entropy_tmle_1 <- evaluate_tmle$evaluate_cross_entropy()
  df_entropy_tmle_1$metric_name <- "cross_entropy"
  df_mse_tmle_1 <- evaluate_tmle$evaluate_mse()
  df_mse_tmle_1$metric_name <- "mse"
  evaluate_km <- evaluate_metric$new(
    survival = km_fit_1, survival_truth = survival_truth_1
  )
  df_entropy_km_1 <- evaluate_km$evaluate_cross_entropy()
  df_entropy_km_1$metric_name <- "cross_entropy"
  df_mse_km_1 <- evaluate_km$evaluate_mse()
  df_mse_km_1$metric_name <- "mse"

  df_mse_moss_1$method <- "MOSS"
  df_mse_moss_hazard_1$method <- "MOSS_hazard"
  df_mse_sl_1$method <- "super learner"
  df_mse_ipcw_1$method <- "IPCW"
  df_mse_ee_1$method <- "EE"
  df_mse_tmle_1$method <- "TMLE"
  df_mse_km_1$method <- "KM"
  df_entropy_moss_1$method <- "MOSS"
  df_entropy_moss_hazard_1$method <- "MOSS_hazard"
  df_entropy_sl_1$method <- "super learner"
  df_entropy_ipcw_1$method <- "IPCW"
  df_entropy_ee_1$method <- "EE"
  df_entropy_tmle_1$method <- "TMLE"
  df_entropy_km_1$method <- "KM"
  df_plot <- plyr::rbind.fill(
    df_mse_moss_1,
    df_mse_moss_hazard_1,
    df_mse_sl_1,
    df_mse_ipcw_1,
    df_mse_ee_1,
    df_mse_tmle_1,
    df_mse_km_1,
    df_entropy_moss_1,
    df_entropy_moss_hazard_1,
    df_entropy_sl_1,
    df_entropy_ipcw_1,
    df_entropy_ee_1,
    df_entropy_tmle_1,
    df_entropy_km_1
  )
  return(df_plot)
}



#Plot KM curve
n.data <- nrow(data_out)
km.fit <- survfit(Surv(time = T.tilde, event = Delta) ~ A, data = data_out[data_out$T.tilde<=100,])
plot(km.fit, col=c('palegreen3','darkgreen'), lty = 2,main = paste('n=', n.data, '\n # of nuisance covariate = ?'),xlab = 'Time')


lines(onestepfit1$Psi.hat, type="l", cex=.2, col = 'red',lty=1) #This show the bias
lines(onestepfit0$Psi.hat, type="l", cex=.2, col = 'blue',lty=1) #This show the bias

lines(colMeans(onestepdiff$MOSS_A1$Qn.A1.t_initial), type="l", cex=.2, col = 'red',lty=4) #This show the bias
lines(colMeans(onestepdiff$MOSS_A1$Qn.A1.t_full), type="l", cex=.2, col = 'red',lty=1) #This show the bias
lines(colMeans(onestepdiff$MOSS_A0$Qn.A1.t_initial), type="l", cex=.2, col = 'blue',lty=4) #This show the bias
lines(colMeans(onestepdiff$MOSS_A0$Qn.A1.t_full), type="l", cex=.2, col = 'blue',lty=1) #This show the bias




get_result_2 <- function(scenario, scenario_tmle, scenario_EIC,scenario_var = data_out){
  TMLE_estimation <- plyr::ldply(scenario_tmle, function(x) x)
  EIC_estimation <- plyr::ldply(scenario_EIC, function(x) x)
  Estimation <- plyr::ldply(scenario, function(x) x)
  SD <- plyr::ldply(scenario, function(x) x) %>% sapply(sd) %>% as.vector()
  SD_t <- plyr::ldply(scenario_tmle, function(x) x) %>% sapply(sd) %>% as.vector()

  Estimation_average <-  plyr::ldply(scenario, function(x) colMeans(x))
  Estimation_average_t <- plyr::ldply(scenario_tmle, function(x) colMeans(x))

  return(list(
    TMLE_estimation= TMLE_estimation,
    EIC_estimation=EIC_estimation,
    Estimation=Estimation,
    scenario_var = scenario_var[,grep("W",names(scenario_var),value = T)],
    Summary=data.frame(time=1:100,SD=SD,SD_t=SD_t,
                       Mean = colMeans(Estimation_average) ,
                       TMLE = colMeans(Estimation_average_t))
  ))
}

scenario_p <- get_result_2(list(Estimation),list(TMLE_estimation),list(EIC),data_out)
raw_scenario <- data.frame(group = "Initial",scenario_p$Estimation)
tmle_scenario <- data.frame(group = "TMLE",scenario_p$TMLE_estimation)


plot(scenario_p$Summary$SD,type='l',lty=4)
lines(scenario_p$Summary$SD_t,type='l')

#---------------------Indentify HTE ---------------------------

IdentifyHTE <- function(scenario_p,raw_scenario,SL.library = c("SL.glm.interaction","SL.glmnet","SL.gam")){
  var <- scenario_p$scenario_var
  W_names <- grep('W',names(var),value = T)

  fitdat <- data.frame(Y=raw_scenario$X50,var)

  g.SL.1  <- SuperLearner(Y=fitdat$Y, X=fitdat[,W_names], SL.library=SL.library, family="gaussian")
  eval(parse(text=paste0("sample1 <- data.frame(",W_names[1],"=c(1:10000))")))
  eval(parse(text=paste0("sd1 <- data.frame(",W_names[1],"=c(1))")))
  nW <- length(W_names)
  for (i in 1:nW){
    pred.dat <- var
    pred.dat <- as.data.frame(t(rep(0,nW))) %>% dplyr::slice(rep(1:n(), each = 10000))
    names(pred.dat) <- W_names
    eval(parse(text=paste0("pred.dat$",W_names[i]," <- seq(-4,4,0.0008)[1:10000]")))
    g.pred <- predict(g.SL.1, pred.dat)
    eval(parse(text=paste0("sample1$",W_names[i]," <- g.pred$pred")))
    eval(parse(text=paste0("sd1$",W_names[i], " <- sample1$",W_names[i],"%>% IQR()")))
  }
  orderSD <- data.frame(SD=t(sd1),Covariate = names(sd1))
  orderSD <- orderSD[order(orderSD$SD),]

  plot(orderSD$SD,type="l",ylab = "SD")

  secondDerivative <- orderSD$SD
  for(i in 2:(nW-1)){
    secondDerivative[i] = 0
    secondDerivative[i] = orderSD$SD[i+1] + orderSD$SD[i-1] - 2 * orderSD$SD[i]
  }
  secondDerivative[1]  <- 0
  secondDerivative[nW] <- 0
  elbow.idx <- which.max(secondDerivative)
  elbow <- orderSD$SD[elbow.idx]
  orderSD$HTE <- ifelse(orderSD$SD>elbow,"Yes","No")
  HTE.data <-  gather(sample1)
  HTE.data$HTE <- orderSD[match(HTE.data$key,orderSD$Covariate),'HTE']

  return(list(HTE.data = HTE.data,
              orderSD = orderSD,
              elbow = elbow
  ))
}

Identified.HTE <- IdentifyHTE(scenario_p,raw_scenario)

Identified.HTE$orderSD

HTE.data <- Identified.HTE$HTE.data
HTE.data$key <- gsub("W","X",HTE.data$key)
orderSD <- Identified.HTE$orderSD
orderSD$Covariate <- gsub("W","X",orderSD$Covariate)

g = ggplot(data=HTE.data,aes(x=key,y=value)) + labs(x = "Covariate",y = "")+
  geom_boxplot(data=subset(HTE.data,HTE == 'No'),aes(fill = "#2F3131"),width=.5, alpha = 0.3,lwd=.3,outlier.colour = NA) +
  geom_boxplot(data=subset(HTE.data,HTE == 'Yes'),aes(fill = "#426E86"),width=.5, alpha = 0.7,lwd=.3,outlier.colour = NA) +
  scale_fill_manual("Contribute to HTE?",values= c( alpha(c("#2F3131","#426E86"), 0.3)),labels = c("No","Yes"))+ coord_flip()+
  scale_x_discrete(limits=orderSD$Covariate)+theme_hc() +
  geom_text(data=data.frame(), aes(x=orderSD$Covariate, y=0.05,label=round(orderSD$SD,3)), size=3)+
  theme(legend.spacing.y = unit(0, "mm"),
        legend.title=element_text(size=10),
        panel.border = element_rect(colour = "white", fill=NA),
        aspect.ratio = .9, axis.text = element_text(colour = 1, size = 10),
        plot.title = element_text(size = 9, family="Arial Bold"))




#Subgroup

subgroupCI <- function(scenario_p ,convery.idx,ii){
  sample.ii <- scenario_p$TMLE_estimation[convery.idx==ii,]
  psi.i <-  colMeans(sample.ii)
  sd_EIC <- apply(scenario_p$EIC_estimation[convery.idx==ii,], 2, sd)
  sd_EIC <- zoo::na.locf(sd_EIC, option = 'nocb')
  upper_CI <- psi.i + 1.95 * sd_EIC/sqrt(nrow(sample.ii))
  lower_CI <- psi.i - 1.95 * sd_EIC/sqrt(nrow(sample.ii))
  out <- data.frame(psi.i, sd_EIC, upper_CI, lower_CI)
  return(out)
}

compute_group <- function(scenario_p ,step,wname){
  controls <- scenario_p$scenario_var
  W_names <- grep('W',names(var),value = T)
  eval(parse(text=paste0("selection <- controls$",wname)))
  n1<-seq(0,max(selection)/2,max(selection)/2/step)
  convery.idx <- selection
  maEnbk1 <- list();
  for (ii in 1:step){
    I <- which(selection > n1[ii] & selection < n1[ii+1])
    maEnbk1[[ii]] <- I
    convery.idx[I] <- ii
  }
  table(convery.idx)
  convery.idx[which(convery.idx<1)] = round(convery.idx[which(convery.idx<1)])
  convery.idx <- ifelse(convery.idx==0,2,convery.idx)
  step <- max(convery.idx)
  return(list(step=step,
              convery.idx=convery.idx))
}

compute_effect <- function(scenario_p,step,wname,endtime){
  gte_group <- compute_group(scenario_p,step,wname)
  convery.idx <- gte_group$convery.idx
  step <- gte_group$step
  Psi_group <- list()
  for (ii in 1:step){
    sample.i <- subgroupCI(scenario_p,convery.idx,ii)
    sample.i$y <- sample.i$psi.i
    sample.i$x <- c(1:endtime)
    sample.i$group <- ii
    Psi_group[[ii]] <- sample.i
  }
  return(plyr::ldply(Psi_group, function(x) x))
}

kernal_process <- function(scenario,scenario_tmle,scenario_EIC,step,wname='W_mi_count',endtime=100){
  effect_by_sample <- list()
  scenario_p_v2 <- get_result_2(scenario,scenario_tmle,scenario_EIC,data_out)
  effect_by_sample[[i]] <- compute_effect(scenario_p = scenario_p_v2,step=step,wname,endtime)

  effect_by_sample_p <- plyr::ldply(effect_by_sample, function(x) x)
  effect_grouped_psi <- aggregate(effect_by_sample_p[,c("psi.i","sd_EIC","upper_CI","lower_CI")],
                                  by=list(time = effect_by_sample_p$x,group=effect_by_sample_p$group), mean)
  effect_grouped_sd <- aggregate(effect_by_sample_p[,c("psi.i","sd_EIC","upper_CI","lower_CI")],
                                 by=list(time = effect_by_sample_p$x,group=effect_by_sample_p$group), sd)
  effect_grouped_psi$sd_EIC <- effect_grouped_sd$sd_EIC
  effect_grouped_psi$upper_CI <- effect_grouped_psi$psi.i+1.96*effect_grouped_sd$sd_EIC/sqrt(50)
  effect_grouped_psi$lower_CI <- effect_grouped_psi$psi.i-1.96*effect_grouped_sd$sd_EIC/sqrt(50)
  return(effect_grouped_psi)
}

names(data_out)

kernal_p <- kernal_process(scenario = list(Estimation),
                             scenario_tmle = list(TMLE_estimation),
                             scenario_EIC = list(EIC),step = 10,
                             wname = 'W_age')

kernal_p <- data.frame(group = "Test",kernal_p)

kernal_plot <- rbind(kernal_p)
kernal_plot$group.1 <- kernal_plot$group.1/2



g=ggplot(kernal_plot) +
  geom_line(data=subset(kernal_plot,time == '12'), aes(x = group.1, y = psi.i,lty = group), size=.7)+
  geom_ribbon(data=subset(kernal_plot,time == '12'), aes(ymin=lower_CI,ymax=upper_CI,x=group.1,group=group),
              alpha=0.2)+
  scale_linetype_manual("",values=c("twodash", "solid","longdash","dotdash","dotted"))+
  theme_hc()+xlab("X\U2082") +ylab("CATE") +
  guides(lty = guide_legend(nrow = 2, byrow = TRUE))+
  theme(legend.spacing.y = unit(0, "mm"),
        panel.border = element_rect(colour = "white", fill=NA),
        aspect.ratio = .8, axis.text = element_text(colour = 1, size = 8),
        strip.background = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = NA,fill = "#FFFFFF"))

ggsave("CH3_kernal.png", plot = g, device = 'png', path = paste0(getwd(),'/GFX'),
       scale = 1, width = 150, height = 100, units = c("mm"), dpi = 250)




