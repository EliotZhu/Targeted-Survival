library(dplyr,abind)
library(tidyverse)
library(survival,simsurv)
library(survminer)
library(simcausal)
library(here,usethis)

case_simu <- function(p,ratDiv){
  n_sim = 1000
  simulated <- get.data(iti=1234,samplesize=n_sim, conmode ="scenario 3",ratDiv=ratDiv,confoundlevel = 1,p=p)
  df <- simulated$dat
  df <- df[complete.cases(df),]
  table(df$Delta)/nrow(df)
  
  # hist(df$W1[df$A==1],col=rgb(0,0,1,0.5))
  # hist(df$W1[df$A==0],col=rgb(1,0,0,0.5),add = T)
  # hist(df$W2[df$A==1],col=rgb(0,0,1,0.5))
  # hist(df$W2[df$A==0],col=rgb(1,0,0,0.5),add = T)
  
  
  adjustVars <- simulated$wnames
  sl_lib_g <- c("SL.mean", "SL.glm", "SL.gam")
  sl_lib_censor <- c("SL.mean", "SL.glm", "SL.gam")
  sl_lib_failure <- c("SL.mean", "SL.glm", "SL.gam")
  
  
  #df$T.tilde <- df$T.tilde + 1
  k_grid <- 1:max(df$T.tilde)
  k_grid
  
  library(MOSS)
  #SL
  sl_fit <- initial_sl_fit(
    ftime = df$T.tilde,
    ftype = df$Delta,
    trt = df$A,
    adjustVars = data.frame(df[, adjustVars]),
    #adjustVars = df[,c('W','W1')],
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
  
  eic_1 <- moss_fit$eic_out
  
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
  eic_0 <- moss_fit$eic_out
  
  moss_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_1)
  moss_fit_0 <- survival_curve$new(t = k_grid, survival = psi_moss_0)
  
  
  controls <- simulated$dat2[,grep("W",names( simulated$dat2),value = T)]
  true.1 <- t(apply(controls, 1, function(x) simulated$true_surv(x,max(k_grid),A=1)))[,1:max(k_grid)]
  true.0 <- t(apply(controls, 1, function(x) simulated$true_surv(x,max(k_grid),A=0)))[,1:max(k_grid)]
  
  
  # #%Bias for treatment curve 
  # a <-  data.frame(MOSS= round((moss_fit_1$survival %>% t()-colMeans(true.1)),4),
  #                  Sl = round((sl_density_failure_1_marginal$survival %>% t()-colMeans(true.1)),4),
  #                  TR = colMeans(true.1)) 
  # #%Bias for control curve 
  # b <- data.frame(MOSS= round((moss_fit_0$survival %>% t()-colMeans(true.0)),4),
  #            Sl = round((sl_density_failure_0_marginal$survival %>% t()-colMeans(true.0)),4)) 
  # 
  # par(mfrow=c(3,1))
  # plot(abs(a$MOSS),type = 'l',lty=2,col='green')
  # lines(abs(a$Sl),type = 'l',lty=2)
  # 
  # plot(abs(b$MOSS),type = 'l',lty=2,col='green')
  # lines(abs(b$Sl),type = 'l',lty=2)
  
  
  
  moss_hazard_ate_fit <- MOSS_hazard_ate$new(
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_1,
    density_censor = sl_fit$density_censor_1,
    density_failure_0 = sl_fit$density_failure_0,
    density_censor_0 = sl_fit$density_censor_0,
    g1W = sl_fit$g1W,
    k_grid = k_grid
  )
  psi_moss_hazard_ate_1 <- moss_hazard_ate_fit$iterate_onestep(epsilon = 1e-1 / n_sim, 
                                                               max_num_interation = 2e1, verbose = T)
  moss_hazard_ate_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_ate_1)
  
  
  TMLE_diff = psi_moss_hazard_ate_1
  SL_diff = (sl_density_failure_1_marginal$survival-sl_density_failure_0_marginal$survival) %>% t()
  true_diff = (colMeans(true.1)-colMeans(true.0))
  
  # #%Bias for difference curve
  # c=data.frame(MOSS= round(TMLE_diff-true_diff,4),SL = round(SL_diff-true_diff,4))
  # plot(abs(c$MOSS),type = 'l',lty=2,col='green')
  # lines(abs(c$SL),type = 'l',lty=3)
  
  
  out <- list(moss_fit_1 = moss_fit_1,
              moss_fit_0 = moss_fit_0,
              sl_fit_1 = sl_fit$density_failure_1,
              sl_fit_0 = sl_fit$density_failure_0,
              true.1 = true.1,
              true.0 = true.0,
              eic_1 = eic_1,
              eic_0 = eic_0,
              TMLE_diff = TMLE_diff,
              SL_diff = sl_fit$density_failure_1$survival-sl_fit$density_failure_0$survival,
              true_diff = true.1-true.0,
              eic_diff = moss_hazard_ate_fit$eic_out
  )
  
  return(out)
}  


case1 <- case_simu(1,140)
case2 <- case_simu(5,300)
case3 <- case_simu(15,600)
case4 <- case_simu(30,800)