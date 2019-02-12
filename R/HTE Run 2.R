library(dplyr,abind)
library(tidyverse)
library(survival,simsurv)
library(survminer)
library(simcausal)

n_sim = 2000
simulated <- get.data(iti=1234,samplesize=n_sim, conmode ="scenario 3",ratDiv=100,confoundlevel = 1,confoundlevel_cen=1)
df <- simulated$dat
df <- df[complete.cases(df),]
adjustVars <- simulated$wnames

sl_lib_g <- c("SL.mean", "SL.glm", "SL.gam")
sl_lib_censor <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
sl_lib_failure <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")

df$T.tilde <- df$T.tilde + 1
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

#MOSS fitting
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
  epsilon = 10 / n_sim,
  #epsilon = 1e-2,
  max_num_interation = 5e1,
  verbose = T
)
moss_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_1)

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
  epsilon = 10 / n_sim,
  # epsilon = 1e-5,
  max_num_interation = 5e1,
  verbose = T
)
moss_fit_0 <- survival_curve$new(t = k_grid, survival = psi_moss_0)

controls <- simulated$dat2[,grep("W",names( simulated$dat2),value = T)]
true.1 <- t(apply(controls, 1, function(x) simulated$true_surv(x,(length(k_grid)-1)*10,A=1)))
true.0 <- t(apply(controls, 1, function(x) simulated$true_surv(x,(length(k_grid)-1)*10,A=0)))


#%Bias for treatment curve 
data.frame(TMLE= round((moss_fit_1$survival %>% t()-colMeans(true.1))/colMeans(true.1),4)*100,
           Sl = round((sl_density_failure_1_marginal$survival %>% t()-colMeans(true.1))/colMeans(true.1),4)*100) 
#%Bias for control curve 
data.frame(TMLE= round((moss_fit_0$survival %>% t()-colMeans(true.0))/colMeans(true.0),4)*100,
           Sl = round((sl_density_failure_0_marginal$survival %>% t()-colMeans(true.0))/colMeans(true.0),4)*100) 


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
psi_moss_hazard_ate_1 <- moss_hazard_ate_fit$iterate_onestep(
  epsilon = 1e-3, max_num_interation = 1e1, verbose = T)
moss_hazard_ate_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_ate_1)




TMLE_diff = psi_moss_hazard_ate_1
SL_diff = (sl_density_failure_1_marginal$survival-sl_density_failure_0_marginal$survival) %>% t()
true_diff = (colMeans(true.1)-colMeans(true.0))
#%Bias for difference curve
data.frame(TMLE= round(TMLE_diff-true_diff,4)/true_diff*100,Sl = round(SL_diff-true_diff,4)/true_diff*100)



