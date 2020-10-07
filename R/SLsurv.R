#Customised function

case_simu <- function(df, n_sim = 500, max_time = 30){
  require(dplyr)
  require(MOSS)
  require(survival)
  require(simcausal)
  require(abind)
  require(tidyverse)

  df <- df[df$T.tilde<= max_time & df$T.tilde>0,]
  df <- df[complete.cases(df),]


  adjustVars <- grep('w', colnames(df), value = T)
  sl_lib_g <- c( "SL.earth","SL.gam") #choose your own esemble algorithm here 
  sl_lib_censor <- c( "SL.earth","SL.gam")
  sl_lib_failure <- c( "SL.earth","SL.gam")

  #df$T.tilde <- df$T.tilde + 1
  k_grid <- 1:max(df$T.tilde)

  #SL
  sl_fit <- initial_sl_fit(
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    A = df$A,
    W = data.frame(df[, adjustVars]),
    #adjustVars = df[,c('W','W1')],
    t_max = max(df$T.tilde),
    sl_treatment = sl_lib_g,
    sl_censoring = sl_lib_censor,
    sl_failure = sl_lib_failure
  )


  sl_fit$density_failure_1$hazard_to_survival()
  sl_fit$density_failure_0$hazard_to_survival()
  sl_fit$density_failure_1$t <- k_grid
  sl_fit$density_failure_0$t <- k_grid

  sl_density_failure_1_marginal <- sl_fit$density_failure_1$clone(deep = TRUE)
  sl_density_failure_0_marginal <- sl_fit$density_failure_0$clone(deep = TRUE)
  sl_density_failure_1_marginal$survival <- matrix(colMeans(sl_density_failure_1_marginal$survival), nrow = 1)
  sl_density_failure_0_marginal$survival <- matrix(colMeans(sl_density_failure_0_marginal$survival), nrow = 1)


  out <- list(sl_fit_1 = sl_fit$density_failure_1$survival,
              sl_fit_0 = sl_fit$density_failure_0$survival,
              SL_diff = sl_fit$density_failure_1$survival-sl_fit$density_failure_0$survival)
  return(out)
}
