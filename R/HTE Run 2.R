
simulated <- get.data(iti=51234,samplesize=1000, conmode ="scenario 3", endtime=50,ratDiv=400)
df <- simulated$dat
true_surv <- simulated$true_surv1
adjustVars <- simulated$wnames
nsim = 1000

sl_lib_g <- c("SL.mean", "SL.glm", "SL.gam")
sl_lib_censor <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
sl_lib_failure <- c("SL.mean", "SL.glm", "SL.gam", "SL.earth")
eventrate <- table(data_out$Delta[data_out$T.tilde<12])/nrow(data_out)
df$T.tilde <- df$T.tilde + 1
k_grid <- 1:max(df$T.tilde)

library(MOSS)

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
  adjustVars = data.frame(df[, adjustVars]),
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




