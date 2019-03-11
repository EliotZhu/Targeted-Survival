library(dplyr,abind)
library(tidyverse)
library(here,usethis)

case_simu <- function(n_sim,ratDiv,confoundlevel){
  require(MOSS)
  require(survival)
  require(simcausal)

  simulated <- get.data(iti=1234,samplesize=n_sim, conmode ="scenario 3",ratDiv=ratDiv,confoundlevel = 1)
  df <- simulated$dat
  df <- df[complete.cases(df),]
  cat(table(df$Delta)/nrow(df))
  
  # create function for fitting propensity score model

  
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
  
  
  controls <- simulated$dat2[,grep("W",names( simulated$dat2),value = T)]
  true.1 <- t(apply(controls, 1, function(x) simulated$true_surv(x,max(k_grid),A=1)))[,1:max(k_grid)]
  true.0 <- t(apply(controls, 1, function(x) simulated$true_surv(x,max(k_grid),A=0)))[,1:max(k_grid)]
  
  plot(true.1 %>% colMeans(),type = 'l')
  lines(  sl_density_failure_1_marginal$survival %>% t(),type = 'l',lty=2)
  plot(true.0 %>% colMeans(),type = 'l')
  lines(  sl_density_failure_0_marginal$survival %>% t(),type = 'l',lty=2)
  
  
  
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
  moss_fit1 <- moss_fit
  
  
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
  moss_fit0 <- moss_fit

  #moss_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_1)
  #moss_fit_0 <- survival_curve$new(t = k_grid, survival = psi_moss_0)
  
  

  
  
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
  psi_moss_hazard_ate_1 <- moss_hazard_ate_fit$iterate_onestep(epsilon = 1e-1 / n_sim, max_num_interation = 2e1, verbose = T)
  moss_hazard_ate_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_ate_1)
  
  
  TMLE_diff = psi_moss_hazard_ate_1
  SL_diff = (sl_density_failure_1_marginal$survival-sl_density_failure_0_marginal$survival) %>% t()
  true_diff = (colMeans(true.1)-colMeans(true.0))
  
  # #%Bias for difference curve
  # c=data.frame(MOSS= round(TMLE_diff-true_diff,4),SL = round(SL_diff-true_diff,4))
  # plot(abs(c$MOSS),type = 'l',lty=2,col='green')
  # lines(abs(c$SL),type = 'l',lty=3)
  
  out <- list(moss_fit_1 = moss_fit1,
              moss_fit_0 = moss_fit0,
              sl_fit_1 = sl_fit$density_failure_1,
              sl_fit_0 = sl_fit$density_failure_0,
              true.1 = true.1,
              true.0 = true.0,
              TMLE_diff = TMLE_diff,
              SL_diff = sl_fit$density_failure_1$survival-sl_fit$density_failure_0$survival,
              true_diff = true.1-true.0,
              eic_diff = moss_hazard_ate_fit$eic_fit,
              df = df,
              adjustVars = adjustVars
  )
  return(out)
}  
prop.func <- function(x, trt)
{
  # fit propensity score model
  propens.model <- SuperLearner(Y = trt, X = x, SL.library = c("SL.mean", "SL.glm", "SL.glmnet"),
                                family = "binomial")
  pi.x <- predict(propens.model, type = "response")
  pi.x$pred
}

case <- case1
check.overlap(x = data.frame(case$df[, case1$adjustVars]),
              trt = case$df$A,
              propensity.func = prop.func)

case1 <- case_simu(2000,220,0)
case2 <- case_simu(2000,220,0.5)
case3 <- case_simu(2000,220,1)

case1.1 <- case_simu(2000,500,0)
case2.1 <- case_simu(2000,500,0.5)
case3.1 <- case_simu(2000,500,1)
compose_result(case1.1)$Summary
compose_result(case2.1)$Summary
compose_result(case3.1)$Summary

case1.2 <- case_simu(2000,2000,0)
case2.2 <- case_simu(2000,2000,0.5)
case3.2 <- case_simu(2000,2000,1)

case1.3 <- case_simu(22000,0)
case2.3 <- case_simu(22000,1)
case3.3 <- case_simu(22000,2)
case4.3 <- case_simu(22000,3)



compose_result <- function(case){
    Estimation <- case$sl_fit_1$survival - case$sl_fit_0$survival
    SD <- Estimation %>% apply(2,sd) %>% as.vector()
    Estimation_average <-  colMeans(Estimation)

    TMLE_estimation <- case$moss_fit_1$density_failure$survival-case$moss_fit_0$density_failure$survival
    SD_t <- TMLE_estimation %>%apply(2,sd) %>% as.vector()
    Estimation_average_t <- colMeans(TMLE_estimation)

    EIC_estimation <- case$eic_diff
    True.eff <- case$true.1-case$true.0
    True_average <- colMeans(True.eff, na.rm = T)
    
    Bias <- Estimation-True.eff
    RMSE <- colMeans(Bias^2, na.rm = T) %>% sqrt()
    NRMSE <- RMSE/abs(Estimation_average)
    
    Bias_t <- TMLE_estimation-True.eff
    RMSE_t <- colMeans(Bias_t^2, na.rm = T) %>% sqrt()
    NRMSE_t <- RMSE_t/abs((Estimation_average_t))
    

    return(list(
      Estimation=Estimation,
      TMLE_estimation= TMLE_estimation,
      True.eff = True.eff,
      EIC_estimation=EIC_estimation,

      Summary=data.frame(Time=1:ncol(True.eff), SD=SD,SD_t=SD_t,
                         RMSE=NRMSE,RMSE_t=NRMSE_t,
                         Mean = Estimation_average ,
                         TMLE = Estimation_average_t ,
                         pBias = paste0(round(((Estimation_average-True_average)/Estimation_average),4)*100,"%"),
                         pBias_t = paste0(round(((Estimation_average_t-True_average)/Estimation_average_t),4)*100,"%"),
                         true = True_average
      )
    ))
}
scenario1 <- compose_result(case1) 
scenario2 <- compose_result(case2) 
scenario3 <- compose_result(case3) 

scenario1.1 <- compose_result(case1.1)
scenario2.1 <- compose_result(case2.1)
scenario3.1 <- compose_result(case3.1)

scenario1.2 <- compose_result(case1.2)
scenario2.2 <- compose_result(case2.2)
scenario3.2 <- compose_result(case3.2)

scenario1.3 <- compose_result(case1.3)
scenario2.3 <- compose_result(case2.3)
scenario3.3 <- compose_result(case3.3)
scenario4.3 <- compose_result(case4.3)

diagnose.data <- rbind(scenario1$Summary,scenario2$Summary,scenario3$Summary,scenario4$Summary,
                       scenario1.1$Summary,scenario2.1$Summary,scenario3.1$Summary,scenario4.1$Summary,
                       scenario1.2$Summary,scenario2.2$Summary,scenario3.2$Summary,scenario4.2$Summary,
                       scenario1.3$Summary,scenario2.3$Summary,scenario3.3$Summary,scenario4.3$Summary)
diagnose.data.out.2 <-diagnose.data[diagnose.data$Time=="2",]
diagnose.data.out.9 <-diagnose.data[diagnose.data$Time=="9",]
diagnose.data.out.18 <-diagnose.data[diagnose.data$Time=="18",]








write.csv(diagnose.data, file = paste0(getwd(),'/Result/diagnose.data.csv') , row.names = FALSE)




