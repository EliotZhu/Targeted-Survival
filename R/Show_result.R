#---------------------Get data---------------------

raw_CH_3_By_Size_100 <- data.frame(group = "100",read_csv("result/raw_CH_3_By_Size_100.csv"))
raw_CH_3_By_Size_500 <- data.frame(group = "500",read_csv("result/raw_CH_3_By_Size_500.csv"))
raw_CH_3_By_Size_1000 <- data.frame(group = "1000",read_csv("result/raw_CH_3_By_Size_1000.csv"))
raw_CH_3_By_Size_2500 <- data.frame(group = "2500",read_csv("result/raw_CH_3_By_Size_2500.csv"))
raw_CH_3_By_Size_5000 <- data.frame(group = "5000",read_csv("result/raw_CH_3_By_Size_5000.csv"))
raw_CH_3_By_Size_10000 <- data.frame(group = "10000",read_csv("result/raw_CH_3_By_Size_10000.csv"))
raw_CH_3_By_Size_20000 <- data.frame(group = "20000",read_csv("result/raw_CH_3_By_Size_20000.csv"))


scenario_3 <- scenario_3_nohte <- scenario_3_r10 <- scenario_3_r01 <- scenario_3_r001 <- scenario_3_noc <-list()
for (i in 1:50){
  scenario_3 <- suppressMessages(append(scenario_3, list(read_csv(paste0("result/coverage/raw_CH_3_scenario_3_r30_",i,".csv"))) ))
  scenario_3_nohte <- suppressMessages(append(scenario_3_nohte, list(read_csv(paste0("result/coverage/raw_CH_3_scenario_3_nohte_",i,".csv"))) ))
  scenario_3_r10 <- suppressMessages(append(scenario_3_r10, list(read_csv(paste0("result/coverage/raw_CH_3_scenario_3_r10_",i,".csv"))) ))
  scenario_3_r01 <- suppressMessages(append(scenario_3_r01, list(read_csv(paste0("result/coverage/raw_CH_3_scenario_3_r01_",i,".csv"))) ))
  scenario_3_r001 <- suppressMessages(append(scenario_3_r001, list(read_csv(paste0("result/coverage/raw_CH_3_scenario_3_r001_",i,".csv"))) ))
  scenario_3_noc <- suppressMessages(append(scenario_3_noc, list(read_csv(paste0("result/coverage/raw_CH_3_scenario_3_noc_",i,".csv"))) ))
}

scenario_3_t <- scenario_3_nohte_t <- scenario_3_r10_t <- scenario_3_r01_t <-scenario_3_r001_t <-  scenario_3_noc_t <- list()
for (i in 1:50){
  scenario_3_t <- suppressMessages(append(scenario_3_t, list(read_csv(paste0("result/coverage/testdat_CH_3_scenario_3_r30_",i,".csv"))) ))
  scenario_3_r10_t <- suppressMessages(append(scenario_3_r10_t, list(read_csv(paste0("result/coverage/testdat_CH_3_scenario_3_r10_",i,".csv"))) ))
  scenario_3_r01_t <- suppressMessages(append(scenario_3_r01_t, list(read_csv(paste0("result/coverage/testdat_CH_3_scenario_3_r01_",i,".csv"))) ))
  scenario_3_r001_t <- suppressMessages(append(scenario_3_r001_t, list(read_csv(paste0("result/coverage/testdat_CH_3_scenario_3_r001_",i,".csv"))) ))
  scenario_3_nohte_t <- suppressMessages(append(scenario_3_nohte_t, list(read_csv(paste0("result/coverage/testdat_CH_3_scenario_3_nohte_",i,".csv"))) ))
  scenario_3_noc_t <- suppressMessages(append(scenario_3_noc_t, list(read_csv(paste0("result/coverage/testdat_CH_3_scenario_3_noc_",i,".csv"))) ))

}
scenario_3_tmle <- scenario_3_nohte_tmle <- scenario_3_r10_tmle <- scenario_3_r01_tmle <-scenario_3_r001_tmle <-  scenario_3_noc_tmle <- list()
for (i in 1:50){
  scenario_3_tmle <- suppressMessages(append(scenario_3_tmle, list(read_csv(paste0("result/coverage/tmle_CH_3_scenario_3_r30_",i,".csv"))) ))
  scenario_3_r10_tmle <- suppressMessages(append(scenario_3_r10_tmle, list(read_csv(paste0("result/coverage/tmle_CH_3_scenario_3_r10_",i,".csv"))) ))
  scenario_3_r01_tmle <- suppressMessages(append(scenario_3_r01_tmle, list(read_csv(paste0("result/coverage/tmle_CH_3_scenario_3_r01_",i,".csv"))) ))
  scenario_3_r001_tmle <- suppressMessages(append(scenario_3_r001_tmle, list(read_csv(paste0("result/coverage/tmle_CH_3_scenario_3_r001_",i,".csv"))) ))
  scenario_3_nohte_tmle <- suppressMessages(append(scenario_3_nohte_tmle, list(read_csv(paste0("result/coverage/tmle_CH_3_scenario_3_nohte_",i,".csv"))) ))
  scenario_3_noc_tmle <- suppressMessages(append(scenario_3_noc_tmle, list(read_csv(paste0("result/coverage/tmle_CH_3_scenario_3_noc_",i,".csv"))) ))
}
scenario_3_EIC <- scenario_3_nohte_EIC <- scenario_3_r10_EIC <- scenario_3_r01_EIC <- scenario_3_r001_EIC <- scenario_3_noc_EIC <- list()
for (i in 1:50){
  scenario_3_EIC <- suppressMessages(append(scenario_3_EIC, list(read_csv(paste0("result/coverage/sd_EIC_CH_3_scenario_3_r30_",i,".csv"))) ))
  scenario_3_r10_EIC <- suppressMessages(append(scenario_3_r10_EIC, list(read_csv(paste0("result/coverage/sd_EIC_CH_3_scenario_3_r10_",i,".csv"))) ))
  scenario_3_r01_EIC <- suppressMessages(append(scenario_3_r01_EIC, list(read_csv(paste0("result/coverage/sd_EIC_CH_3_scenario_3_r01_",i,".csv"))) ))
  scenario_3_r001_EIC <- suppressMessages(append(scenario_3_r001_EIC, list(read_csv(paste0("result/coverage/sd_EIC_CH_3_scenario_3_r001_",i,".csv"))) ))
  scenario_3_nohte_EIC <- suppressMessages(append(scenario_3_nohte_EIC, list(read_csv(paste0("result/coverage/sd_EIC_CH_3_scenario_3_nohte_",i,".csv"))) ))
  scenario_3_noc_EIC <- suppressMessages(append(scenario_3_noc_EIC, list(read_csv(paste0("result/coverage/sd_EIC_CH_3_scenario_3_noc_",i,".csv"))) ))
}


compute_true_effect <- function(x,ratDiv=70){
  rate <- as.data.frame(x[,5]+x[,6]+x[,7]+x[,8]+2*(x[,1]+x[,2]+x[,3]+x[,4]))
  s_diff_true_1 <-  apply(rate,1,function(x) exp(-x/ratDiv*seq(0,11,1))) %>% t()
  rate <- as.data.frame(x[,5]+x[,6]+x[,7]+x[,8]+x[,1]+x[,2]+x[,3]+x[,4])
  s_diff_true_0 <-  apply(rate,1,function(x) exp(-x/ratDiv*seq(0,11,1))) %>% t()
  return(s_diff_true_1-s_diff_true_0)
}



get_result <- function(scenario, scenario_t, scenario_tmle, scenario_EIC,size=50,ratDiv=50){
  TMLE_estimation <- plyr::ldply(scenario_tmle, function(x) x)
  EIC_estimation <- plyr::ldply(scenario_EIC, function(x) x)
  Estimation <- plyr::ldply(scenario, function(x) x)
  True.eff <- plyr::ldply(scenario_t, function(x) compute_true_effect(x[,grep("W",names(x),value = T)],ratDiv))

  SD <- plyr::ldply(scenario, function(x) x) %>% sapply(sd) %>% as.vector()

  Estimation_average <-  plyr::ldply(scenario, function(x) colMeans(x))
  Estimation_average_t <- plyr::ldply(scenario_tmle, function(x) colMeans(x))
  True_average <- plyr::ldply(scenario_t, function(x) colMeans(compute_true_effect(x[,grep("W",names(x),value = T)],ratDiv)))

  #Estimation_average <-  apply(Estimation,2,mean)
  #Estimation_average_t <-  apply(TMLE_estimation,2,mean)
  #True_average <- apply(True.eff,2,mean)

  #True.eff_average <-  Reduce("+", lapply(scenario_t, function(x) compute_true_effect(x[,grep("W",names(x),value = T)],ratDiv)))/size
  Bias <- Estimation-True.eff
  RMSE <- colMeans(Bias^2) %>% sqrt()
  RMSE <- RMSE/abs(colMeans(Estimation_average))

  Bias_t <- TMLE_estimation-True.eff
  RMSE_t <- colMeans(Bias_t^2) %>% sqrt()
  RMSE_t <- RMSE_t/abs(colMeans(Estimation_average_t))

  return(list(
    True.eff = True.eff,
    TMLE_estimation= TMLE_estimation,
    EIC_estimation=EIC_estimation,
    Estimation=Estimation,
    scenario_var=scenario_t,
    Summary=data.frame(time=1:12,SD=SD,RMSE=RMSE,RMSE_t=RMSE_t,
                       Mean = colMeans(Estimation_average) ,
                       TMLE = colMeans(Estimation_average_t) ,
                       pBias_t = paste0(round(colMeans((Estimation_average_t-True_average)/Estimation_average_t),4)*100,"%"),
                       pBias = paste0(round(colMeans((Estimation_average-True_average)/Estimation_average),4)*100,"%"),
                       true = colMeans(True_average)
                      )
  ))
}
scenario_3_p <- get_result(scenario_3,scenario_3_t, scenario_3_tmle, scenario_3_EIC,ratDiv=115)
scenario_3_r01_p <- get_result(scenario_3_r01,scenario_3_r01_t, scenario_3_r01_tmle, scenario_3_r01_EIC,ratDiv=4200)
scenario_3_r001_p <- get_result(scenario_3_r001,scenario_3_r001_t, scenario_3_r001_tmle, scenario_3_r001_EIC,ratDiv=25000)
scenario_3_r10_p <- get_result(scenario_3_r10,scenario_3_r10_t, scenario_3_r10_tmle, scenario_3_r10_EIC,ratDiv=430)
scenario_3_noc_p <- get_result(scenario_3_noc,scenario_3_noc_t, scenario_3_noc_tmle, scenario_3_noc_EIC,ratDiv=115)
scenario_3_nohte_p <- get_result(scenario_3_nohte,scenario_3_nohte_t, scenario_3_nohte_tmle, scenario_3_nohte_EIC,ratDiv=115)


raw_scenario_3 <- data.frame(group = "Event rate 30%",scenario_3_p$Estimation)
raw_scenario_3_r001 <- data.frame(group = "Event rate 0.1%",scenario_3_r001_p$Estimation)
raw_scenario_3_r01 <- data.frame(group = "Event rate 1%",scenario_3_r01_p$Estimation)
raw_scenario_3_r10 <- data.frame(group = "Event rate 10%",scenario_3_r10_p$Estimation)
raw_scenario_3_noc <- data.frame(group = "No confounding",scenario_3_noc_p$Estimation)
raw_scenario_3_nohte <- data.frame(group = "No HTE",scenario_3_nohte_p$Estimation)


true_scenario_3 <- data.frame(group = "Event rate 30%",scenario_3_p$True.eff)
true_scenario_3_r10 <- data.frame(group = "Event rate 10%",scenario_3_r10_p$True.eff)
true_scenario_3_r01 <- data.frame(group = "Event rate 1%",scenario_3_r01_p$True.eff)
true_scenario_3_r001 <- data.frame(group = "Event rate 0.1%",scenario_3_r001_p$True.eff)
true_scenario_3_noc <- data.frame(group = "No confounding",scenario_3_noc_p$True.eff)
true_scenario_3_nohte <- data.frame(group = "No HTE",scenario_3_nohte_p$True.eff)


tmle_scenario_3 <- data.frame(group = "Event rate 30%",scenario_3_p$TMLE_estimation)
tmle_scenario_3_r10 <- data.frame(group = "Event rate 10%",scenario_3_r10_p$TMLE_estimation)
tmle_scenario_3_r01 <- data.frame(group = "Event rate 1%",scenario_3_r01_p$TMLE_estimation)
tmle_scenario_3_r001 <- data.frame(group = "Event rate 0.1%",scenario_3_r001_p$TMLE_estimation)
tmle_scenario_3_noc <- data.frame(group = "No confounding",scenario_3_noc_p$TMLE_estimation)
tmle_scenario_3_nohte <- data.frame(group = "No HTE",scenario_3_nohte_p$TMLE_estimation)


