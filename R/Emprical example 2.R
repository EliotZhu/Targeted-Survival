
data_out <- Blanca[,c('female',"age","imd2015","prev_vka","af_hospital","mi_count",
                      "chf_count","pvd_count","cvd_count","dem_count","copd_count","rhe_count","dtm_u_count","dtm_c_count",
                      "ckd_count","can_count","charlson_index_2y","vas_count","hyp_count",
                      "chadvasc_index_2y","bleed_count","alch_count","nsaid_asp_count","l_inr_count","hasbled_index_2y")]

Blanca <- Blanca_v1
adjustVars <- Blanca[,c('female',"age","imd2015","prev_vka","af_hospital","mi_count","paroxysmal","chronic",
                      "typical_flutter","atypical_flutter","esrf_count","mi_count","giu_count","liv_m_count",
                      "chf_count","pvd_count","cvd_count","dem_count","copd_count","rhe_count","dtm_u_count","dtm_c_count",
                      "ckd_count","can_count","charlson_index_2y","vas_count","hyp_count",
                      "stroke_count","Coronary.artery.operations_count","Antiarrythmics_count","Antidiabetics_count",
                      "Antihypertensives_count","Aspirin_count","Betablocker.Antiarrythmics_count","Betablockers_count",
                      "Clopidogrel_count","CP450_Inhibitors_count","Dipyridamole_count","Fibrinolytic.Drugs_count",
                      "NOAC_count","NSAID_count","Parenteral.Anticoagulants_count","Prasugrel_count","Proton.Pump.Inhibitors_count",
                      "Rifampicin_count","Selected.Anticonvulsants_count","Selective.Serotonin.Re.uptake.Inhibitors_count",
                      "Statins_count","Ticagrelor_count","Ticlopidine_count","VKA_count",
                      "chadvasc_index_2y","bleed_count","alch_count","nsaid_asp_count","l_inr_count","hasbled_index_2y")]

group_choice <-  Blanca[,c("age","prev_vka","esrf_count","Coronary.artery.operations_count","Aspirin_count",
                           "Prasugrel_count","chadvasc_index_2y")]

normalise <-  function(x){
  if(length(unique(x))<=2){
    x <- x
  }else if(length(which(x==0))/length(x)>=0.9){
    x <- ifelse(x==0,0,1)
  }else{
    x <- (x-min(x))/(max(x)-min(x))
  }
}


adjustVars[is.na(adjustVars)] <- 0
adjustVars <- as.data.frame(apply(adjustVars, 2, as.numeric))
apply(adjustVars, 2,function(x) round(length(which(x==0))/length(x)*100,4))
table(adjustVars$charlson_index_2y)

adjustVars <- as.data.frame(apply(adjustVars,2,function(x) normalise(x)))
names(adjustVars) <- paste("W",names(adjustVars),sep = "_")
adjustVars_name <- names(adjustVars)
data_out <- adjustVars


#Get time 
data_out$af_end <- as.Date(Blanca$af_end, "%Y-%m-%d")
data_out$af_start <- as.Date(Blanca$af_start, "%Y-%m-%d")
data_out$appendage_occlusion_date <- as.Date(Blanca$appendage_occlusion_date, "%Y-%m-%d")
data_out$cardioversion_date<- as.Date(Blanca$cardioversion_date, "%Y-%m-%d")
data_out$ablation_date <- as.Date(Blanca$ablation_date, "%Y-%m-%d")

data_out$appendage_occlusion_date[which(is.na(data_out$appendage_occlusion_date))] <- Blanca$af_end[which(is.na(data_out$appendage_occlusion_date))]
data_out$cardioversion_date[which(is.na(data_out$cardioversion_date))] <- Blanca$af_end[which(is.na(data_out$cardioversion_date))]
data_out$ablation_date[which(is.na(data_out$ablation_date))] <- Blanca$af_end[which(is.na(data_out$ablation_date))]

data_out$af_end[data_out$af_end>data_out$appendage_occlusion_date] <- Blanca$appendage_occlusion_date[data_out$af_end>data_out$appendage_occlusion_date] 
data_out$af_end[data_out$af_end>data_out$cardioversion_date] <- Blanca$cardioversion_date[data_out$af_end>data_out$cardioversion_date] 
data_out$af_end[data_out$af_end>data_out$ablation_date] <- Blanca$ablation_date[data_out$af_end>data_out$ablation_date] 

data_out$af_end <-  as.Date(Blanca$af_end, "%Y-%m-%d")
follow_up <- data_out$af_end-data_out$af_start


case_estimate <- function(df,adjustVars){
  require(dplyr)
  df$T.tilde <- round(df$T.tilde/180)
  df <- df[df$T.tilde<=200 & df$T.tilde>0,]
  
  n_sim <- nrow(df)
  k_grid <- 1:max(df$T.tilde)
  require(MOSS)
  message("SL")
  sl_lib_g <- c("SL.mean", "SL.glm","SL.earth")
  sl_lib_censor <- c("SL.mean", "SL.glm","SL.earth")
  sl_lib_failure <- c("SL.mean", "SL.glm","SL.earth")
  sl_fit <- initial_sl_fit(
    ftime = df$T.tilde,
    ftype = df$Delta,
    trt = df$A,
    adjustVars = data.frame(df[, names(adjustVars)]),
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
  eic_0 <- moss_fit$eic_out
  moss_fit0 <- moss_fit
  
  
  #plot(sl_density_failure_1_marginal$survival %>% t(), lty = 1,type = 'l',col = 'blue')
  #lines(sl_density_failure_0_marginal$survival %>% t(), lty = 1,type = 'l',col = 'red')
  #lines(moss_fit1$density_failure$survival %>% colMeans(), lty = 2,type = 'l',col = 'blue')
  #lines(moss_fit0$density_failure$survival %>% colMeans(), lty = 2,type = 'l',col = 'red')


  out <- list(moss_fit_1 = moss_fit1,
              moss_fit_0 = moss_fit0,
              sl_fit_1 = sl_fit$density_failure_1,
              sl_fit_0 = sl_fit$density_failure_0,
              sl_fit =sl_fit,
              df = df, 
              eic_1 = eic_1,
              eic_0 = eic_0)
  return(out)
} 
case_estimate_avg <- function(df,sl_fit){
  require(dplyr)
  
  message("moss diff")
  moss_hazard_ate_fit <- MOSS_hazard_ate$new(
    A = df$A,
    T_tilde = df$T.tilde,
    Delta = df$Delta,
    density_failure = sl_fit$density_failure_1,
    density_censor = sl_fit$density_censor_1,
    density_failure_0 = sl_fit$density_failure_0,
    density_censor_0 = sl_fit$density_censor_0,
    g1W = sl_fit$g1W,
    k_grid = 1:max(df$T.tilde)
  )
  psi_moss_hazard_ate_1 <- moss_hazard_ate_fit$iterate_onestep(epsilon = 1e-1 / n_sim, 
                                                               max_num_interation = 1e1, verbose = T)
  moss_hazard_ate_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_ate_1)
  
  out <- list(TMLE_diff = moss_hazard_ate_fit$density_failure$survival-moss_hazard_ate_fit$density_failure_0$survival,
              SL_diff = sl_fit$density_failure_1$survival-sl_fit$density_failure_0$survival,
              eic_diff = moss_hazard_ate_fit$eic_out)
  return(out)
}  







#Event combine
combodate <- as.Date(Blanca$combodate, "%Y-%m-%d")
event_time <- combodate-data_out$af_start
data_out$T.tilde <- ifelse(is.na(event_time),follow_up,ifelse(event_time>follow_up,follow_up,event_time))
data_out$Delta <- ifelse(is.na(event_time),0,ifelse(event_time>follow_up,0,1))
data_out$A <- ifelse(Blanca$anticoagulant=="NOAC",1,0)
df <- data_out[complete.cases(data_out),]

case_blanca_comb <- case_estimate(df)



#Event stroke
issedate <- as.Date(Blanca$issedate, "%Y-%m-%d")
event_time <- issedate-data_out$af_start
data_out$T.tilde <- ifelse(is.na(event_time),follow_up,ifelse(event_time>follow_up,follow_up,event_time))
data_out$Delta <- ifelse(is.na(event_time),0,ifelse(event_time>follow_up,0,1))
data_out$A <- ifelse(Blanca$anticoagulant=="NOAC",1,0)
df <- data_out[complete.cases(data_out),]

case_blanca_issedate <- case_estimate(df)


#Event stroke by age subgroups
'%notin%' <- function(x,y)!('%in%'(x,y))
issedate <- as.Date(Blanca$issedate, "%Y-%m-%d")
event_time <- issedate-data_out$af_start
data_out_sg <- data_out[,names(data_out)%notin%c('W_age')]
adjustVars2 <- adjustVars[,names(adjustVars)%notin%c('W_age')]
data_out_sg$T.tilde <- ifelse(is.na(event_time),follow_up,ifelse(event_time>follow_up,follow_up,event_time))
data_out_sg$Delta <- ifelse(is.na(event_time),0,ifelse(event_time>follow_up,0,1))
data_out_sg$A <- ifelse(Blanca$anticoagulant=="NOAC",1,0)
df_age1 <- data_out_sg[complete.cases(data_out_sg)&Blanca$age>=65,]
df_age2 <- data_out_sg[complete.cases(data_out_sg)&Blanca$age< 65,]

case_blanca_issedate_age1 <- case_estimate(df_age1,adjustVars)
case_blanca_issedate_age2 <- case_estimate(df_age2,adjustVars)



data_out_sg <- data_out[,names(data_out)%notin%c('W_age')]
adjustVars2 <- adjustVars[,names(adjustVars)%notin%c('W_age')]
data_out_sg$T.tilde <- ifelse(is.na(event_time),follow_up,ifelse(event_time>follow_up,follow_up,event_time))
data_out_sg$Delta <- ifelse(is.na(event_time),0,ifelse(event_time>follow_up,0,1))
data_out_sg$A <- ifelse(Blanca$anticoagulant=="NOAC",1,0)
df <- data_out_sg[complete.cases(data_out_sg),]
df_chadvasc1 <- df[Blanca$chadvasc_index_2y>=1 & Blanca$female ==0,]
df_chadvasc2 <- df[Blanca$chadvasc_index_2y>=2 & Blanca$female ==1,]


case_blanca_issedate_chadvasc1 <- case_estimate(df_chadvasc1,adjustVars)
case_blanca_issedate_chadvasc2 <- case_estimate(df_chadvasc2,adjustVars)



#Event stroke by prevka subgroups after
prevka_0 <- ITE[case_blanca_issedate$df$W_prev_vka==0,] %>% colMeans()
prevka_1 <- ITE[case_blanca_issedate$df$W_prev_vka==1,] %>% colMeans()
plot(prevka_1, lty = 1,type = 'l',col = 'red')
lines(prevka_0, lty = 1,type = 'l',col = 'blue')





#Event stroke by subgroups
'%notin%' <- function(x,y)!('%in%'(x,y))
issedate <- as.Date(Blanca$issedate, "%Y-%m-%d")
event_time <- issedate-data_out$af_start
data_out_sg <- data_out[,names(data_out)%notin%c('W_Coronary.artery.operations_count')]
adjustVars2 <- adjustVars[,names(adjustVars)%notin%c('W_Coronary.artery.operations_count')]
data_out_sg$T.tilde <- ifelse(is.na(event_time),follow_up,ifelse(event_time>follow_up,follow_up,event_time))
data_out_sg$Delta <- ifelse(is.na(event_time),0,ifelse(event_time>follow_up,0,1))
data_out_sg$A <- ifelse(Blanca$anticoagulant=="NOAC",1,0)
df_cao1 <- data_out_sg[complete.cases(data_out_sg)&Blanca$Coronary.artery.operations_count==1,]
df_cao0 <- data_out_sg[complete.cases(data_out_sg)&Blanca$Coronary.artery.operations_count==0,]

case_blanca_issedate_cao1 <- case_estimate(df_cao1,adjustVars2)
case_blanca_issedate_cao0 <- case_estimate(df_cao0,adjustVars2)



#CAO_1 <- ITE[case_blanca_issedate$df$W_Coronary.artery.operations_count==1,] %>% colMeans()
#CAO_0 <- ITE[case_blanca_issedate$df$W_Coronary.artery.operations_count==0,] %>% colMeans()
#plot(CAO_1, lty = 1,type = 'l',col = 'red')
#lines(CAO_0, lty = 1,type = 'l',col = 'blue')


#Event stroke by subgroups

issedate <- as.Date(Blanca$issedate, "%Y-%m-%d")
event_time <- issedate-data_out$af_start
data_out_sg <- data_out[,names(data_out)%notin%c('W_Prasugrel_count','W_Aspirin_count')]
adjustVars2 <- adjustVars[,names(adjustVars)%notin%c('W_Prasugrel_count','W_Aspirin_count')]
data_out_sg$T.tilde <- ifelse(is.na(event_time),follow_up,ifelse(event_time>follow_up,follow_up,event_time))
data_out_sg$Delta <- ifelse(is.na(event_time),0,ifelse(event_time>follow_up,0,1))
data_out_sg$A <- ifelse(Blanca$anticoagulant=="NOAC",1,0)
df <- data_out_sg[complete.cases(data_out_sg),]
df_ant1 <- df[Blanca$Prasugrel_count==1 | Blanca$Aspirin_count >0,]
df_ant0 <- df[Blanca$Prasugrel_count==0 & Blanca$Aspirin_count ==0,]

case_blanca_issedate_ant1 <- case_estimate(df_ant1,adjustVars2)
case_blanca_issedate_ant0 <- case_estimate(df_ant0,adjustVars2)




#ant_1 <- ITE[case_blanca_issedate$df$W_Prasugrel_count==1 | case_blanca_issedate$df$W_Aspirin_count>0,] %>% colMeans()
#ant_0 <- ITE[case_blanca_issedate$df$W_Prasugrel_count==0 | case_blanca_issedate$df$W_Aspirin_count==0,] %>% colMeans()
#plot(ant_1, lty = 1,type = 'l',col = 'red')
#lines(ant_0, lty = 1,type = 'l',col = 'blue')


rivaroxaban <- ifelse( is.na(Blanca$rivaroxaban),0,1)
dabigatran <- ifelse( is.na(Blanca$dabigatran),0,1)
apixaban <- ifelse( is.na(Blanca$apixaban),0,1)
edoxaban <- ifelse( is.na(Blanca$edoxaban),0,1)
rivaroxaban_ITE <- ITE[rivaroxaban[as.numeric(rownames(case_blanca_issedate$df))]==1 ,] %>% colMeans()
dabigatran_ITE <- ITE[dabigatran[as.numeric(rownames(case_blanca_issedate$df))]==1 ,] %>% colMeans()
apixaban_ITE <- ITE[apixaban[as.numeric(rownames(case_blanca_issedate$df))]==1 ,] %>% colMeans()
edoxaban_ITE <- ITE[edoxaban[as.numeric(rownames(case_blanca_issedate$df))]==1 ,] %>% colMeans()



#Event bleeding
mbdate <- as.Date(Blanca$mbdate, "%Y-%m-%d")
event_time <- mbdate-data_out$af_start
data_out$T.tilde <- ifelse(is.na(event_time),follow_up,ifelse(event_time>follow_up,follow_up,event_time))
data_out$Delta <- ifelse(is.na(event_time),0,ifelse(event_time>follow_up,0,1))
data_out$A <- ifelse(Blanca$anticoagulant=="NOAC",1,0)
df <- data_out[complete.cases(data_out) ,]

case_blanca_bleed <- case_estimate(df)




#Event death
dod <- as.Date(Blanca$dod, "%Y-%m-%d")
event_time <- dod-data_out$af_start
data_out$T.tilde <- ifelse(is.na(event_time),follow_up,ifelse(event_time>follow_up,follow_up,event_time))
data_out$Delta <- ifelse(is.na(event_time),0,ifelse(event_time>follow_up,0,1))
data_out$A <- ifelse(Blanca$anticoagulant=="NOAC",1,0)
df <- data_out[complete.cases(data_out),]

case_blanca_death <- case_estimate(df)




table(data_out$T.tilde)
hist(data_out$T.tilde)
table(data_out$Delta[data_out$T.tilde<100])/nrow(data_out)
table(data_out$A)/nrow(data_out)

#data_out <- data_out[idx,]








compose_result <- function(case = case_blanca_death){
  Estimation <- case$sl_fit_1$survival - case$sl_fit_0$survival
  SD <- Estimation %>% apply(2,sd) %>% as.vector()
  Estimation_average <-  colMeans(Estimation)
  
  TMLE_estimation <- case$moss_fit_1$density_failure$survival-case$moss_fit_0$density_failure$survival
  SD_t <- TMLE_estimation %>%apply(2,sd) %>% as.vector()
  Estimation_average_t <- colMeans(TMLE_estimation)
  
  EIC_estimation <- case$eic_1-case$eic_0

  
  return(list(
    Estimation=Estimation,
    TMLE_estimation= TMLE_estimation,
    EIC_estimation=EIC_estimation,
    Summary=data.frame(Time=1:ncol(Estimation), SD=SD,SD_t=SD_t,
                       Mean = Estimation_average ,
                       TMLE = Estimation_average_t
    )
  ))
}


issedate_chadvasc1 <- compose_result(case_blanca_issedate_chadvasc1)
issedate_chadvasc2 <- compose_result(case_blanca_issedate_chadvasc2)
issedate_age1 <- compose_result(case_blanca_issedate_age1)
issedate_age2 <- compose_result(case_blanca_issedate_age2)
issedate_cao1 <- compose_result(case_blanca_issedate_cao1)
issedate_cao0 <- compose_result(case_blanca_issedate_cao0)
issedate_ant1 <- compose_result(case_blanca_issedate_ant1)
issedate_ant0 <- compose_result(case_blanca_issedate_ant0)

issedate <- compose_result(case_blanca_issedate)
comb <- compose_result(case_blanca_comb)
bleed <- compose_result(case_blanca_bleed)
death <- compose_result(case_blanca_death)



diagnose.data <- rbind(issedate_chadvasc1$Summary,issedate_chadvasc2$Summary,
                       issedate_age1$Summary,issedate_age2$Summary,
                       issedate_cao1$Summary,issedate_cao0$Summary,
                       issedate_ant1$Summary,issedate_ant0$Summary,
                       issedate$Summary,comb$Summary,
                       bleed$Summary,death$Summary)



#---------------------Plot Large Violin---------------------
violin.data <-raw_scenario
names(violin.data) <- gsub("X","Time:",names(violin.data))
violin.data <- reshape2::melt(violin.data, id.vars=c("group"), value.name = "value")
violin.data$type <- "SL"

violin.tmle <- tmle_scenario
names(violin.tmle) <- gsub("X","Time:",names(violin.tmle))
violin.tmle <- reshape2::melt(violin.tmle, id.vars=c("group"), value.name = "value")
violin.tmle$type <- "MOSS"
names(violin.tmle) <- names(violin.data)

violin.data <- rbind(violin.tmle,violin.data)
violin.data.1 <- violin.data #%>% group_by(variable) %>% sample_n(1000)
violin.data.1 <- violin.data.1[violin.data.1$variable=="Time:12"|violin.data.1$variable=="Time:2"|violin.data.1$variable=="Time:6",]
table(violin.data.1$type)



vl.avg <- aggregate(violin.data$value, by=list(group = violin.data$group,variable=violin.data$variable,type = violin.data$type), mean)
vl.avg$value <- vl.avg$x
vl.avg <- vl.avg[vl.avg$variable=="Time:12"|vl.avg$variable=="Time:2"|vl.avg$variable=="Time:6",]

g = ggplot(data=violin.data.1,aes(x=group,y=value)) + labs(x = "",y = "Treatment Effect")+theme_hc()+ ylim(c(0,0.04))+
  geom_violin(data=subset(violin.data.1,type == 'MOSS'),aes(fill = "#F9BA32"),  alpha = 0.4,color=NA,lwd=.1,scale="width")+
  geom_violin(data=subset(violin.data.1,type == 'SL'),aes(fill = "#426E86"), alpha = 0.4, color=NA,lwd=.1,scale="width")+
  geom_point(data = subset(vl.avg,type == 'MOSS'),  alpha = .8, size = 2,color="#F9BA32",shape=3,stroke = 1)+
  geom_point(data = subset(vl.avg,type == 'SL'), alpha = .8, size = 2, color="#426E86",shape=3,stroke = 1)+
  scale_fill_manual("",values= c(alpha(c("#426E86","#F9BA32"),.4)),labels = c("ITE(SL)","ITE(TMLE)"))+
  facet_wrap(~variable,ncol=1,strip.position = "bottom")+coord_flip()+
  theme(legend.spacing.y = unit(0, "mm"),
        panel.border = element_rect(colour = "white", fill=NA),
        aspect.ratio = .3, axis.text = element_text(colour = 1, size = 10),
        strip.background = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = NA,fill = "#FFFFFF"))
ggsave("CH3_large_violin.png", plot = g, device = 'png', path = paste0(getwd(),'/GFX'),
       scale = 1, width = 130, height = 150, units = c("mm"), dpi = 300)




#---------------------Indentify HTE ---------------------------

IdentifyHTE <- function(scenario_p,raw_scenario,SL.library = c("SL.glm.interaction","SL.glmnet","SL.gam"),point){
  var <- scenario_p$scenario_var
  W_names <- grep('W',names(var),value = T)
  
  fitdat <- data.frame(Y=raw_scenario[,point],var)
  
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

Identified.HTE <- IdentifyHTE(scenario_p,raw_scenario,point=8)

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
  n1<-seq(min(selection),max(selection),length.out = step)
  convery.idx <- selection
  maEnbk1 <- list();
  for (ii in 1:step){
    I <- which(selection > n1[ii] & selection < n1[ii+1])
    maEnbk1[[ii]] <- I
    convery.idx[I] <- ii
  }
  table(convery.idx)
  step <- max(convery.idx)
  #plot(data.frame(X=controls$W_age,Y= convery.idx))
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

kernal_process <- function(scenario,scenario_tmle,scenario_EIC,step,wname='W_age',endtime=15,data_out){
  effect_by_sample <- list()
  scenario_p_v2 <- compose_result(scenario,scenario_tmle,scenario_EIC,data_out)
  effect_by_sample <- compute_effect(scenario_p = scenario_p_v2,step=step,wname,endtime)
  # effect_grouped_psi <- aggregate(effect_by_sample[,c("psi.i","sd_EIC","upper_CI","lower_CI")],
  #                                 by=list(time = effect_by_sample$x,group=effect_by_sample$group), mean)
  # 
  # effect_grouped_sd <- aggregate(effect_by_sample[,c("psi.i","sd_EIC","upper_CI","lower_CI")],
  #                                by=list(time = effect_by_sample$x,group=effect_by_sample$group), sd)
  effect_by_sample$time <- effect_by_sample$x
  
  effect_grouped_psi <- effect_by_sample
  effect_grouped_sd <- effect_by_sample
  
  apply(effect_by_sample,2,sd)
  effect_grouped_psi$sd_EIC <- effect_grouped_sd$sd_EIC
  effect_grouped_psi$upper_CI <- effect_grouped_psi$psi.i+1.96*effect_grouped_sd$sd_EIC/sqrt(50)
  effect_grouped_psi$lower_CI <- effect_grouped_psi$psi.i-1.96*effect_grouped_sd$sd_EIC/sqrt(50)
  return(effect_grouped_psi)
}


kernal_p <- kernal_process(scenario = list(Estimation),
                           scenario_tmle = list(TMLE_estimation),
                           scenario_EIC = list(EIC),step = 5,
                           wname = 'W_copd_count',
                           data_out=data_out)

kernal_p <- data.frame(group = "Test",kernal_p)

table(data_out$W_copd_count)

kernal_plot <- rbind(kernal_p)
kernal_plot$group.1 <- kernal_plot$group.1/10

plot(data.frame(X=data_out$W_age,Y=Estimation[,10]))
plot(data.frame(X=data_out$W_gender,Y=Estimation[,10]))
plot(data.frame(X=data_out$W_age,Y=scenario_p_v2$Estimation$`15`))


plot(data.frame(X=controls$W_age,Y=  Estimation[,14]))


g=ggplot(kernal_plot) +
  geom_line(data=subset(kernal_plot,time == '10'), aes(x = group.1, y = psi.i,lty = group), size=.7)+
  #geom_ribbon(data=subset(kernal_plot,time == '10'), aes(ymin=lower_CI,ymax=upper_CI,x=group.1,group=group),
  #            alpha=0.2)+
  scale_linetype_manual("",values=c("twodash", "solid","longdash","dotdash","dotted"))+
  theme_hc()+xlab("W_copd_count") +ylab("CATE") +
  guides(lty = guide_legend(nrow = 2, byrow = TRUE))+
  theme(legend.spacing.y = unit(0, "mm"),
        panel.border = element_rect(colour = "white", fill=NA),
        aspect.ratio = .8, axis.text = element_text(colour = 1, size = 8),
        strip.background = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = NA,fill = "#FFFFFF"))

ggsave("CH3_kernal.png", plot = g, device = 'png', path = paste0(getwd(),'/GFX'),
       scale = 1, width = 150, height = 100, units = c("mm"), dpi = 250)




