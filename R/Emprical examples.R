data_out <- Elliott[,c('gender',"age","imd2015","paroxysmal","persistent","chronic",
                      "typical_flutter","atypical_flutter","af_hospital","mi_count","chf_count","pvd_count","cvd_count","dem_count",
                      "copd_count","rhe_count","giu_count","hp_count","ckd_count",
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
adjustVars <- names(data_out)
Elliott$combodate <- as.Date(Elliott$combodate, "%Y-%m-%d")
Elliott$af_start <- as.Date(Elliott$af_start, "%Y-%m-%d")
Elliott$af_end <- as.Date(Elliott$af_end, "%Y-%m-%d")
follow_up <- Elliott$af_end-Elliott$af_start
event_time <- Elliott$combodate-Elliott$af_start
data_out$T.tilde <- ifelse(is.na(event_time),follow_up,ifelse(event_time>follow_up,follow_up,event_time))
data_out$Delta <- ifelse(is.na(event_time),0,ifelse(event_time>follow_up,0,1))
data_out$A <- Elliott$anticoagulant
data_out$T.tilde <- round(data_out$T.tilde/180)
data_out <- data_out[data_out$T.tilde<=200 & data_out$T.tilde>0,]

table(data_out$T.tilde)
hist(data_out$T.tilde)


table(data_out$Delta[data_out$T.tilde<100])/nrow(data_out)
table(data_out$A)/nrow(data_out)

#data_out <- data_out[idx,]
data_out <- data_out[complete.cases(data_out),]
df <- data_out
sl_lib_g <- c("SL.mean", "SL.glm","SL.earth")
sl_lib_censor <- c("SL.mean", "SL.glm","SL.earth")
sl_lib_failure <- c("SL.mean", "SL.glm","SL.earth")
df$T.tilde <- df$T.tilde + 1
k_grid <- 1:max(df$T.tilde)
n_sim <- nrow(df)


library(MOSS)
message("SL")
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
ITE <- sl_fit$density_failure_1$survival-sl_fit$density_failure_0$survival


eic_fit_1 <- eic$new(
  A = df$A,
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  density_failure = sl_fit$density_failure_1,
  density_censor = sl_fit$density_censor_1,
  g1W = sl_fit$g1W,
  psi = colMeans(sl_fit$density_failure_1$survival),
  A_intervene = 1)$all_t(k_grid = k_grid)
mean_eic_1 <- colMeans(eic_fit_1)

eic_fit_0 <- eic$new(
  A = df$A,
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  density_failure = sl_fit$density_failure_0,
  density_censor = sl_fit$density_censor_0,
  g1W = sl_fit$g1W,
  psi = colMeans(sl_fit$density_failure_0$survival),
  A_intervene = 1)$all_t(k_grid = k_grid)
mean_eic_0 <- colMeans(eic_fit_0)



plot(sl_density_failure_0_marginal$survival %>% t(), lty = 1,type = 'l',col = 'blue')
lines(sl_density_failure_1_marginal$survival %>% t(), lty = 1,type = 'l',col = 'red')

compose_result <- function(scenario, scenario_tmle, scenario_EIC,scenario_var = data_out){
  if (exists('scenario')){
    Estimation <- plyr::ldply(scenario, function(x) x)
    SD <- plyr::ldply(scenario, function(x) x) %>% sapply(sd) %>% as.vector()
    Estimation_average <-  plyr::ldply(scenario, function(x) colMeans(x))
  }else{
    message('No scenario found')
  }
  if (exists('scenario_tmle')){
    TMLE_estimation <- plyr::ldply(scenario_tmle, function(x) x)
    SD_t <- plyr::ldply(scenario_tmle, function(x) x) %>% sapply(sd) %>% as.vector()
    Estimation_average_t <- plyr::ldply(scenario_tmle, function(x) colMeans(x))
  }else{
    message('No TMLE scenario found, set as initial fitting')
    scenario_tmle <- scenario
    TMLE_estimation <- plyr::ldply(scenario_tmle, function(x) x)
    SD_t <- plyr::ldply(scenario_tmle, function(x) x) %>% sapply(sd) %>% as.vector()
    Estimation_average_t <- plyr::ldply(scenario_tmle, function(x) colMeans(x))
  }
  if (exists('scenario_EIC')){
    EIC_estimation <- plyr::ldply(scenario_EIC, function(x) x)
  }else{
    message('No EIC found')
  }

  
  
  return(list(
    TMLE_estimation= TMLE_estimation,
    EIC_estimation=EIC_estimation,
    Estimation=Estimation,
    scenario_var = scenario_var[,grep("W",names(scenario_var),value = T)],
    Summary=data.frame(time=1:length(colMeans(Estimation_average)),SD=SD,SD_t=SD_t,
                       Mean = colMeans(Estimation_average) ,
                       TMLE = colMeans(Estimation_average_t))
  ))
}
Estimation <- ITE
TMLE_estimation <- ITE
EIC <- eic_fit_1-eic_fit_0

scenario_p <- compose_result(list(Estimation),list(TMLE_estimation),list(EIC),data_out)
raw_scenario <- data.frame(group = "Initial",scenario_p$Estimation)
tmle_scenario <- data.frame(group = "TMLE",scenario_p$TMLE_estimation)


plot(scenario_p$Summary$SD,type='l',lty=4)
lines(scenario_p$Summary$SD_t,type='l')




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

g = ggplot(data=violin.data.1,aes(x=group,y=value)) + labs(x = "",y = "Treatment Effect")+theme_hc()+
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

kernal_process <- function(scenario,scenario_tmle,scenario_EIC,step,wname='W_age',endtime=15){
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

names(data_out)

kernal_p <- kernal_process(scenario = list(Estimation),
                            scenario_tmle = list(TMLE_estimation),
                             scenario_EIC = list(EIC),step = 10,
                             wname = 'W_age')

kernal_p <- data.frame(group = "Test",kernal_p)

kernal_plot <- rbind(kernal_p)
kernal_plot$group.1 <- kernal_plot$group.1/10

plot(data.frame(X=data_out$W_age,Y=Estimation[,10]))
plot(data.frame(X=data_out$W_gender,Y=Estimation[,10]))
plot(data.frame(X=data_out$W_age,Y=scenario_p_v2$Estimation$`15`))


plot(data.frame(X=controls$W_age,Y=  Estimation[,14]))


g=ggplot(kernal_plot) +
  geom_line(data=subset(kernal_plot,time == '5'), aes(x = group.1, y = psi.i,lty = group), size=.7)+
  geom_ribbon(data=subset(kernal_plot,time == '5'), aes(ymin=lower_CI,ymax=upper_CI,x=group.1,group=group),
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




