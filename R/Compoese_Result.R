
#---------------------Large Violin 1---------------------
violin.true <- rbind(true_scenario_3,
                     true_scenario_3_noc,
                     true_scenario_3_r10,
                     true_scenario_3_r01,
                     true_scenario_3_r001)
names(violin.true) <- gsub("X","Time:",names(violin.true))
violin.true <- reshape2::melt(violin.true, id.vars=c("group"), value.name = "value")
violin.true$type <- "True"

violin.data <- rbind(raw_scenario_3,
                     raw_scenario_3_noc,
                     raw_scenario_3_r10,
                     raw_scenario_3_r01,
                     raw_scenario_3_r001)
names(violin.data) <- gsub("V","Time:",names(violin.data))
violin.data <- reshape2::melt(violin.data, id.vars=c("group"), value.name = "value")
violin.data$type <- "Raw"

violin.tmle <- rbind(tmle_scenario_3,
                     tmle_scenario_3_noc,
                     tmle_scenario_3_r10,
                     tmle_scenario_3_r01,
                     tmle_scenario_3_r001)
names(violin.tmle) <- gsub("V","Time:",names(violin.tmle))
violin.tmle <- reshape2::melt(violin.tmle, id.vars=c("group"), value.name = "value")
violin.tmle$type <- "tmle"


#---------------------Plot Large Violin---------------------
violin.data <- rbind(violin.data,violin.tmle,violin.true)
violin.data.1 <- violin.data %>% group_by(variable) #%>% sample_n(1000)
violin.data.1 <- violin.data.1[violin.data.1$variable=="Time:12"|violin.data.1$variable=="Time:1"|violin.data.1$variable=="Time:6",]

table(violin.data.1$group)


vl.avg <- aggregate(violin.data$value, by=list(group = violin.data$group,variable=violin.data$variable,type = violin.data$type), mean)
vl.avg$value <- vl.avg$x
vl.avg <- vl.avg[vl.avg$variable=="Time:12"|vl.avg$variable=="Time:1"|vl.avg$variable=="Time:6",]


g = ggplot(data=violin.data.1,aes(x=group,y=value)) + labs(x = "",y = "Treatment Effect")+theme_hc()+ylim(-.3,0.05)+
  geom_violin(data=subset(violin.data.1,type == 'True'),aes(fill = "#2F3131"),  alpha = 0.4, color=NA,lwd=.1,scale="width")+
  geom_violin(data=subset(violin.data.1,type == 'tmle'),aes(fill = "#F9BA32"),  alpha = 0.4,color=NA,lwd=.1,scale="width")+
  geom_violin(data=subset(violin.data.1,type == 'Raw'),aes(fill = "#426E86"), alpha = 0.4, color=NA,lwd=.1,scale="width")+
  geom_point(data = subset(vl.avg,type == 'True'),size = 2, color="#2F3131",shape=3,stroke = 1)+
  geom_point(data = subset(vl.avg,type == 'tmle'),  alpha = .8, size = 2,color="#F9BA32",shape=3,stroke = 1)+
  geom_point(data = subset(vl.avg,type == 'Raw'), alpha = .8, size = 2, color="#426E86",shape=3,stroke = 1)+
  scale_fill_manual("",values= c(alpha(c("#2F3131","#426E86","#F9BA32"),.4)),labels = c("ITE(True)","ITE(SL)","ITE(TMLE)"))+
  scale_x_discrete(limits = c("No confounding","Event rate 30%","Event rate 10%","Event rate 1%","Event rate 0.1%"))+
  facet_wrap(~variable,ncol=1,strip.position = "bottom")+coord_flip()+
  theme(legend.spacing.y = unit(0, "mm"),
        panel.border = element_rect(colour = "white", fill=NA),
        aspect.ratio = .3, axis.text = element_text(colour = 1, size = 10),
        strip.background = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = NA,fill = "#FFFFFF"))
ggsave("CH3_large_violin.png", plot = g, device = 'png', path = paste0(getwd(),'/GFX'),
       scale = 1, width = 130, height = 150, units = c("mm"), dpi = 300)



# g = ggplot(data=violin.data.1,aes(x=group,y=value)) + labs(x = "",y = "Treatment Effect")+theme_hc()+ylim(-.4,0.1)+
#   geom_violin(data=subset(violin.data.1,type == 'tmle'),aes(fill = "#FFB11B"),  alpha = 0.3,color="#FFB11B",lwd=.1,scale="width")+
#   geom_boxplot(data=subset(violin.data.1,type == 'tmle'),aes(fill = "#FFB11B"),width=.25, alpha = 0.1,color="#FFB11B",lwd=.1,outlier.colour = NA) +
#   geom_violin(data=subset(violin.data.1,type == 'True'),aes(fill = "#2B5F75"),  alpha = 0.1, color="#2B5F75",lwd=.1,scale="width")+
#   geom_boxplot(data=subset(violin.data.1,type == 'True'),aes(fill = "#2B5F75"),width=.25, alpha = 0.1,color="#2B5F75",lwd=.1,outlier.colour = NA) +
#   #scale_x_discrete(limits=c("Scenario 1","Scenario 2","Scenario 3","Scenario 4","Scenario 5"))+
#   scale_fill_manual("",values= c( alpha(c("#2B5F75","#FFB11B"), 0.1)),labels = c("ITE(True)","ITE(TMLE)"))+ coord_flip()+
#   facet_wrap(~variable,ncol=3)+coord_flip()+
#   theme(legend.spacing.y = unit(0, "mm"),
#         panel.border = element_rect(colour = "white", fill=NA),
#         aspect.ratio = 1, axis.text = element_text(colour = 1, size = 10),
#         legend.background = element_blank(),
#         legend.box.background = element_rect(colour = NA,fill = "#FFFFFF"))
# ggsave("CH3_violin_tmle.png", plot = g, device = 'png', path = paste0(getwd(),'/GFX'),
#        scale = 1, width = 300, height = 120, units = c("mm"), dpi = 300)
#

#---------------------RMSE by SIZE---------------------
summary_3_r10 <- data.frame(group = 'Event rate 10%',scenario_3_r10_p$Summary)
summary_3_r01 <- data.frame(group = 'Event rate 1%',scenario_3_r01_p$Summary)
summary_3_r001 <- data.frame(group = 'Event rate 0.1%',scenario_3_r001_p$Summary)
summary_3_r30 <- data.frame(group = 'Event rate 30%',scenario_3_p$Summary)
summary_3_noc <- data.frame(group = 'Event rate 30%, no confounding',scenario_3_noc_p$Summary)

diagnose.data <- rbind(summary_3_r10,summary_3_r30,summary_3_r01,summary_3_r001,summary_3_noc)
diagnose.data.out <-diagnose.data[diagnose.data$time=="12",]
write.csv(diagnose.data, file = paste0(getwd(),'/result/diagnose.data.csv') , row.names = FALSE)

#---------------------RMSE---------------------
g=ggplot(diagnose.data, aes(x = time, y = RMSE, lty = group)) +
  geom_line(aes(colour = group),size = 1)+
  theme_hc()+
  xlab("Time point") +
  ylab("RMSE") +
  scale_x_discrete(limits=c(1:12))+
  scale_colour_economist()+
  theme_hc()+
  theme(legend.spacing.y = unit(0, "mm"),
        panel.border = element_rect(colour = "white", fill=NA),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 10),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = NA,fill = "#FFFFFF"))+
  guides(col = guide_legend(nrow = 3, byrow = TRUE))
ggsave("CH3_rmse_time.png", plot = g, device = 'png', path = paste0(getwd(),'/GFX'),
       scale = 1, width = 100, height = 100, units = c("mm"), dpi = 300)


#---------------------SD---------------------
g=ggplot(diagnose.data, aes(x = time, y = SD, lty = group)) +
  geom_line(aes(colour = group))+
  theme_hc()+
  xlab("Time point") +
  ylab("Standard Deviation") +
  scale_x_discrete(limits=c(1:12))+
  scale_colour_economist()+
  theme_hc()+
  theme(legend.spacing.y = unit(0, "mm"),
        panel.border = element_rect(colour = "white", fill=NA),
        aspect.ratio = 1, axis.text = element_text(colour = 1, size = 10),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = NA,fill = "#FFFFFF"))+
  guides(col = guide_legend(nrow = 3, byrow = TRUE))
ggsave("CH3_SE_time.png", plot = g, device = 'png', path = paste0(getwd(),'/GFX'),
       scale = 1, width = 100, height = 100, units = c("mm"), dpi = 300)


#---------------------Indentify HTE ---------------------------
scenario_3_HTE1 <- scenario_3_HTE2 <- scenario_3_HTE3 <-list()
for (i in 1:15){
  scenario_3_HTE1 <- suppressMessages(append(scenario_3_HTE1, list(read_csv(paste0("result/coverage/raw_CH_3_scenario_3_HTE1_",i,".csv"))) ))
  #scenario_3_HTE2 <- suppressMessages(append(scenario_3_HTE2, list(read_csv(paste0("result/coverage/raw_CH_3_scenario_3_HTE2_",i,".csv"))) ))
  #scenario_3_HTE3 <- suppressMessages(append(scenario_3_HTE3, list(read_csv(paste0("result/coverage/raw_CH_3_scenario_3_HTE3_",i,".csv"))) ))
}
scenario_3_HTE1_t <- scenario_3_HTE2_t <- scenario_3_HTE3_t <- list()
for (i in 1:15){
  scenario_3_HTE1_t <- suppressMessages(append(scenario_3_HTE1_t, list(read_csv(paste0("result/coverage/testdat_CH_3_scenario_3_HTE1_",i,".csv"))) ))
  #scenario_3_HTE2_t <- suppressMessages(append(scenario_3_HTE2_t, list(read_csv(paste0("result/coverage/testdat_CH_3_scenario_3_HTE2_",i,".csv"))) ))
  #scenario_3_HTE3_t <- suppressMessages(append(scenario_3_HTE3_t, list(read_csv(paste0("result/coverage/testdat_CH_3_scenario_3_HTE3_",i,".csv"))) ))
}
scenario_3_HTE1_tmle <- scenario_3_HTE2_tmle <- scenario_3_HTE3_tmle <- list()
for (i in 1:15){
  scenario_3_HTE1_tmle <- suppressMessages(append(scenario_3_HTE1_tmle, list(read_csv(paste0("result/coverage/tmle_CH_3_scenario_3_HTE1_",i,".csv"))) ))
  #scenario_3_HTE2_tmle <- suppressMessages(append(scenario_3_HTE2_tmle, list(read_csv(paste0("result/coverage/tmle_CH_3_scenario_3_HTE2_",i,".csv"))) ))
  #scenario_3_HTE3_tmle <- suppressMessages(append(scenario_3_HTE3_tmle, list(read_csv(paste0("result/coverage/tmle_CH_3_scenario_3_HTE3_",i,".csv"))) ))
}
scenario_3_HTE1_EIC <- scenario_3_HTE2_EIC <- scenario_3_HTE3_EIC <- list()
for (i in 1:15){
  scenario_3_HTE1_EIC <- suppressMessages(append(scenario_3_HTE1_EIC, list(read_csv(paste0("result/coverage/sd_EIC_CH_3_scenario_3_HTE1_",i,".csv"))) ))
  #scenario_3_HTE2_EIC <- suppressMessages(append(scenario_3_HTE2_EIC, list(read_csv(paste0("result/coverage/sd_EIC_CH_3_scenario_3_HTE2_",i,".csv"))) ))
  #scenario_3_HTE3_EIC <- suppressMessages(append(scenario_3_HTE3_EIC, list(read_csv(paste0("result/coverage/sd_EIC_CH_3_scenario_3_HTE3_",i,".csv"))) ))
}
scenario_3_HTE1 <- get_result(scenario_3_HTE1,scenario_3_HTE1_t, scenario_3_HTE1_tmle, scenario_3_HTE1_EIC,ratDiv=250)
scenario_3_HTE2 <- get_result(scenario_3_HTE2,scenario_3_HTE2_t, scenario_3_HTE2_tmle, scenario_3_HTE2_EIC,ratDiv=250)
scenario_3_HTE3 <- get_result(scenario_3_HTE3,scenario_3_HTE3_t, scenario_3_HTE3_tmle, scenario_3_HTE3_EIC,ratDiv=250)

raw_scenario_3_HTE3 <- data.frame(group = "Event rate 30%",scenario_3_HTE3$Estimation)
raw_scenario_3_HTE2 <- data.frame(group = "Event rate 30%",scenario_3_HTE2$Estimation)
raw_scenario_3_HTE1 <- data.frame(group = "Event rate 30%",scenario_3_HTE1$Estimation)


scenario_p <- scenario_3_HTE1
raw_scenario <- raw_scenario_3_HTE1

IdentifyHTE <- function(scenario_p,raw_scenario,SL.library = c("SL.glm.interaction","SL.glmnet","SL.gam")){
  var <-plyr::ldply(scenario_p$scenario_var, function(x) x)
  W_names <- grep('W',names(var),value = T)
  
  index <- sample(1:25000,10000)
  fitdat <- data.frame(Y=raw_scenario$V12,var)[index,]
  
  g.SL.1  <- SuperLearner(Y=fitdat$Y, X=fitdat[,W_names], SL.library=SL.library, family="gaussian")
  sample1 <- data.frame(W1=c(1:10000))
  sd1 <- data.frame(W1=c(1))
  nW <- length(W_names)
  for (i in 1:nW){
    pred.dat <- var
    pred.dat <- as.data.frame(t(rep(0,nW))) %>% dplyr::slice(rep(1:n(), each = 10000))
    names(pred.dat) <- W_names
    eval(parse(text=paste0("pred.dat$W",i," <- seq(-4,4,0.0008)[1:10000]")))
    g.pred <- predict(g.SL.1, pred.dat)
    eval(parse(text=paste0("sample1$W",i," <- g.pred$pred")))
    eval(parse(text=paste0("sd1$W",i, " <- sample1$W",i,"%>% IQR()")))
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


#plot(as.numeric(orderSD$SD),type="l")
#plot(as.numeric(sd1^2-sd0^2),type="l")

IdentifiedHTE.glmnet <- IdentifyHTE(scenario_3_p,raw_scenario_3,SL.library = c("SL.glmnet"))
IdentifiedHTE.glm <- IdentifyHTE(scenario_3_p,raw_scenario_3,SL.library = c("SL.glm"))
IdentifiedHTE.gam <- IdentifyHTE(scenario_3_p,raw_scenario_3,SL.library = c("SL.gam"))
IdentifiedHTE.earth <- IdentifyHTE(scenario_3_p,raw_scenario_3,SL.library = c("SL.earth"))
write.csv(diagnose.data, file = paste0(getwd(),'/result/diagnose.data.csv') , row.names = FALSE)



IdentifiedHTE.large <- IdentifyHTE(scenario_3_HTE3,raw_scenario_3_HTE3)

IdentifiedHTE.large$orderSD

IdentifiedHTE <- IdentifyHTE(scenario_3_p,raw_scenario_3)
HTE.data.3 <- IdentifiedHTE$HTE.data
HTE.data.3$key <- gsub("W","X",HTE.data.3$key)
orderSD.3 <- IdentifiedHTE$orderSD
orderSD.3$Covariate <- gsub("W","X",orderSD.3$Covariate)
HTE.data <- HTE.data.3
orderSD <- orderSD.3



IdentifiedHTE <- IdentifyHTE(scenario_3_nohte_p,raw_scenario_3_nohte)
HTE.data.nohte <- IdentifiedHTE$HTE.data
HTE.data.nohte$key <- gsub("W","X",HTE.data.nohte$key)
orderSD.nohte <- IdentifiedHTE$orderSD
orderSD.nohte$Covariate <- gsub("W","X",orderSD.nohte$Covariate)
HTE.data <- HTE.data.nohte
orderSD <- orderSD.nohte


g = ggplot(data=HTE.data,aes(x=key,y=value)) + labs(x = "Covariate",y = "")+ ylim(-.06,0.06)+
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

ggsave("CH3_HTE_box.png", plot = g, device = 'png', path = paste0(getwd(),'/GFX'),
       scale = 1, width = 110, height = 100, units = c("mm"), dpi = 300)

ggsave("CH3_HTE_nohte.png", plot = g, device = 'png', path = paste0(getwd(),'/GFX'),
       scale = 1, width = 110, height = 100, units = c("mm"), dpi = 300)


#purmuted KS
ComputeKS <- function(scenario_p){
  var <-plyr::ldply(scenario_p$scenario_var, function(x) x)
  W_names <- grep('W',names(var),value = T)
  
  index <- sample(seq(1:10000),5000)
  sample1 <- data.frame(scenario_p$TMLE_estimation,A=var$A)[index,] %>% filter(A==1)
  sample0 <- data.frame(scenario_p$TMLE_estimation,A=var$A)[index,] %>% filter(A==0)
  diff <- colMeans(sample1)-colMeans(sample0)
  ks.score <- list()
  plot(density(sample1[,3]-diff[3]))
  lines(density(sample0[,3]))
  for (i in 1:length(diff)){
    ks.score <- append(ks.score, ks.test(sample1[,i]-diff[i],sample0[,i])$p.value)
  }
  ks.score <- plyr::ldply(ks.score, function(x) x)
  return(ks.score)
}
KS_hte <- ComputeKS(scenario_3_p)
KS_nohte <- ComputeKS(scenario_3_nohte_p)

#diff <- colMeans(sample1)-colMeans(sample0)
#combined_data <- rbind(sample1-diff,sample0)
# ks.score <- list()
# for (i in 1:1999){
#   index.1 <- sample(seq(1,20000,1),10000)
#   combine.1 <- combined_data[index.1]
#   combine.0 <- combined_data[!seq(1,20000,1)%in%index.1]
#   ks.score <- append(ks.score, ks.test(combine.1,combine.0)$p.value)
# }
# ks.score <- plyr::ldply(ks.score, function(x) x)
# sum(ks.score<=one.ks$p.value)/2000



#---------------------Coverage ---------------------------
get_kernal_estimation <- function(onesample,tmle.estimation,tmle.EIC){
  var <-plyr::ldply(scenario_p$scenario_var, function(x) x)
  W_names <- grep('W',names(var),value = T)
  controls <- onesample[,grep("W",names(onesample),value = T)]
  rtn <- list()
  n1<-seq(0,max(controls[,j]),max(controls[,j])/50)
  for (ii in 1:50){
    I <- which(controls[,j] > n1[ii] & controls[,j] < n1[ii+1])
    psi <- apply(tmle.estimation[I,],2,mean)
    sd_EIC <- apply(tmle.EIC[I,], 2, sd)
    sd_EIC <- zoo::na.locf(sd_EIC, option = 'nocb')
    out <- data_frame(as.numeric(psi), se=as.numeric(sd_EIC/sqrt(length(I))))
  }
  
  return(rtn)
}

library("cluster")
library("factoextra")


onesample <-  plyr::ldply(scenario_3_t, function(x) x)
one.est <- get_kernal_estimation(onesample,scenario_3_p$TMLE_estimation,scenario_3_p$EIC_estimation)
t.est <- get_kernal_estimation(onesample,scenario_3_p$True.eff,scenario_3_p$EIC_estimation)

combined.est <- list()
for (i in 1:length(one.est)){
  names(one.est[[i]]) <- c(1:length(one.est[[i]]))
  combined.est[[i]] <- plyr::ldply(one.est[[i]],
                                   function(x) data.frame(x[,1] %>% t(),x[,2] %>% t()))
}
for (i in 1:length(t.est)){
  names(t.est[[i]]) <- c(1:length(t.est[[i]]))
  t.est[[i]] <- plyr::ldply(t.est[[i]],
                            function(x) data.frame(x[,1] %>% t()))
}



#---------------------HTE by TIME---------------------
plotResults <- function(X1=controls[,1],X2=controls[,2],ptps,tau.hat,tau=NULL,start,end){
  # # Interpolate
  tau.hat <-  tau.hat[,ptps]
  if(!is.null(tau)){
    tau <- tau[,ptps]
    RMSE <- sqrt(mean((tau.hat - tau)^2))
    RMSE <- paste0(" / RMSE: ",round(RMSE,4))
  }else{
    RMSE <- ""
  }
  RMSE <- ""
  
  if(is.null(start)){
    start <- round(min(zi))
    end <- round(max(zi))
  }
  
}

controls <-plyr::ldply(scenario_3_p$scenario_var, function(x) x)
W_names <- grep('W',names(controls),value = T)
nW <- length(W_names)
index <- sample(seq(1:100000),10000)
fitdat <- data.frame(Y=raw_scenario_3$V12,controls)[index,]
tau.hat <- data.frame(Y=scenario_3_p$True.eff$`12`,controls)[index,]
RMSE <- sqrt(mean((tau.hat$Y - fitdat$Y)^2))
RMSE <- paste0(" / RMSE: ",round(RMSE,4))

#Noc
controls <-plyr::ldply(scenario_3_noc_p$scenario_var, function(x) x)
fitdat <- data.frame(Y=raw_scenario_3_noc$V12,controls)[index,]
tau.hat <- data.frame(Y=scenario_3_noc_p$True.eff$`12`,controls)[index,]
RMSE <- sqrt(mean((tau.hat$Y - fitdat$Y)^2))
RMSE <- paste0(" / RMSE: ",round(RMSE,4))

#True
controls <-plyr::ldply(scenario_3_noc_p$scenario_var, function(x) x)
tau.hat <- data.frame(Y=scenario_3_noc_p$True.eff$`12`,controls)[index,]
fitdat <- tau.hat
RMSE <- sqrt(mean((tau.hat$Y - fitdat$Y)^2))
RMSE <- paste0(" / RMSE: ",round(RMSE,4))


g.SL.1  <- SuperLearner(Y=fitdat$Y, X=fitdat[,W_names], SL.library=c("SL.earth","SL.glm.interaction","SL.glmnet"), family="gaussian")
sample1 <- data.frame(W1=c(1:10000))
pred.dat <- var
pred.dat <- as.data.frame(t(rep(0,nW))) %>% dplyr::slice(rep(1:n(), each = 10000))
names(pred.dat) <- W_names
eval(parse(text=paste0("pred.dat$W",1," <- t(rbind(rep(0,5000),rep(1,5000)))[1:10000]")))
eval(parse(text=paste0("pred.dat$W",2," <- rep(seq(0,1,0.0002),2)[1:10000]")))
g.pred <- predict(g.SL.1, pred.dat)

pred.func <- function(x1,x2){
  dat <- data.frame(W1=x1,W2=x2,W3=0,W4=0,W5=0,W6=0,W7=0,W8=0,W9=0,W10=0)
  g.pred <- predict(g.SL.1, dat)
  return(g.pred$pred)
}
W1 <- as.numeric(t(rbind(rep(0,100),rep(1,100)))[1:200])
W2 <- as.numeric(rep(seq(0,1,0.005),2)[1:200])
prediction <-  outer(W1, W2, pred.func)
rownames(prediction) <- W1
colnames(prediction) <- W2
mtrx.melt <- reshape::melt(prediction, id.vars = c("W1", "W2"), measure.vars = "Z")

g = ggplot(mtrx.melt, aes(x = X1, y = X2, z = value )) +
  geom_tile(aes(fill = value)) +
  geom_contour(color = "white", alpha = 0.5)+
  ggtitle(paste0("Time: ",12, RMSE))+
  scale_fill_distiller(palette="YlGnBu", na.value="white",limits = c(-0.16,-0.01)) +
  xlab("X\U2081") +
  ylab("X\U2082") +
  guides(fill = guide_colorbar(title =sprintf("ITE")))+
  theme_hc()+
  theme(legend.position = "right")+
  theme(text=element_text(size=8, family="Arial Unicode MS"),
        plot.title = element_text(size = 8, family="Arial Bold"))

ggsave("CH3_time_12.png", plot = g, device = 'png', path = paste0(getwd(),'/GFX'),
       scale = 1, width = 70, height = 50, units = c("mm"), dpi = 300)
ggsave("CH3_time_12_noc.png", plot = g, device = 'png', path = paste0(getwd(),'/GFX'),
       scale = 1, width = 70, height = 50, units = c("mm"), dpi = 300)
ggsave("CH3_time_12_true.png", plot = g, device = 'png', path = paste0(getwd(),'/GFX'),
       scale = 1, width = 70, height = 50, units = c("mm"), dpi = 300)



plotTime <- function(X1,X2,tau.hat,tau=NULL,
                     start=-.5,end=0,time_start,time_end){
  p <- list()
  p.names <- ""
  for (i in time_start:time_end){
    p[[i]] <-  plotResults(X1,X2,ptps=i,
                           tau.hat=tau.hat,
                           tau=true.eff, start=start,end=end)
    p[[i]] <- p[[i]]+
      theme(axis.title.x = element_blank())+
      theme(axis.title.y = element_blank())
    p.names <- c(p.names,paste0("p[[",i,"]]"))
  }
  p.names <- p.names[1:time_end+1]
  p <- eval(parse(text= paste0("ggarrange(",paste(p.names,collapse = ","),
                               ",nrow=3,ncol=4,common.legend = TRUE,legend = 'bottom')")))
  
  # eval( parse(text= paste("subplot(",paste(p.names,collapse = ","),
  #                         paste(",shareX = T,shareY = T,nrows = 3)"))))
}
p <- plotTime(X1=controls[,1],X2=controls[,2],tau.hat=prediction,
              tau = true.eff,start=-.30,end=0,
              time_start=8,time_end=12)
ggsave("CH3_modeeate_confound(5000-10000).png", plot = p, device = 'png', path = paste0(getwd(),'/GFX'),
       scale = 1, width = 200, height = 200, units = c("mm"), dpi = 300)





#---------------------kernal grouping for HTE---------------------
subgroupCI <- function(scenario_p = scenario_3_p,convery.idx,ii){
  sample.ii <- scenario_p$TMLE_estimation[convery.idx==ii,]
  psi.i <-  colMeans(sample.ii)
  sd_EIC <- apply(scenario_p$EIC_estimation[convery.idx==ii,], 2, sd)
  sd_EIC <- zoo::na.locf(sd_EIC, option = 'nocb')
  upper_CI <- psi.i + 1.95 * sd_EIC/sqrt(nrow(sample.ii))
  lower_CI <- psi.i - 1.95 * sd_EIC/sqrt(nrow(sample.ii))
  out <- data.frame(psi.i, sd_EIC, upper_CI, lower_CI)
  return(out)
}

compute_group <- function(scenario_p = scenario_3_p,step){
  controls <-plyr::ldply(scenario_p$scenario_var, function(x) x)
  W_names <- grep('W',names(var),value = T)
  selection <- controls$W2
  n1<-seq(0,max(selection),max(selection)/step)
  convery.idx <- selection
  maEnbk1 <- list();
  for (ii in 1:step){
    I <- which(selection > n1[ii] & selection < n1[ii+1])
    maEnbk1[[ii]] <- I
    convery.idx[I] <- ii
  }
  convery.idx[which(convery.idx<1)] = round(convery.idx[which(convery.idx<1)])
  convery.idx <- ifelse(convery.idx==0,2,convery.idx)
  return(convery.idx)
}

compute_effect <- function(scenario_p,step){
  convery.idx <- compute_group(scenario_p,step)
  Psi_group <- list()
  for (ii in 1:step){
    true.i <- data.frame(y=colMeans(scenario_p$True.eff[convery.idx==ii,]), x=c(1:12))
    sample.i <- subgroupCI(scenario_p,convery.idx,ii)
    sample.i$y <- sample.i$psi.i
    sample.i$x <- c(1:12)
    sample.i$true <- true.i$y
    sample.i$group <- ii
    Psi_group[[ii]] <- sample.i
  }
  return(plyr::ldply(Psi_group, function(x) x))
}

kernal_process <- function(scenario,scenario_t,scenario_tmle,scenario_EIC,ratDiv,step){
  effect_by_sample <- list()
  for (i in 1:50){
    cat(i)
    scenario_p_v2 <- get_result(scenario[i],scenario_t[i], scenario_tmle[i], scenario_EIC[i],ratDiv=ratDiv)
    effect_by_sample[[i]] <- compute_effect(scenario_p = scenario_p_v2,step=step)
  }
  
  effect_by_sample_p <- plyr::ldply(effect_by_sample, function(x) x)
  effect_grouped_psi <- aggregate(effect_by_sample_p[,c("psi.i","sd_EIC","upper_CI","lower_CI","true")],
                                  by=list(time = effect_by_sample_p$x,group=effect_by_sample_p$group), mean)
  effect_grouped_sd <- aggregate(effect_by_sample_p[,c("psi.i","sd_EIC","upper_CI","lower_CI","true")],
                                 by=list(time = effect_by_sample_p$x,group=effect_by_sample_p$group), sd)
  effect_grouped_psi$sd_EIC <- effect_grouped_sd$sd_EIC
  effect_grouped_psi$upper_CI <- effect_grouped_psi$psi.i+1.96*effect_grouped_sd$sd_EIC/sqrt(50)
  effect_grouped_psi$lower_CI <- effect_grouped_psi$psi.i-1.96*effect_grouped_sd$sd_EIC/sqrt(50)
  return(effect_grouped_psi)
}

kernal_3_p <- kernal_process(scenario_3,scenario_3_t, scenario_3_tmle, scenario_3_EIC,ratDiv=115,step = 15)
kernal_3_noc_p <- kernal_process(scenario_3_noc,scenario_3_noc_t, scenario_3_noc_tmle, scenario_3_noc_EIC,ratDiv=115,step = 15)

kernal_3_p <- data.frame(group = "Event rate 30%",kernal_3_p)
scenario_3_p$Summary
kernal_3_noc_p <- data.frame(group = "No confounding",kernal_3_noc_p)
scenario_3_noc_p$Summary


kernal_plot <- rbind(kernal_3_p,kernal_3_noc_p)
kernal_plot$group.1 <- kernal_plot$group.1/15
#write.csv(kernal_plot, file = paste0(getwd(),'/result/kernal_plot.csv') , row.names = FALSE)
kernal_plot$Bias <- abs((kernal_plot$psi.i-kernal_plot$true)/kernal_plot$psi.i)
kernal_plot$Bias[kernal_plot$group=="Event rate 0.1%"] <- NA
kernal_plot$Bias[kernal_plot$group=="Event rate 1%"] <- NA



# kernal_plot <- kernal_plot[,c("group","time","group.1","psi.i","Bias","upper_CI",'lower_CI')]
# kernal_plot <- reshape2::melt(kernal_plot, id.vars=c("group","time","group.1","lower_CI","upper_CI"), value.name = "value")
# kernal_plot$variable = gsub("psi.i","(a) Estimated CATE",kernal_plot$variable)
# kernal_plot$variable = gsub("Bias","(b) Estimation %Bias",kernal_plot$variable)


#
g=ggplot(kernal_plot) +
  geom_line(data=subset(kernal_plot,time == '12'), aes(x = group.1, y = psi.i,lty = group), size=.7)+
  geom_ribbon(data=subset(kernal_plot,time == '12'), aes(ymin=lower_CI,ymax=upper_CI,x=group.1,group=group),
              alpha=0.2)+
  #facet_wrap(~variable,ncol=2,scales="free", strip.position = "bottom")+
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

library(scales)
g=ggplot(kernal_plot) +
  geom_line(data=subset(kernal_plot,time == '12'), aes(x = group.1, y = Bias, lty = group), size=.7)+
  #facet_wrap(~variable,ncol=2,scales="free", strip.position = "bottom")+
  scale_linetype_manual("",values=c("twodash", "solid","longdash","dotdash","dotted"))+
  theme_hc()+xlab("X\U2082") +ylab("%Bias") +
  guides(lty = guide_legend(nrow = 2, byrow = TRUE))+
  scale_y_continuous(labels=percent)+
  theme(legend.spacing.y = unit(0, "mm"),
        panel.border = element_rect(colour = "white", fill=NA),
        aspect.ratio = .8, axis.text = element_text(colour = 1, size = 8),
        strip.background = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = NA,fill = "#FFFFFF"))
ggsave("CH3_kernal_bias.png", plot = g, device = 'png', path = paste0(getwd(),'/GFX'),
       scale = 1, width = 150, height = 100, units = c("mm"), dpi = 250)









#A histgram check
# hist(true.0[,ptps], col=rgb(0,0,1,1/4),breaks=100,lty="blank")
# hist(Y.hat.0[,ptps], col=rgb(1,0,0,1/4),breaks=100,add=T,lty="blank")
# hist(true.1[,ptps], col=rgb(0,0,1,1/4),breaks=100,lty="blank")
# hist(Y.hat.1[,ptps], col=rgb(1,0,0,1/4),breaks=100,add=T,lty="blank")
#
# hist.data <- data.frame(x= c(Y.hat.0[,ptps],true.0[,ptps]),
#                         y = c(rep("est",length(Y.hat.0[,ptps])),rep("true",length(true.0[,ptps]))))
# g = ggplot(data=hist.data,aes(x=x)) + labs(x = "Survival Probablity",y = "Count")+theme_hc()+
#   geom_histogram(data=subset(hist.data,y == 'est'),aes(fill = "#D70026"),  alpha = 0.2,bins=50, show.legend = T)+
#   geom_histogram(data=subset(hist.data,y == 'true'),aes(fill = "#828282"), alpha = 0.2,bins=50, show.legend = T)+
#   scale_fill_manual("Colour",values= c( "#D70026","#828282"),labels = c("Estimation", "True"))+
#   theme(legend.spacing.y = unit(0, "mm"),
#         panel.border = element_rect(colour = "white", fill=NA),
#         aspect.ratio = 1, axis.text = element_text(colour = 1, size = 10),
#         legend.background = element_blank(),
#         legend.box.background = element_rect(colour = NA,fill = "#FFFFFF"),
#         legend.position=c(.15, .9))



