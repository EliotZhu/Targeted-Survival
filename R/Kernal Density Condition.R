
#Kernal Density (uniform) Estim1te of Expected by time
h2<-1;
n2<-seq(0,max(selection),max(selection)/step);
m2<-matrix(0,nrow= length(n2),ncol=3);
maEnbk1<-matrix(0,nrow= length(n2),ncol=length(tpts));
sAEnbk1<-matrix(0,nrow= length(n2),ncol=length(tpts));
naEnbk1<-matrix(0,nrow= length(n2),ncol=length(tpts));
for (ii in 1:length(tpts)){
  reponse <- score_prop[,ii]
  for (jj in 1:length(n2)){
    I<-which( ((selection-n2[jj])/h2)^2 <= 1)
    maEnbk1[jj,ii]<-mean(reponse[I],na.rm = T);
    sAEnbk1[jj,ii]<-sd(reponse[I],na.rm = T);
    naEnbk1[jj,ii]<-length(I);
  }
}
plot_ly(x=~time, y=~n2, z=~maEnbk1, type='contour',contours = list(start = -.3,end = 0,size = .02))


maEnbk1<-matrix(0,nrow= length(n2),ncol=length(tpts)+1);
for (ii in 1:(length(tpts)+1)){
  temp <- 0
  for (jj in 1:length(n2)){
    maEnbk1[jj,ii]<-cond.func(n2[jj])[ii];
  }
}
plot_ly(x=~time, y=~n2, z=~maEnbk1, type='contour',contours = list(start = -.3,end = 0,size = .02))




library(plotly)

plot_ly(x=~time, y=~n2, z=~maEnbk1, type='contour')
plot_ly(x=~time, y=~n2, z=~maEnbk1-1.96*sAEnbk1/sqrt(naEnbk1), type='surface')
plot_ly(x=~time, y=~n2, z=~maEnbk1+1.96*sAEnbk1/sqrt(naEnbk1), type='surface')
plot_ly(x=~time, y=~n2, z=~naEnbk1, type='contour')

#for a line
predframe <- data.frame(maEnbk1,n1,naEnbk1,sAEnbk1)
ggplot(predframe, aes(n1,maEnbk1))+
  geom_line()+
  geom_ribbon(aes(ymin=maEnbk1-1.96*sAEnbk1,
                  ymax=maEnbk1+1.96*sAEnbk1),alpha=0.3)+
  labs(x = "W",y = "Probablity")+
  theme_hc()





h2<-5;
h1<-10;
if(length(selection)>1){
  n1<-seq(0,max(selection[,1]),max(selection[,1])/step);
  m1<-matrix(0,nrow= length(n1),ncol=3);
  n2<-seq(0,max(selection[,2]),max(selection[,2])/step);
  m2<-matrix(0,nrow= length(n2),ncol=3);
  maEnbk1<-matrix(0,nrow= length(n1),ncol=length(n2));
  sAEnbk1<-matrix(0,nrow= length(n1),ncol=length(n2));
  naEnbk1<-matrix(0,nrow= length(n1),ncol=length(n2));
  for (ii in 1:length(n1)){
    for (jj in 1:length(n2)){
      I<-which(((selection[,1]-n1[ii])/h1)^2 + ((selection[,2]-n2[jj])/h2)^2 <= 1)
      maEnbk1[ii,jj]<-mean(reponse[I]);
      sAEnbk1[ii,jj]<-sd(reponse[I]);
      naEnbk1[ii,jj]<-length(I);
    }
  }
}else{
  n1<-seq(0,max(selection),max(selection)/step);
  m1<-matrix(0,nrow= length(n1),ncol=3);
  maEnbk1<-matrix(0,nrow= length(n1));
  sAEnbk1<-matrix(0,nrow= length(n1));
  naEnbk1<-matrix(0,nrow= length(n1));
  for (ii in 1:length(n1)){
    I <- which(((selection-n1[ii])/h1)^2 <= 1)
    maEnbk1[ii]<-mean(reponse[I]);
    sAEnbk1[ii]<-sd(reponse[I]);
    naEnbk1[ii]<-length(I);
  }
}
library(plotly)
predframe <- data.frame(maEnbk1,n1,naEnbk1,sAEnbk1)
# maEnbk1 <- maEnbk1[,12]
# naEnbk1 <- naEnbk1[,12]
# sAEnbk1 <- sAEnbk1[,12]
# for (jj in 1:length(n2)){
#   maEnbk1[jj]<-cond.func(n2[jj])[ii];
# }
ggplot(predframe, aes(n2,maEnbk1))+
  geom_line()+
  geom_ribbon(aes(ymin=maEnbk1-1.96*sAEnbk1/sqrt(naEnbk1),
                  ymax=maEnbk1+1.96*sAEnbk1/sqrt(naEnbk1)),alpha=0.3)+
  labs(x = "W1",y = "Survival Probablity")+
  theme_hc()




plot_ly(x = ~n2, y=~n1, z=~maEnbk1, type='contour')
plot_ly(x=~n2, y=~n1, z=~maEnbk1-1.96*sAEnbk1/sqrt(naEnbk1), type='surface')
plot_ly(x=~n2, y=~n1, z=~maEnbk1+1.96*sAEnbk1/sqrt(naEnbk1), type='surface')
plot_ly(x=~n2, y=~n1, z=~naEnbk1, type='contour')








#score_prop <- dat_b$A*(dat_b$Y-g1.pred)/ex - (1-dat_b$A)*(dat_b$Y-g0.pred)/(1-ex)+g1.pred-g0.pred;
#score_prop <- dat_b$A*(g1.pred)/ex - (1-dat_b$A)*(g0.pred)/(1-ex)

#Leave-on2-out CV to Choose Bandwidth
I <- seq(1,length(selection[,1]),1); #specify subsample
nnb <- length(I);
Inb <- (1:nnb);
reponse_s<-reponse[I];
selection_1<-selection[,1];
selection_2<-selection[,2];

CV<-matrix(1e6,nrow= length(n1),ncol=length(n2));
for (ii in 1:length(n1)){
  for (jj in 1:length(n2)){
    mkk<- matrix(0,nrow= nnb);
    for (kk in 1:nnb){
      if (ii == 1 & jj == 1){
        Ikk <- which((Inb != kk) & (selection_1 == selection_1[kk]) & (selection_2 == selection_2[kk]))
      }else if(ii == 1 & jj > 1){
        Ikk<- which((Inb != kk) & (selection_1 == selection_1[kk]) & (abs(selection_2-selection_2[kk]) < jj))
      }else if(ii > 1 & jj == 1){
        Ikk<- which((Inb != kk) & (abs(selection_1-selection_1[kk]) < ii) & (selection_2 == selection_2[kk]))
      }else{
        Ikk<- which((Inb != kk) & ((((selection_1 - selection_1[kk])/(ii-1))^2 +
                                      ((selection_2 - selection_2[kk])/(jj-1))^2) <= 1))
      }
      mkk[kk]<-mean(reponse_s[Ikk]);
    }
    CV[ii,jj] <- sum((mkk-reponse_s)^2);
  }
}



#small test
# data <- catheter
# death <- data[,1]  # outcome
# rhc <- data[,2]  # treatment
# sex <- data[,3]
# race <- data[,4]
# income <- data[,5]
# cat1 <- data[,6]
# cat2 <- data[,7]
# ninsclas <- data[,8]
# age <- data[,9]

# n <- length(death)
#
# controls <-  data.frame(W1 = sex, W2 = race, W3 = income, W4 = cat1, W5 = cat2, W6 = ninsclas, W7 = age)
# dat <- data.frame(T.tilde=100,delta=death, A = rhc, controls)
# dat <- data_out[data_out$T.tilde<=24,]
# controls <- dat[,grep("W",names(dat),value = T)]

