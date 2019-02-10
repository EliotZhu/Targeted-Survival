library(here,usethis)
library(dplyr,abind)
library(tidyverse)
library(survival,simsurv)
library(survminer)
library(simcausal)

get.data <- function(iti,samplesize, conmode, endtime=50,ratDiv){
  D <- DAG.empty()
  D <- D +
    node("W1", distr ="rbinom", prob = .5,size=1)+
    node("W2", distr ="runif", min = 0, max = 1)+
    node("W3", distr ="rbinom", prob = .5,size=1)+
    node("W4", distr ="runif", min = 0, max = 1)+
    node("W5", distr ="rbinom", prob = .5,size=1)+
    node("W6", distr ="runif", min = 0, max = 1)+
    node("W7", distr ="rbinom", prob = .5,size=1)+
    node("W8", distr ="runif", min = 0, max = 1)+
    node("W9", distr ="runif", min = 0, max = 1)+
    node("W10", distr ="runif", min = 0, max = 1)

  if (conmode == "scenario 1"){
    D <- D+ node("odds",distr = "rconst", const = 0.2+log(1.2)*W1+log(1.2)*W2)+
      node("A", distr = "rbinom", size = 1, prob = odds / (1 + odds)) +
      node("rate",distr = "rconst", const = (W5+W6+W7+W8+W1+W2+W3+W4)+(W1+W2+W3+W4)*A)
    wnames <- c('W1','W2','W3','W4','W5','W6','W7','W8','W9','W10')
  }else if(conmode == "scenario 2"){
    D <- D+ node("odds",distr = "rconst", const = 0.2+log(1.2)*W1+log(1.2)*W2+log(2)*W3+log(2)*W4)+
      node("A", distr = "rbinom", size = 1, prob = odds / (1 + odds)) +
      node("rate",distr = "rconst", const = (W5+W6+W7+W8+W1+W2+W3+W4)+(W1+W2+W3+W4)*A)
    wnames <- c('W1','W2','W3','W4','W5','W6','W7','W8','W9','W10')
  }else if(conmode == "scenario 3"){
    D <- D+ node("odds",distr = "rconst", const = 0.2+W1+W2+2*W3+2*W4+W5+W6+W7+W8)+
      node("A", distr = "rbinom", size = 1, prob = odds / (1 + odds)) +
      node("rate",distr = "rconst", const = (W5+W6+W7+W8+W1+W2+W3+W4)+(W1+W2+W3+W4)*A)
    wnames <- c('W1','W2','W3','W4','W5','W6','W7','W8','W9','W10')
  }else if(conmode == "no"){
    D <- D+ node("odds",distr = "rconst", const = 1)+
      node("A", distr = "rbinom", size = 1, prob = .5) +
      node("rate",distr = "rconst", const = (W5+W6+W7+W8+W1+W2+W3+W4)+(W1+W2+W3+W4)*A)
    wnames <- c('W1','W2','W3','W4','W5','W6','W7','W8','W9','W10')
  }else if(conmode == "scenario 3.1"){
    D <- D+ node("odds",distr = "rconst", const = 0.2+log(1.2)*W1+log(1.2)*W2+log(2)*W3+log(2)*W4+log(1.2)*W5+log(1.2)*W6+log(2)*W7+log(2)*W8)+
      node("A", distr = "rbinom", size = 1, prob = odds / (1 + odds)) +
      node("rate",distr = "rconst", const = (W5+W6+W7+W8+W1+W2+W3+W4)+A)
    wnames <- c('W1','W2','W3','W4','W5','W6','W7','W8','W9','W10')
  }else if(conmode == "scenario HTE 1"){
    for (i in 10:30){
      D <- D +eval(parse(text= paste0("node('W",i,"', distr ='runif', min = 0, max = 1)")))
    }
    D <- D+ node("odds",distr = "rconst", const = 0.2+W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12+W13+W14+W15)+
      node("A", distr = "rbinom", size = 1, prob = odds / (1 + odds))
    wnames <- grep('W',names(D),value = T)
    D <- D+ eval(parse(text= paste0("node('rate',distr = 'rconst', const =",paste(wnames, collapse ="+"),
                                    "+(W10+W11+W12+W13+W14+W15)*A+A)")))
    wnames <- grep('W',names(D),value = T)

  }else if(conmode == "scenario HTE 2"){
    for (i in 10:30){
      D <- D +eval(parse(text= paste0("node('W",i,"', distr ='rbinom', prob = .5,size=1)")))
    }
    D <- D+ node("odds",distr = "rconst", const = 0.2+W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12+W13+W14+W15)+
      node("A", distr = "rbinom", size = 1, prob = odds / (1 + odds))
    wnames <- grep('W',names(D),value = T)
    D <- D+ eval(parse(text= paste0("node('rate',distr = 'rconst', const =",paste(wnames, collapse ="+"),
                                    "+(W10+W11+W12+W13+W14+W15)*A*log(1.2)+A)")))

  }else if(conmode == "scenario HTE 3"){
    for (i in 10:60){
      D <- D +eval(parse(text= paste0("node('W",i,"', distr ='runif', min = 0, max = 1)")))
    }
    D <- D+ node("odds",distr = "rconst", const = 0.2+W1+W2+W3+W4+W5+W6+W7+W8+W9+W10+W11+W12+W13+W14+W15)+
      node("A", distr = "rbinom", size = 1, prob = odds / (1 + odds))
    wnames <- grep('W',names(D),value = T)
    D <- D+ eval(parse(text= paste0("node('rate',distr = 'rconst', const =",paste(wnames, collapse ="+"),
                                    "+(W10+W11+W12+W13+W14+W15)*A*log(1.2)+A)")))
  }

  D <- D+
    node("rate1", distr = "rconst", const = rate/ratDiv) +
    node("Trexp", distr = "rexp", rate = rate1) +
    node("Cweib", distr = "rweibull", shape = .8 - .01*(rnorm(1)), scale = 20) +
    node("T", distr = "rconst", const = round(Trexp,0)) +
    node("C", distr = "rconst", const = round(Cweib,0)) +
    #node("C", distr = "rconst", const = T+1) +
    node("T.tilde", distr = "rconst", const = ifelse(T <= C , T, C)) +
    node("Delta", distr = "rconst", const = ifelse(T <= C , 1, 0))
  setD <- set.DAG(D)

  dat <- sim(setD,n=samplesize,rndseed= iti)

  data_out <- dat[,names(dat) %in% c("ID",wnames,"A","T.tilde","Delta" )]
  
  dat2 <- sim(setD,n=20000,rndseed= 12345)
  
  true.func <- function(x,points,ratDiv=ratDiv,A){
    x <- as.matrix(x,nrow=1)
    rate <- as.numeric(x[5]+x[6]+x[7]+x[8]+x[1]+x[2]+x[3]+x[4]+A*(x[1]+x[2]+x[3]+x[4]))
    s_diff_true <-  exp(-rate/ratDiv*seq(0,points,1))
    return(s_diff_true)
  }
  

  return(list(dat = as.data.frame(data_out), wnames = wnames,
              true_surv = true.func, 
              dat2 = dat2))
}





