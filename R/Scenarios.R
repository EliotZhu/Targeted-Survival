#' get the simulation scenario
#' @export
get.data <- function(iti=1234,samplesize=1000, conmode="scenario 3",ratDiv=1,confoundlevel=4,p=1){
  mypoly <- function(x){
    out <- data.frame("X1"=x)
    for (i in 1:p){
      eval(parse(text= paste0("out$X",i," <- x^i")))
    }
    return(apply(out,1,sum))
  }


  D <- DAG.empty()
  D <- D +
    node("W1", distr ="runif", min = 0, max = 2)+
    node("W2", distr ="runif", min = 0, max = 2)+
    node("W3", distr ="rbinom", prob = .5,size=1)+
    node("W4", distr ="rbinom", prob = .5,size=1)
  if(conmode == "scenario 3"){
    for (i in 5:20){
      D <- D +eval(parse(text= paste0("node('W",i,"', distr ='rbinom', prob = .1,size=1)")))
    }
    D <- D+ node("odds",distr = "rconst", const = exp(0.25 + confoundlevel*(W1 - W2 + W3 - W4)))+
        #node("A", distr = "rbinom", size = 1, prob = ifelse(W1==1,1,0)) +
        node("A", distr = "rbinom", size = 1, prob = odds / (1 + odds)) +
        node("rate",distr = "rconst", const = ((mypoly(W1)+W2+mypoly(W3)+W4)*A+
                                                mypoly(W1)+W2+mypoly(W3)+W4)/ratDiv)+
        node("Cweib", distr = "rweibull", shape = 1+W5/5, scale = 50)+
        node("Trexp", distr = "rexp", rate = rate) +
        node("T", distr = "rconst", const = round(Trexp/11)) +
        node("C", distr = "rconst", const = round(Cweib/11)) 
    wnames <- grep('W',names(D),value = T)
    true_surv <- function(x,tgrid,A){
      x <- as.matrix(x,nrow=1)
      rate <- as.numeric((mypoly(x[1])+cos(x[2])+mypoly(x[3])+x[4])*A+mypoly(x[1])+cos(x[2])+mypoly(x[3])+x[4])
      s_diff_true <-    1 - pexp(seq(0,tgrid,1)*11, rate = rate/ratDiv)
      return(s_diff_true)
    }
  }else if(conmode == "scenario s"){
    D <- D+node("A", distr = "rbinom", size = 1, prob = .15 + .5 * as.numeric(W2 > .75)) +
      node("rate",distr = "rconst", const = 1 + .7 * W2^2 - .8 * A)+
      node("Cweib", distr = "rweibull", shape = 1 + .5 * W2, scale = 75)+
      node("Trexp", distr = "rexp", rate = rate) +
      node("T", distr = "rconst", const = round(Trexp*2)) +
      node("C", distr = "rconst", const = round(Cweib*2)) 
    wnames <- c('W1','W2')
    true_surv <- function(x,tgrid,A){
      x <- as.matrix(x,nrow=1)
      rate <- as.numeric(1 + .7 * x[2]^2 - .8 * A)
      s_diff_true <-    1 - pexp(seq(0,tgrid,1)/2, rate = rate)
      return(s_diff_true)
    }
    
  }else if(conmode == "no"){
    for (i in 5:20){
      D <- D +eval(parse(text= paste0("node('W",i,"', distr ='rbinom', prob = .1,size=1)")))
    }
    
    D <- D+ node("odds",distr = "rconst", const = 1)+
      #node("A", distr = "rbinom", size = 1, prob = ifelse(W1==1,1,0)) +
      node("A", distr = "rbinom", size = 1, prob = odds / (1 + odds)) +
      node("rate",distr = "rconst", const = ((W1^2+W2+W3^2+W4)*A+W1^2+W2+W3^2+W4)/ratDiv)+
      node("Cweib", distr = "rweibull", shape = 1+W5/5, scale = 50)+
      node("Trexp", distr = "rexp", rate = rate) +
      node("T", distr = "rconst", const = round(Trexp/10)) +
      node("C", distr = "rconst", const = round(Cweib/10)) 
    wnames <- grep('W',names(D),value = T)
    true_surv <- function(x,tgrid,A){
      x <- as.matrix(x,nrow=1)
      rate <- as.numeric((x[1]^2+x[2]+x[3]^2+x[4])*A+x[1]^2+x[2]+x[3]^2+x[4])
      s_diff_true <-    1 - pexp(seq(0,tgrid,1)*10, rate = rate/ratDiv)
      return(s_diff_true)
    }
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

  }

  D <- D+node("T.tilde", distr = "rconst", const = ifelse(T <= C , T, C)) +
         node("Delta", distr = "rconst", const = ifelse(T <= C , 1, 0))
  setD <- set.DAG(D)

  dat <- sim(setD,n=samplesize,rndseed= iti)
  dat2 <- sim(setD,n=10,rndseed= iti)
  data_out <- dat[,names(dat) %in% c("ID",wnames,"A","T.tilde","Delta" )]
  data_out2 <- dat[,names(dat2) %in% c("ID",wnames,"A","T.tilde","Delta" )]
  

  return(list(dat = as.data.frame(data_out), wnames = wnames,
              true_surv = true_surv, 
              dat2 = as.data.frame(data_out2)))
}





