---
title: Targeted estimation of heterogeneous treatment effect in observational survival
  analysis
output:
  html_document:
    df_print: paged
  md_document:
    variant: markdown_github
---
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.1016/j.jbi.2020.103474.svg)](https://doi.org/10.1016/j.jbi.2020.103474)


This paper proposed a three-stage modular design for estimating the treatment effect heterogeneity in observational survival analysis. The method provides monotonic survival curves with adjustment for selection and censoring bias. We automate the identification of features contributing to the effect heterogeneity. We avoid the ad-hoc subgroup selection by non-parametrically estimating the conditional treatment effect. We provide evidence that oral anticoagulants confer protection against stroke and death on newly diagnosed non-valvular atrial fibrillation patients.



```r
require(dplyr)
require(MOSS) #devtools::install_github('wilsoncai1992/MOSS')
require(survival)
#require(simcausal) #if you want complicate simulation, please install from local directory install.packages("~/simcausal_0.5.5.tar", repos = NULL)
require(abind)
require(tidyverse)
```

Generate some sample 


```r
n <- 100
W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.5))
A <- rbinom(n, 1, 0.5)
EventTime <- rgeom(n,plogis(-4 + W$W1 * W$W2 - A)) + 1
CensorTime <- rgeom(n, plogis(-6 + W$W1)) + 1
T.tilde <- pmin(EventTime, CensorTime)
Delta <- as.numeric(T.tilde == EventTime)
df = data.frame(A = A, T.tilde = T.tilde, Delta = Delta, W1 = W$W1, W2=W$W2)
df$ID <- seq.int(nrow(df))
max_time = 30
head(df)
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["A"],"name":[1],"type":["int"],"align":["right"]},{"label":["T.tilde"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["Delta"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["W1"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["W2"],"name":[5],"type":["int"],"align":["right"]},{"label":["ID"],"name":[6],"type":["int"],"align":["right"]}],"data":[{"1":"1","2":"70","3":"1","4":"0.9409179","5":"0","6":"1","_rn_":"1"},{"1":"0","2":"22","3":"1","4":"0.2297090","5":"1","6":"2","_rn_":"2"},{"1":"0","2":"123","3":"0","4":"0.7850615","5":"0","6":"3","_rn_":"3"},{"1":"0","2":"20","3":"1","4":"0.9759315","5":"1","6":"4","_rn_":"4"},{"1":"0","2":"219","3":"1","4":"0.2648601","5":"0","6":"5","_rn_":"5"},{"1":"0","2":"7","3":"1","4":"0.3117311","5":"0","6":"6","_rn_":"6"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

Run the method proposed in the step 1 and 2 of paper 


```r
  df <- df[df$T.tilde<= max_time & df$T.tilde>0,]
  df <- df[complete.cases(df),]


  adjustVars <- grep('W', colnames(df), value = T)
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
```

```
## <survival_curve>
##   Public:
##     ci: function (A, T_tilde, Delta, density_failure, density_censor, 
##     clone: function (deep = FALSE) 
##     create_ggplot_df: function (W = NULL) 
##     display: function (type, W = NULL) 
##     hazard: 0.0451819736967532 0.0436533823808446 0.0489955579854037 ...
##     hazard_to_pdf: function () 
##     hazard_to_survival: function () 
##     initialize: function (t, hazard = NULL, survival = NULL, pdf = NULL) 
##     n: function () 
##     pdf: NULL
##     pdf_to_hazard: function () 
##     pdf_to_survival: function () 
##     survival: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1  ...
##     survival_to_hazard: function () 
##     survival_to_pdf: function () 
##     t: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ...
```

```r
  sl_fit$density_failure_0$hazard_to_survival()
```

```
## <survival_curve>
##   Public:
##     ci: function (A, T_tilde, Delta, density_failure, density_censor, 
##     clone: function (deep = FALSE) 
##     create_ggplot_df: function (W = NULL) 
##     display: function (type, W = NULL) 
##     hazard: 0.0475411182334166 0.0455994921160187 0.0523617507987143 ...
##     hazard_to_pdf: function () 
##     hazard_to_survival: function () 
##     initialize: function (t, hazard = NULL, survival = NULL, pdf = NULL) 
##     n: function () 
##     pdf: NULL
##     pdf_to_hazard: function () 
##     pdf_to_survival: function () 
##     survival: 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1  ...
##     survival_to_hazard: function () 
##     survival_to_pdf: function () 
##     t: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 ...
```

```r
  sl_fit$density_failure_1$t <- k_grid
  sl_fit$density_failure_0$t <- k_grid

  sl_density_failure_1_marginal <- sl_fit$density_failure_1$clone(deep = TRUE)
  sl_density_failure_0_marginal <- sl_fit$density_failure_0$clone(deep = TRUE)
  sl_density_failure_1_marginal$survival <- matrix(colMeans(sl_density_failure_1_marginal$survival), nrow = 1)
  sl_density_failure_0_marginal$survival <- matrix(colMeans(sl_density_failure_0_marginal$survival), nrow = 1)


  out <- list(sl_fit_1 = sl_fit$density_failure_1$survival,
              sl_fit_0 = sl_fit$density_failure_0$survival,
              SL_diff = sl_fit$density_failure_1$survival-sl_fit$density_failure_0$survival)
```

Individual difference in survival probabilities


```r
head(out$SL_diff)
```

```
##      [,1]        [,2]        [,3]        [,4]        [,5]        [,6]
## [1,]    0 0.002359145 0.004497207 0.006436149 0.008203044 0.009827627
## [2,]    0 0.001946110 0.003716567 0.005328653 0.006804153 0.008167389
## [3,]    0 0.003366193 0.006387993 0.009100609 0.011545264 0.013765537
## [4,]    0 0.002637547 0.005021719 0.007177871 0.009136799 0.010931966
## [5,]    0 0.002384616 0.004545251 0.006504171 0.008288781 0.009929158
## [6,]    0 0.003811482 0.007218212 0.010262176 0.012991508 0.015456271
##             [,7]       [,8]       [,9]      [,10]      [,11]      [,12]
## [1,] 0.011337991 0.01275748 0.01410449 0.01539323 0.01663610 0.01784411
## [2,] 0.009441638 0.01064649 0.01179761 0.01290743 0.01398703 0.01504674
## [3,] 0.015801386 0.01768509 0.01944117 0.02108767 0.02263917 0.02410741
## [4,] 0.012594720 0.01415090 0.01562064 0.01701929 0.01835995 0.01965397
## [5,] 0.011453691 0.01288595 0.01424450 0.01554366 0.01679587 0.01801222
## [6,] 0.017701800 0.01976430 0.02167093 0.02344123 0.02509057 0.02663071
##           [,13]      [,14]      [,15]      [,16]      [,17]      [,18]
## [1,] 0.01902703 0.02019475 0.02135683 0.02252124 0.02369351 0.02487641
## [2,] 0.01609610 0.01714536 0.01820513 0.01928558 0.02039590 0.02154456
## [3,] 0.02550122 0.02682793 0.02809236 0.02929475 0.03042919 0.03148225
## [4,] 0.02091102 0.02214046 0.02335081 0.02454812 0.02573500 0.02690988
## [5,] 0.01920244 0.02037640 0.02154357 0.02271176 0.02388624 0.02506940
## [6,] 0.02806973 0.02941344 0.03066411 0.03181823 0.03286469 0.03378321
##           [,19]      [,20]      [,21]      [,22]      [,23]     [,24]     [,25]
## [1,] 0.02606765 0.02725571 0.02841098 0.02947377 0.15092907 0.2407442 0.2460838
## [2,] 0.02273787 0.02397710 0.02525113 0.02652517 0.05514450 0.1033213 0.2346520
## [3,] 0.03242920 0.03322828 0.03381156 0.03407675 0.12312549 0.1967486 0.1958203
## [4,] 0.02806421 0.02917764 0.03020855 0.03108242 0.08011228 0.1709877 0.2191278
## [5,] 0.02625838 0.02744085 0.02858609 0.02963296 0.06488300 0.1318424 0.2215471
## [6,] 0.03454033 0.03508374 0.03533446 0.03518209 0.07687517 0.1477393 0.1716609
##          [,26]     [,27]      [,28]       [,29]
## [1,] 0.2267264 0.1994208 0.03762036 0.005657391
## [2,] 0.2459229 0.2265690 0.19915492 0.167877704
## [3,] 0.1744332 0.1486427 0.03174058 0.003834908
## [4,] 0.2083863 0.1838137 0.15577120 0.036182515
## [5,] 0.2206802 0.1980732 0.17041200 0.123817172
## [6,] 0.1550720 0.1312250 0.10652947 0.016409125
```

```r
hist(out$SL_diff)
```

<img src="README_files/figure-html/unnamed-chunk-4-1.png" width="672" />
The average on the ITEs will be the TMLE adjusted ATE. 

```r
#ATE
mean(out$SL_diff)
```

```
## [1] 0.04917922
```

Please refer to the paper for BART variable importance measure and kernal grouping to get the CATE for subgroups.  




