
survtmle_single_t <- function(dat,
                              tk,
                              dW = rep(1, nrow(dat)),
                              T.cutoff = NULL,
                              SL.ftime = c("SL.glm","SL.mean","SL.step", "SL.earth"),
                              SL.ctime = c("SL.glm","SL.mean"),
                              SL.trt = c("SL.glm","SL.mean","SL.step", "SL.earth")){
  # ===================================================================================
  # preparation
  # ===================================================================================
  dat <- dat
  dW <- dW
  n.data <- nrow(dat)
  W_names <- grep('W', colnames(dat), value = TRUE)

  T.uniq <- unique(sort(dat$T.tilde))

  # create function inputs
  ftime <- dat$T.tilde
  ftype <- dat$Delta

  if(all(dW == 0)) {
    trt <- 1 - dat$A # when dW is all zero, flip observed A
  }else if(all(dW == 1)){
    trt <- dat$A
  }else{
    stop('not implemented!')
  }

  adjustVars <- as.data.frame(dat[,W_names])
  # ====================================================================================
  # compute values for all time points
  # ====================================================================================
  fit_max_time <-survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
                           t0 = tk,
                           SL.ftime = SL.ftime,
                           SL.ctime = SL.ctime,
                           SL.trt = SL.trt,
                           # glm.ftime = paste(c('trt', W_names), collapse = ' + '),
                           # glm.trt = paste(W_names, collapse = ' + '),
                           method="hazard", returnModels = TRUE,
                           verbose = FALSE)
  # 7.8min
  # allTimes <- timepoints(object = fit_max_time, times = T.uniq, returnModels = FALSE)
  est <- 1 - fit_max_time$est['1 1',]
  var <- fit_max_time$var['1 1', '1 1']

  return(list(est = est,
              var = var,
              meanIC = fit_max_time$meanIC,
              ic = fit_max_time$ic))
}

