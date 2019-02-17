#' TMLE update 
#' @export
fit_survtmle <- function(T.tilde, Delta, A, W_df) {
  t_0 <- max(T.tilde)
  fit <- survtmle(
    ftime = T.tilde,
    ftype = Delta,
    trt = A,
    t0 = max(T.tilde),
    adjustVars = W_df,
    SL.trt =  c("SL.mean", "SL.glm", "SL.gam"),
    SL.ftime = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
    SL.ctime = c("SL.mean", "SL.glm", "SL.gam", "SL.earth"),
    method = "hazard",
    returnIC = TRUE,
    verbose = FALSE
  )
  # extract cumulative incidence at each timepoint
  tpfit <- timepoints(fit, times = seq_len(t_0))
  len_groups <- as.numeric(
    unique(lapply(lapply(tpfit, FUN = `[[`, "est"), FUN = length))
  )
  names_groups <- unique(
    lapply(lapply(tpfit, FUN = `[[`, "est"), FUN = rownames)
  )[[1]]
  est_only <- t(matrix(unlist(
    lapply(tpfit, FUN = `[[`, "est")
  ), ncol = len_groups, byrow = TRUE))
  est_only <- as.data.frame(est_only)
  rownames(est_only) <- names_groups
  colnames(est_only) <- paste0("t", seq_len(ncol(est_only)))
  
  s_0 <- 1 - as.numeric(est_only[1, ])
  s_1 <- 1 - as.numeric(est_only[2, ])
  return(data.frame(time = 1:t_0, s_0 = s_0, s_1 = s_1))
}
