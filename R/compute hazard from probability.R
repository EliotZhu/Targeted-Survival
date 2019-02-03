# ===================================================================================
# compute hazard
# ===================================================================================

finite.differences <- function(x, y) {
  if (length(x) != length(y)) {
    stop('x and y vectors must have equal length')
  }
  n <- length(x)
  # Initialize a vector of length n to enter the derivative approximations
  fdx <- vector(length = n)
  # Iterate through the values using the forward differencing method
  for (i in 2:n) {
    fdx[i-1] <- (y[i-1] - y[i]) / (x[i-1] - x[i])
  }
  # For the last value, since we are unable to perform the forward differencing method
  # as only the first n values are known, we use the backward differencing approach
  # instead. Note this will essentially give the same value as the last iteration
  # in the forward differencing method, but it is used as an approximation as we
  # don't have any more information
  fdx[n] <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])
  return(fdx)
}

h.hat.t_full_1 <- matrix(NA, nrow = n.data, ncol = ncol(Qn.A1.t_full_1))
h.hat.t_full_0 <- matrix(NA, nrow = n.data, ncol = ncol(Qn.A1.t_full_0))
for (i in 1:n.data) {
  H <- -log(Qn.A1.t_full_1[i,])
  h.hat.t_full_1[i,] <- finite.differences(rep(seq(1,ncol(Qn.A1.t_1)),interval1),H)
  h.hat.t_full_1[i,] <- ifelse(is.na(h.hat.t_full_1[i,]),0,h.hat.t_full_1[i,])
}
for (i in 1:n.data) {
  H <- -log(Qn.A1.t_full_0[i,])
  h.hat.t_full_0[i,] <- finite.differences(rep(seq(1,ncol(Qn.A1.t_0)),interval1),H)
  h.hat.t_full_0[i,] <- ifelse(is.na(h.hat.t_full_0[i,]),0,h.hat.t_full_0[i,])
}
