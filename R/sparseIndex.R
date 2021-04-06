#' Compute a sparse index on the results of a sparse GSVD
#'
#' @param X the data matrix
#' @param res.gsvd the result of a sparse GSVD of X
#' @param tol a tolerance parameter indicating when a small value should be considered equal to 0
#'
#' @return
#' @export
#'
#' @examples
sparseIndex <- function(X, res.gsvd, tol = 1e-16) {
  R <- length(res.gsvd$d)
  I <- nrow(X)
  J <- ncol(X)
  rdsLeft <- res.gsvd$rdsLeft
  rdsRight <- res.gsvd$rdsRight
  # Compute the fit part of the index
  d0_full <- svd(X, nu = 0, nv = 0)$d
  d0 <- d0_full[1:R]
  dsparse <- res.gsvd$d
  r1 <- sum(dsparse^2) / sum(d0^2)

  # Compute the sparsity part of the index
  U <- res.gsvd$U
  V <- res.gsvd$V
  n0inU <- sum(abs(U) <= tol)
  n0inV <- sum(abs(V) <= tol)
  radiusIndexLeftG <- gmean(rdsLeft / sqrt(I))
  radiusIndexRightG <- gmean(rdsRight / sqrt(J))
  radiusIndexLeftA <- mean(rdsLeft / sqrt(I))
  radiusIndexRightA <- mean(rdsRight / sqrt(J))

  r2 <- n0inU / (I * R)
  r3 <- n0inV / (J * R)
  r4 <- (n0inU + n0inV) / ((I + J) * R)
  # Combine
  SI1 <- gmean(c(r1, mean(c(r2, r3))))
  SI2 <- gmean(c(r1, r4))
  SI3 <- r1 * mean(c(r2, r3))
  SI4 <- r1 * r4
  SI5 <- gmean(c(r1, 1 - radiusIndexLeftG, 1 -  radiusIndexRightG))
  SI6 <- prod(c(r1, 1 - radiusIndexLeftG, 1 - radiusIndexRightG))
  SI7 <- mean(c(r1, 1 - radiusIndexLeftA, 1 - radiusIndexRightA))

  return(list(SI1 = SI1, SI2 = SI2, SI3 = SI3, SI4 = SI4,
              SI5 = SI5, SI6 = SI6, SI7 = SI7,
              r1 = r1, r2 = r2, r3 = r3, r4 = r4,
              n0inU = n0inU, n0inV = n0inV,
              rdsLeft = rdsRight, rdsLeft = rdsRight,
              radiusIndexLeftG = radiusIndexLeftG,
              radiusIndexRightG = radiusIndexRightG,
              radiusIndexLeftA = radiusIndexLeftA,
              radiusIndexRightA = radiusIndexRightA))
}


#' Geometric mean
#'
#' @param x vector of positive values
#' @param na.rm should missing values be removed (default to TRUE)
#'
#' @return the geometric mean of x
#' @export
#'
#' @examples
#' gmean(c(0.5, 0.5))
gmean <- function(x, na.rm = TRUE) {
  if (any(abs(x) < 2*.Machine$double.eps)) return(0)
  if (any(x < 0)) return(0)
  exp(mean(log(x), na.rm = na.rm))
}
