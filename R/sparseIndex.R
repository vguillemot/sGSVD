#' Compute sparse indices on the results of a sparse GSVD
#'
#' @param res.sgsvd The result of a sparse Generalized Singular Value Decomposition of a data-matrix, usually obtained with sGSVD::sparseGSVD or one of the companion functions in the SPAFAC package ;
#' @param singularValues The singular values of the original data matrix (used to compute the fit of the sparse analysis)
#' @param correction The type of correction for proportion of explained variance (r1), e.g., "gsvd" (no correction), "mca", "mfa"
#' @param tol a tolerance parameter indicating when a small value should be considered equal to 0
#'
#' @return various sparsity indices and their components.
#' @export
#'
#' @examples
#' set.seed(2045)
#' X <- matrix(rnorm(20), 4, 5)
#' res.svd <- svd(X)
#' res.sgsvd <- sparseGSVD(X, k = 2L)
#' sparseIndex(res.sgsvd, res.svd$d)
sparseIndex <- function(res.sgsvd, singularValues, correction = "gsvd", tol = 1e-16) {
  R <- length(res.sgsvd$d)
  U <- res.sgsvd$U
  V <- res.sgsvd$V
  U.sq <- U^2
  V.sq <- V^2
  if (is.null(res.sgsvd$grpLeft)){
    ctrLeft <- U.sq
  }else{
    ctrLeft <- apply(U.sq, 2, function(x) tapply(x, res.sgsvd$grpLeft, FUN = sum))
  }
  if (is.null(res.sgsvd$grpRight)){
    ctrRight <- V.sq
  }else{
    ctrRight <- apply(V.sq, 2, function(x) tapply(x, res.sgsvd$grpRight, FUN = sum))
  }
  I <- NROW(ctrLeft)
  J <- NROW(ctrRight)
  rdsLeft <- res.sgsvd$rdsLeft
  rdsRight <- res.sgsvd$rdsRight

  # Compute the fit part of the index
  # d0 <- singularValues[1:R]
  # dsparse <- res.sgsvd$d
  # r1 <- sum(dsparse^2) / sum(d0^2)
  r1 <- compute.fit(singularValues, res.sgsvd$d, J, correction = correction)

  # Compute the sparsity part of the index

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
  SI1left <- gmean(c(r1, r2))
  SI1right <- gmean(c(r1, r3))
  SI2 <- gmean(c(r1, r4))
  SI2left <- gmean(c(r1, r2))
  SI2right <- gmean(c(r1, r3))
  SI3 <- r1 * mean(c(r2, r3))
  SI3left <- r1 * r2
  SI3right <- r1 * r3
  SI4 <- r1 * r4
  SI4left <- r1 * r2
  SI4right <- r1 * r3
  SI5 <- gmean(c(r1, 1 - radiusIndexLeftG, 1 -  radiusIndexRightG))
  SI5left <- gmean(c(r1, 1 - radiusIndexLeftG))
  SI5right <- gmean(c(r1, 1 -  radiusIndexRightG))
  SI6 <- prod(c(r1, 1 - radiusIndexLeftG, 1 - radiusIndexRightG))
  SI6left <- prod(c(r1, 1 - radiusIndexLeftG))
  SI6right <- prod(c(r1, 1 - radiusIndexRightG))
  SI7 <- mean(c(r1, 1 - radiusIndexLeftA, 1 - radiusIndexRightA))
  SI7left <- mean(c(r1, 1 - radiusIndexLeftA))
  SI7right <- mean(c(r1, 1 - radiusIndexRightA))

  return(list(
    SI1 = SI1, SI1left = SI1left, SI1right = SI1right,
    SI2 = SI2, SI2left = SI2left, SI2right = SI2right,
    SI3 = SI3, SI3left = SI3left, SI3right = SI3right,
    SI4 = SI4, SI4left = SI4left, SI4right = SI4right,
    SI5 = SI5, SI5left = SI5left, SI5right = SI5right,
    SI6 = SI6, SI6left = SI6left, SI6right = SI6right,
    SI7 = SI7, SI7left = SI7left, SI7right = SI7right,
    r1 = r1, r2 = r2, r3 = r3, r4 = r4,
    n0inU = n0inU, n0inV = n0inV,
    rdsLeft = rdsLeft, rdsRight = rdsRight,
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

#' Compute r1 (the explained variance) for sparsity indices
#'
#' @param d a vector of singular values
#' @param pseudo.d a vector of the pseudo singular values (from sparseGSVD)
#' @param J the number of variables (for MCA correction)
#' @param correction the type of correction, e.g., "gsvd" (no correction), "mca", "mfa"
#'
#' @return the corrected eigenvalues and tau
#' @export
#'
#' @examples
#'
compute.fit <- function(d, pseudo.d, J, correction = "gsvd") {
  if (correction == "mca"){
    lambda <- (J/(J-1)*(d^2-(1/J)))^2
    pseudo.lambda <- (J/(J-1)*(pseudo.d^2-(1/J)))^2
    lambda[lambda < (1/J)] = 0
    pseudo.lambda[pseudo.lambda < (1/J)] = 0
    r1 <- sum(pseudo.lambda)/sum(lambda)
  }else if (correction == "mfa"){
    stop("MFA correction is not available yet.")
  }else{
    r1 <- sum(pseudo.d^2) / sum(d^2)
  }
  return(r1)
}
