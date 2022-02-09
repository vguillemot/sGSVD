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
sparseIndex <- function(res.sgsvd, singularValues, correction = "gsvd", tol = 1e-10) {
  R <- length(res.sgsvd$d)
  singularValues <- singularValues[1:R]
  U <- res.sgsvd$U
  V <- res.sgsvd$V
  U.sq <- U^2
  V.sq <- V^2
  if (is.null(res.sgsvd$grpLeft)) {
    ctrLeft <- U.sq
  } else {
    ctrLeft <- apply(U.sq, 2, function(x) tapply(x, res.sgsvd$grpLeft, FUN = sum))
  }
  if (is.null(res.sgsvd$grpRight)) {
    ctrRight <- V.sq
  } else {
    ctrRight <- apply(V.sq, 2, function(x) tapply(x, res.sgsvd$grpRight, FUN = sum))
  }
  I <- NROW(ctrLeft)
  J <- NROW(ctrRight)
  rdsLeft <- res.sgsvd$rdsLeft
  rdsRight <- res.sgsvd$rdsRight

  # Compute the fit part of the index
  # d0 <- singularValues[1:R]
  # dsparse <- res.sgsvd$d
  # r1 <- cumsum(dsparse^2) / cumsum(d0^2)
  r1 <- compute.fit(singularValues, res.sgsvd$d, J, correction = correction)

  # Compute the sparsity part of the index

  n0inU <- cumsum(colSums(ctrLeft <= tol))
  n0inV <- cumsum(colSums(ctrRight <= tol))
  radiusIndexLeftG <- cumgmean(rdsLeft / sqrt(I))
  radiusIndexRightG <- cumgmean(rdsRight / sqrt(J))
  radiusIndexLeftA <- cummean(rdsLeft / sqrt(I))
  radiusIndexRightA <- cummean(rdsRight / sqrt(J))

  r2 <- n0inU / (I * (1:R))
  r3 <- n0inV / (J * (1:R))
  r4 <- (n0inU + n0inV) / ((I + J) * (1:R))
  # Combine
  SI <- r1 * r4
  SIleft <- r1 * r2
  SIright <- r1 * r3

  return(list(
    SI = SI, SIleft = SIleft, SIright = SIright,
    r1 = r1, r2 = r2, r3 = r3, r4 = r4,
    n0inU = n0inU, n0inV = n0inV,
    rdsLeft = rdsLeft, rdsRight = rdsRight,
    radiusIndexLeftG = radiusIndexLeftG,
    radiusIndexRightG = radiusIndexRightG,
    radiusIndexLeftA = radiusIndexLeftA,
    radiusIndexRightA = radiusIndexRightA))
}

#' @rdname sparseIndex
#' @export

sparseIndexEigen <- function(res.sgevd, eigenValues, correction = "gevd", tol = 1e-10) {
  R <- length(res.sgevd$values)
  eigenValues <- eigenValues[1:R]
  U <- res.sgevd$vectors
  U.sq <- U^2
  if (is.null(res.sgevd$grp)) {
    ctr <- U.sq
  } else {
    ctr <- apply(U.sq, 2, function(x) tapply(x, res.sgevd$grp, FUN = sum))
  }

  I <- NROW(ctr)
  rds <- res.sgevd$rds

  # Compute the fit part of the index
  fitRatio <- compute.fit.eigen(eigenValues, res.sgevd$values, I, correction = correction)
  # Compute the sparsity part of the index
  n0 <- cumsum(colSums(ctr <= tol))
  zeroRatio <- n0 / (I * (1:R))
  # Combine
  SI <- fitRatio * zeroRatio

  return(list(
    SI = SI,
    fitRatio = fitRatio,
    zeroRatio = zeroRatio,
    n0 = n0,
    r1 = fitRatio, r2 = NA, r3 = NA, r4 = zeroRatio, rds = rds))
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
  if (correction == "mca") {
    lambda <- (J / (J - 1) * (d ^ 2 - (1 / J))) ^ 2
    pseudo.lambda <- (J / (J - 1) * (pseudo.d ^ 2 - (1 / J))) ^ 2
    lambda[lambda < (1 / J)] <- 0
    # pseudo.lambda[pseudo.lambda < (1/J)] = 0
    r1 <- cumsum(pseudo.lambda) / cumsum(lambda)
  } else if (correction == "mfa") {
    stop("MFA correction is not available yet.")
  } else {
    r1 <- cumsum(pseudo.d ^ 2) / cumsum(d ^ 2)
  }
  return(r1)
}

#' @rdname compute.fit
#' @export

compute.fit.eigen <- function(ev, pseudo.ev, I, correction = "gevd") {

  if (correction != "gevd") {
    warning("Eigenvalue corrections are not yet supported")
  }

  r1 <- cumsum(pseudo.ev) / cumsum(ev)

  return(r1)
}
#' Cumulative arithmetic mean
#'
#' @param x a vector of numeric values
#'
#' @return the cumulative arithmetic mean
#' @export
#'
#' @examples
cummean <- function(x, na.rm = FALSE) {
  if (na.rm) x <- na.omit(x)
  return(cumsum(x) / seq_along(x))
}

#' Cumulative geometric mean
#'
#' @param x vector of positive values
#' @param na.rm should missing values be removed (default to TRUE)
#'
#' @return the cumulative geometric mean of x
#' @export
#'
#' @examples
#' gmean(c(0.5, 0.5))
cumgmean <- function(x, na.rm = TRUE) {
  if (any(abs(x) < 2*.Machine$double.eps)) return(0)
  if (any(x < 0)) return(0)
  exp(cummean(log(x), na.rm = na.rm))
}

