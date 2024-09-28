#' Calculate Weighted Sparsity Indices from Sparse GSVD Results
#'
#' This function computes weighted sparsity indices based on the results of a
#' sparse Generalized Singular Value Decomposition (sGSVD). It combines measures
#' of sparsity and fit to provide a comprehensive index.
#'
#' @param res.sgsvd A list containing the results of a sparse Generalized
#'   Singular Value Decomposition (sGSVD) applied to a data matrix.
#'   Typically generated using `sGSVD::sparseGSVD` or similar functions in
#'   the sGSVD package.
#' @param singularValues A numeric vector of the singular values from the
#'   original data matrix, used to assess the fit of the sparse GSVD.
#'   These are usually obtained from a standard SVD or GSVD.
#' @param correction A string indicating the type of correction to apply when
#'   calculating the proportion of explained variance (`r1`). Options include
#'   `"gsvd"` (no correction), `"mca"`, or `"mfa"`. Default is `"gsvd"`.
#' @param weight A numeric value (between 0 and 1) specifying the weight
#'   applied to the fit ratio when combining it with the sparsity measure.
#'   The complement `(1 - weight)` is applied to the sparsity ratio. Default
#'   is `0.5`, meaning equal weight for fit and sparsity.
#' @param tol A numeric tolerance value to define when elements should be
#'   considered zero. Any value below this threshold is treated as zero.
#'   Default is `1e-10`.
#'
#' @return A list with three components:
#' \describe{
#'   \item{SI}{A numeric vector containing the weighted sparsity indices
#'     combining fit and sparsity.}
#'   \item{SIleft}{A numeric vector containing the weighted sparsity indices
#'     for the left singular vectors (rows).}
#'   \item{SIright}{A numeric vector containing the weighted sparsity indices
#'     for the right singular vectors (columns).}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(2045)
#' X <- matrix(rnorm(20), 4, 5)
#' res.svd <- svd(X)
#' res.sgsvd <- sparseGSVD(X, k = 2L)
#' weightedSparsityIndex(
#'     res.sgsvd = res.sgsvd,
#'     singularValues = res.svd$d,
#'     correction = "gsvd",
#'     weight = 0.4) # (slightly) favors sparsity over fit
weightedSparsityIndex <- function(
    res.sgsvd,
    singularValues,
    correction = "gsvd",
    weight = 0.5,
    tol = 1e-10) {
  R <- length(res.sgsvd$d)
  singularValues <- singularValues[1:R]
  U <- res.sgsvd$u
  V <- res.sgsvd$v
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
  fitratio <- compute.fit(singularValues, res.sgsvd$d, J, correction = correction)

  # Compute the sparsity part of the index
  n0inU <- cumsum(colSums(ctrLeft <= tol))
  n0inV <- cumsum(colSums(ctrRight <= tol))
  zeroRatioU <- n0inU / (I * (1:R))
  zeroRatioV <- n0inV / (J * (1:R))
  zeroRatio <- (n0inU + n0inV) / ((I + J) * (1:R))

  # Combine
  SI <- fitratio^(weight) * zeroRatio^(1-weight)
  SIleft <- fitratio^(weight) * zeroRatioU^(1-weight)
  SIright <- fitratio^(weight) * zeroRatioV^(1-weight)

  return(list(SI = SI, SIleft = SIleft, SIright = SIright))
}
