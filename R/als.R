#' Alternative least squares method helper for the cGSVD
#'
#' @param X data matrix on which the ALS is performed
#' @param initLeft initial left vector
#' @param initRight initial right vector
#' @param projLeft projection function applied to the left vector at each iteration
#' @param projRight projection function applied to the right vector at each iteration
#' @param rdsLeft radius of the constraint applied to the left vector
#' @param rdsRight radius of the constraint applied to the right vector
#' @param grpLeft vector describing the groups for the left vector, Default: NULL
#' @param grpRight vector describing the groups for the right vector, Default: NULL
#' @param OrthSpaceLeft matrix defining the orthogonal space for the left vector
#' @param OrthSpaceRight matrix defining the orthogonal space for the right vector
#' @param itermaxALS the maximum number of ALS iterations, Default: 1000
#' @param itermaxPOCS the maximum number of POCS iterations, Default: 1000
#' @param epsALS precision for ALS, Default: 1e-10
#' @param epsPOCS precision for POCS, Default: 1e-10
#'
#' @return A list with the singular value \code{d}, the left and right pseudo-singular
#' vectors \code{u} and \code{v}, the number of ALS iterations \code{iterALS}, and
#' the total number of POCS iterations \code{iterTOTAL}.
#' @export
#'
#' @examples
#' X <- matrix(rnorm(20), 5, 4)
#' als(X, initLeft = svd(X)$u[, 1], initRight = svd(X)$v[, 1],
#'     projLeft = projOrth_then_projL1L2, projRight = projOrth_then_projL1L2,
#'     rdsLeft = 1, rdsRight = 1,
#'     OrthSpaceLeft = matrix(0, 5, 1), OrthSpaceRight = matrix(0, 4, 1),
#'     itermaxPOCS = 100, epsPOCS = 1e-8)

als <- function(X, initLeft, initRight, projLeft, projRight,
                rdsLeft, rdsRight, grpLeft = NULL, grpRight = NULL,
                OrthSpaceLeft, OrthSpaceRight,
                itermaxALS = 1000, itermaxPOCS = 1000,
                epsALS = 1e-10, epsPOCS = 1e-10) {
  uold <- unew <- initLeft
  vold <- vnew <- initRight

  iterTOTAL <- 0
  for (iter in 1:itermaxALS) {
    res.projRight <- projRight(vec = crossprod(X, uold), rds = rdsRight, grp = grpRight, OrthSpace = OrthSpaceRight, itermax = itermaxPOCS, eps = epsPOCS)
    vnew <- res.projRight$x
    Xvnew <- X %*% vnew
    res.projLeft <- projLeft(vec = Xvnew, rds = rdsLeft, grp = grpLeft, OrthSpace = OrthSpaceLeft, itermax = itermaxPOCS, eps = epsPOCS)
    unew <- res.projLeft$x
    iterTOTAL <- iterTOTAL + res.projLeft$k + res.projRight$k
    if ( normL2(vnew - vold) < epsALS && normL2(unew - uold) < epsALS ) break
    vold <- vnew
    uold <- unew
  }

  d <- drop(crossprod(unew, Xvnew))
  return(list(d = d, u = unew, v = vnew, iterALS = iter, iterTOTAL = iterTOTAL))
}

#' Power iteration with a constrained projection at each step
#'
#' @param X data matrix on which the power iteration is performed
#' @param init initial vector
#' @param proj projection function applied to the vector at each iteration
#' @param rds radius of the constraint applied to the vector
#' @param grp vector describing the groups, Default: NULL
#' @param OrthSpace matrix defining the orthogonal space
#' @param itermaxALS the maximum number of ALS iterations, Default: 1000
#' @param itermaxPOCS the maximum number of POCS iterations, Default: 1000
#' @param epsALS precision for ALS, Default: 1e-10
#' @param epsPOCS precision for POCS, Default: 1e-10
#'
#' @return A list with the eigenvalue \code{lambda}, the pseudo-eigenvector
#' \code{u}, the number of ALS iterations \code{iterALS}, and the total
#' number of POCS iterations \code{iterTOTAL}.
#' @noRd
powerIteration <- function(
    X, init, proj, rds, grp = NULL, OrthSpace,
    itermaxALS = 1000, itermaxPOCS = 1000,
    epsALS = 1e-10, epsPOCS = 1e-10) {
  uold <- unew <- init

  iterTOTAL <- 0
  for (iter in 1:itermaxALS) {
    res.proj <- proj(vec = X %*% uold, rds = rds, grp = grp,
                     OrthSpace = OrthSpace, itermax = itermaxPOCS,
                     eps = epsPOCS)
    unew <- res.proj$x
    iterTOTAL <- iterTOTAL + res.proj$k
    if ( normL2(unew - uold) < epsALS ) break
    uold <- unew
  }

  lambda <- drop(crossprod(unew, X %*% unew))
  return(list(lambda = lambda, u = unew, iterALS = iter, iterTOTAL = iterTOTAL))
}


