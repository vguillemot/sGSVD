#' Alternative least squares method helper for the cGSVD
#'
#' @param initLeft
#' @param initRight
#' @param projLeft
#' @param projRight
#' @param rdsLeft
#' @param rdsRight
#' @param grpLeft
#' @param grpRight
#' @param OrthSpaceLeft
#' @param OrthSpaceRight
#' @param itermaxALS
#' @param itermaxPOCS
#' @param epsALS
#' @param epsPOCS
#' @param X data matrix on which the ALS is performed
#'
#' @return
#' @export
#'
#' @examples
#'

als <- function(X, initLeft, initRight, projLeft, projRight, rdsLeft, rdsRight, grpLeft = NULL, grpRight = NULL, OrthSpaceLeft, OrthSpaceRight, itermaxALS = 1000, itermaxPOCS = 1000, epsALS = 1e-10, epsPOCS = 1e-10) {
  uold <- unew <- initLeft
  vold <- vnew <- initRight

  iterTOTAL <- 0
  for (iter in 1:itermaxALS) {
    res.projRight <- projRight(vec = t(X) %*% uold, rds = rdsRight, grp = grpRight, OrthSpace = OrthSpaceRight, itermax = itermaxPOCS, eps = epsPOCS)
    vnew <- res.projRight$x
    res.projLeft <- projLeft(vec = X %*% vnew, rds = rdsLeft, grp = grpLeft, OrthSpace = OrthSpaceLeft, itermax = itermaxPOCS, eps = epsPOCS)
    unew <- res.projLeft$x
    iterTOTAL <- iterTOTAL + res.projLeft$k + res.projRight$k
    if ( normL2(vnew - vold) < epsALS && normL2(unew - uold) < epsALS ) break
    vold <- vnew
    uold <- unew
  }
  # print(res.projLeft$lambda)
  # print(res.projRight$lambda)

  d <- drop(t(unew) %*% X %*% vnew)
  return(list(d = d, u = unew, v = vnew, iterALS = iter, iterTOTAL = iterTOTAL))
}

powerIteration <- function(X, init, proj, rds, grp = NULL, OrthSpace, itermaxALS = 1000, itermaxPOCS = 1000, epsALS = 1e-10, epsPOCS = 1e-10) {
  uold <- unew <- init

  iterTOTAL <- 0
  for (iter in 1:itermaxALS) {
    res.proj <- proj(vec = X %*% uold, rds = rds, grp = grp, OrthSpace = OrthSpace, itermax = itermaxPOCS, eps = epsPOCS)
    unew <- res.proj$x
    iterTOTAL <- iterTOTAL + res.proj$k + res.proj$k
    if ( normL2(unew - uold) < epsALS ) break
    uold <- unew
  }

  lambda <- drop(t(unew) %*% X %*% unew)
  return(list(lambda = lambda, u = unew, iterALS = iter, iterTOTAL = iterTOTAL))
}


