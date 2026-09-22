#' Constrained eigen-value decomposition of a symmetric matrix with an additionnal positive constraint
#'
#' @param X a symmetric square (data) matrix;
#' @param k the desired rank of the singular decomposition;
#' @param init How to initialize the algorithm, Default: 'svd'
#' @param rds The radiuses (radii?) (>0) of the L1 or LG constraint; one for each dimension;
#' @param seed a random seed for result reproducibility; if NULL (the default), no random seed will be used
#' @param grp vector describing the groups; default to one group per row
#' @param orthogonality whether the orthogonality constraint is applied on the "loadings" (default)
#' @param OrthSpace matrix defining the orthogonal space, Default: NULL
#' @param projPriority the order in which the projections are applied, Default: 'sparsity'
#' @param itermaxALS the maximum number of ALS iterations, Default: 1000
#' @param itermaxPOCS the maximum number of POCS iterations, Default: 1000
#' @param epsALS precision for ALS, Default: 1e-10
#' @param epsPOCS precision for POCS, Default: 1e-10
#' @return Pseudo-eigen vectors and values
#' @examples
#' U <- matrix(rnorm(20), 5, 4)
#' sparseEIGENpos(U %*% t(U))
#' @author Vincent Guillemot
#' @export
sparseEIGENpos <- function(
  X, k = 2L,
  init = "svd", seed = NULL,
  rds = rep(1, k),
  grp = NULL,
  orthogonality = "loadings",
  OrthSpace = NULL,
  projPriority = "sparsity",
  itermaxALS = 1000, itermaxPOCS = 1000,
  epsALS = 1e-10, epsPOCS = 1e-10) {

  # Test that the arguments are valid
  garb <- runTestsEIGEN(X, k, init, seed,
                        rds, grp,
                        orthogonality, OrthSpace,
                        projPriority)

  I <- nrow(X)

  # Build initialization matrices either with SVD (prefered method)
  # or randomly
  res.init <- initializeEIGEN(X = X, I = I, k = k,
                              init = init, seed = seed)
  U0 <- res.init$U0
  # Build projection based on the arguments
  proj <- makeComposedPositiveProjection(projPriority = projPriority, grp = grp)

  if (is.null(OrthSpace)) OrthSpace <- matrix(0, I, 1)
  U <- matrix(0, I, k)

  iter <- matrix(NA, k, 2,
                 dimnames = list(paste0("Dim. ", 1:k),
                                 c("Total", "ALS")))
  lambda <- rep(NA, k)

  for (r in 1:k) {
    ## Power Iteration with orth projection
    res.powit <- powerIteration(
      X = X,                 # original matrix
      init = U0[,r],  # initialization vectors
      proj = proj,
      rds = rds[r],
      grp = grp,
      OrthSpace = OrthSpace,
      itermaxALS = 1000, itermaxPOCS = 1000,
      epsALS = 1e-10, epsPOCS = 1e-10)

    U[, r] <- res.powit$u

    if (orthogonality == "loadings") {
      OrthSpace <- U
    }

    iter[r,] <- c(res.powit$iterTOTAL, res.powit$iterALS)
    lambda[r] <- res.powit$lambda
  }

  oD <- order(lambda, decreasing = TRUE)
  res <- list(values = lambda[oD], vectors = U[, oD], iter = iter)
  return(res)
}

#' Build the composed positive projection function for a given priority order
#'
#' @param projPriority the order in which the projections are applied (cannot be 'orth')
#' @param grp vector describing the groups, or NULL for a plain L1L2 constraint
#'
#' @return a projection function suitable for \code{\link{powerIteration}}
#' @noRd
makeComposedPositiveProjection <- function(projPriority, grp) {
  if (projPriority == "orth") {
    stop("Cannot have orthogonal priority with a positive constraint... yet!")
  } else {
    if (is.null(grp)) return(projOrth_then_projPos_then_projL1L2)
    return(projOrth_then_projPos_then_projLGL2)
  }
}

