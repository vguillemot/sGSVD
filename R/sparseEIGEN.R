#' Constrained eigen-value decomposition of a symmetric matrix
#'
#' @param X a symmetric square (data) matrix;
#' @param k the desired rank of the singular decomposition;
#' @param init How to initialize the algorithm
#' @param rds The radiuses (radii?) (>0) of the L1 or LG constraint; one for each dimension;
#' @param seed
#' @param grp
#' @param orthogonality
#' @param OrthSpace
#' @param projPriority
#' @param itermaxALS
#' @param itermaxPOCS
#' @param epsALS
#' @param epsPOCS
#' $L_1$ ball for each left vector
#' @return Pseudo-eigen vectors and values
#' @examples
#' U <- matrix(rnorm(20), 5, 4)
#' sparseEIGEN(U %*% t(U))
#' @author Vincent Guillemot
#' @export
sparseEIGEN <- function(X, k = 2L,
                 init = NULL, seed = NULL,
                 rds = rep(1, k),
                 grp = NULL,
                 orthogonality = "loadings",
                 OrthSpace = NULL,
                 projPriority = "orth",
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
  proj <- makeComposedProjection(projPriority = projPriority, grp = grp)

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

runTestsEIGEN <- function(X, k, init, seed,
                         rds, grp,
                         orthogonality, OrthSpace,
                         projPriority) {

  ##### Test X ####
  if (nrow(X)==1 & ncol(X)==1)
    stop("You are attempting a gsGSVD of a scalar.")

  if (nrow(X) != ncol(X))
    stop("X should be a square matrix.")

  if (!isSymmetric(unname(X)))
    stop("X should be symmetric.")

  if (any(is.na(X)))
    stop("X should not contain missing values")

  ##### Test k ####
  if (!is.integer(k)) stop("k should be an integer.")
  if (k <= 1) stop("k should be > 1.")

  ##### Test initialization ####
  if (is.null(init)) {
    stop("init should be either svd or rand.")
  }
  if (! init %in% c("svd", "rand"))
    stop("init should be either svd or rand.")

  return(NULL)
}

initializeEIGEN <- function(X, I, k, init, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  if (init == "svd") {
    svdx <- svd(X, nu=k, nv=k)
    U0 <- svdx$u
  } else if (init == "rand") {
    U0 <- 1/(I-1) * mvrnorm(n = I, mu = rep(0,k),
                            Sigma = diag(k), empirical = TRUE)
  } else {
    stop("Unkown error, please contact support!")
  }

  return(list(U0 = U0))
}
