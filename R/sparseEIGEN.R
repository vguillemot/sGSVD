

#' Constrained eigen-value decomposition of a symmetric matrix
#'
#' @param X a symmetric square (data) matrix;
#' @param k the desired rank of the singular decomposition;
#' @param init How to initialize the algorithm.
#' svd (default) = initializes with the SVD,
#' eig or eigen = initializes with the EVD,
#' rand = initializes with a random vector.
#' @param rds The radiuses (radii?) (>0) of the L1 or LG constraint; one for each dimension;
#' @param seed a random seed for result reproducibility; if NULL (the default), no random seed will be used
#' @param grp vector describing the groups; default to one group per row
#' @param orthogonality whether the orthogonality constraint is applied on the "loadings" (default)
#' @param OrthSpace matrix defining the orthogonal space, Default: NULL
#' @param projPriority the order in which the projections are applied, Default: 'orth'
#' @param compute_sparsity_index whether the sparsity index should be computed, Default: TRUE
#' @param correction4SI correction for the explained variance for sparsity indices, e.g., "gevd" (no correction)
#' @param itermaxALS the maximum number of ALS iterations, Default: 1000
#' @param itermaxPOCS the maximum number of POCS iterations, Default: 1000
#' @param epsALS precision for ALS, Default: 1e-10
#' @param epsPOCS precision for POCS, Default: 1e-10
#' @param tol.si tolerance for the computation of the Sparse Index, set by default to .Machine$double.eps
#' @return Pseudo-eigen vectors and values
#' @examples
#' U <- matrix(rnorm(20), 5, 4)
#' sparseEIGEN(U %*% t(U))
#' @author Vincent Guillemot
#' @export
sparseEIGEN <- function(X, k = 2L,
                 init = "svd", seed = NULL,
                 rds = rep(1, k),
                 grp = NULL,
                 orthogonality = "loadings",
                 OrthSpace = NULL,
                 projPriority = "orth",
                 compute_sparsity_index = TRUE,
                 correction4SI = "gevd",
                 itermaxALS = 1000, itermaxPOCS = 1000,
                 epsALS = 1e-10, epsPOCS = 1e-10,
                 tol.si = .Machine$double.eps) {

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

  # Compute Sparsity Index
  res$rds <- rds
  res$grp <- grp

  if (compute_sparsity_index) {
    res.SI <- sparseIndexEigen(
      res.sgevd = res,
      eigenValues = eigen(X, only.values = TRUE)$values,
      correction = correction4SI,
      tol = tol.si)
    res$SI <- res.SI
  }

  return(res)
}

#' Validate the arguments of \code{\link{sparseEIGEN}} and \code{\link{sparseEIGENpos}}
#'
#' @param X a symmetric square (data) matrix
#' @param k the desired rank of the decomposition
#' @param init how to initialize the algorithm
#' @param seed a random seed for result reproducibility
#' @param rds the radiuses of the L1 or LG constraint
#' @param grp vector describing the groups
#' @param orthogonality whether the orthogonality constraint is applied on the loadings
#' @param OrthSpace matrix defining the orthogonal space
#' @param projPriority the order in which the projections are applied
#'
#' @return NULL, invisibly; called for its side effect of raising an error on invalid input
#' @noRd
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
    stop("init should be either svd, eigen, eig, rand, or a matrix.")
  }

  if (is.character(init)) {
    if (!init %in% c("svd", "eigen", "eig", "rand"))
      stop("init should be either svd, eigen, eig, or rand.")
  }

  if (is.matrix(init)) {
    if (any(dim(init) != c(nrow(X), k))) {
      stop("Wrong dimensions for initalization matrix!")
    }
  }

  return(NULL)
}

#' Build the initial pseudo-eigenvectors for \code{\link{sparseEIGEN}} and \code{\link{sparseEIGENpos}}
#'
#' @param X a symmetric square (data) matrix
#' @param I number of rows/columns of X
#' @param k the desired rank of the decomposition
#' @param init how to initialize the algorithm
#' @param seed a random seed for result reproducibility, Default: NULL
#'
#' @return A list with the initial vectors \code{U0}
#' @noRd
initializeEIGEN <- function(X, I, k, init, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  if (is.character(init)) {
    if (init == "svd" | init == "eigen" | init == "eig" ) {
      svdx <- svd(X, nu=k, nv=0)
      U0 <- svdx$u
    } else if (init == "rand") {
      U0 <- 1/(I-1) * mvrnorm(n = I, mu = rep(0,k),
                              Sigma = diag(k), empirical = TRUE)
    }
  } else if (is.matrix(init)) {
    U0 <- init
  } else {
    stop("Wrong initialization parameters. Should be eig, eigen, svd, or rand.")
  }

  return(list(U0 = U0))
}
