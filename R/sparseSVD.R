#' Constrained SVD of a matrix (wrapper of c++ functions).
#'
#' @param X a (data) matrix;
#' @param R the desired rank of the singular decomposition;
#' @param au The radiuses (radii?) (>0) of the
#' $L_1$ ball for each left vector
#' @param av The radiuses (radii)?
#' (>0) of the $L_1$ balls for each right vector
#' @param itermax The maximum number of iterations
#' @param eps Precision
#' @param init How to initialize the algorithm
#' @return Pseudo-singular vectors and values
#' @examples
#' X <- matrix(rnorm(20), 5, 4)
#' sparseSVD(X)
#' @author Vincent Guillemot
#' @export
sparseSVD <- function(X, R = 2L,
                 init, initLeft = NULL, initRight = NULL, seed = NULL,
                 rdsLeft = rep(1, R), rdsRight = rep(1, R),
                 grpLeft = NULL, grpRight = NULL,
                 orthogonality = "loadings",
                 OrthSpaceLeft = NULL, OrthSpaceRight = NULL,
                 projPriority = "orth",
                 projPriorityLeft = projPriority,
                 projPriorityRight = projPriority,
                 itermaxALS = 1000, itermaxPOCS = 1000,
                 epsALS = 1e-10, epsPOCS = 1e-10) {

  # Test that the arguments are valid
  garb <- runTestsSVD(X, R, init, initLeft, initRight, seed,
                   rdsLeft, rdsRight,
                   grpLeft, grpRight,
                   orthogonality, OrthSpaceLeft, OrthSpaceRight,
                   projPriority,
                   projPriorityLeft,
                   projPriorityRight)

  I <- nrow(X)
  J <- ncol(X)

  # Build initialization matrices either with SVD (prefered method)
  # or randomly
  res.init <- initializeSVD(X, I, J, R, init, initLeft, initRight, seed)
  U0 <- res.init$U0
  V0 <- res.init$V0
  # Build projections based on the arguments
  projLeft <- makeComposedProjection(projPriority = projPriorityLeft, grp = grpLeft)
  projRight <- makeComposedProjection(projPriority = projPriorityRight, grp = grpRight)

  if (is.null(OrthSpaceLeft)) OrthSpaceLeft <- matrix(0, I, 1)
  if (is.null(OrthSpaceRight)) OrthSpaceRight <- matrix(0, J, 1)
  U <- matrix(0, I, R)
  V <- matrix(0, J, R)

  iter <- matrix(NA, R, 2, dimnames = list(paste0("Dim. ", 1:R), c("Total", "ALS")))
  d <- rep(NA, R)

  for (r in 1:R) {
    ## Power Iteration with orth projection
    res.als <- als(
      X = X,                 # original matrix
      initLeft = U0[,r], initRight = V0[,r], # initialization vectors
      projLeft = projLeft, projRight = projRight,
      rdsLeft = rdsLeft[r], rdsRight = rdsRight[r],
      grpLeft = grpLeft, grpRight = grpRight,
      OrthSpaceLeft = OrthSpaceLeft,
      OrthSpaceRight = OrthSpaceRight,
      itermaxALS = 1000, itermaxPOCS = 1000,
      epsALS = 1e-10, epsPOCS = 1e-10)

    U[, r] <- res.als$u
    V[, r] <- res.als$v

    if (orthogonality == "loadings") {
      OrthSpaceLeft <- U
      OrthSpaceRight <- V
    }

    iter[r,] <- c(res.als$iterTOTAL, res.als$iterALS)
    d[r] <- res.als$d
  }

  oD <- order(d, decreasing = TRUE)
  # oD <- 1:R
  res <- list(d = d[oD], U=U[,oD], V=V[,oD], iter=iter)
  return(res)
}



makeComposedProjection <- function(projPriority, grp) {
  if (projPriority == "orth") {
    if (is.null(grp)) return(projL1L2_then_projOrth)
    return(projLGL2_then_projOrth)
  } else {
    if (is.null(grp)) return(projOrth_then_projL1L2)
    return(projOrth_then_projLGL2)
  }
}

runTestsSVD <- function(X, R, init, initLeft, initRight,
                         rdsLeft, rdsRight,
                         grpLeft, grpRight,
                         projPriority,
                         projPriorityLeft,
                         projPriorityRight,
                         itermaxALS, itermaxPOCS,
                         epsALS, epsPOCS) {

  ##### Test X ####
  if (nrow(X)==1 & ncol(X)==1)
    stop("You are attempting a gsGSVD of a scalar.")

  if (any(is.na(X)))
    stop("X should not contain missing values")

  ##### Test R ####
  if (!is.integer(R)) stop("R should be an integer.")
  if (R <= 1) stop("R should be > 1.")

  ##### Test initialization ####
  if (is.null(init)) {
    if (is.null(initLeft) | ! is.matrix(initLeft))
      stop("initLeft should be a matrix.")
    if (is.null(initRight)  | ! is.matrix(initRight))
      stop("initRight should be a matrix.")
  }
  if (! init %in% c("svd", "rand"))
    stop("init should be either svd or rand.")

  return(NULL)
}

initializeSVD <- function(X, I, J, R, init, initLeft, initRight, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  if (any(c(init, initLeft, initRight) == "svd")) {
    svdx <- svd(X, nu=R, nv=R)
  }

  if (is.null(init)) {
    if (initLeft == "svd") {
      U0 <- svdx$u
    } else if (initLeft == "rand") {
      U0 <- 1/(I-1) * mvrnorm(n = I, mu = rep(0,R),
                              Sigma = diag(R), empirical = TRUE)
    } else {
      U0 <- initLeft
    }

    if (initRight == "svd") {
      V0 <- svdx$u
    } else if (initRight == "rand") {
      V0 <- 1/(I-1) * mvrnorm(n = I, mu = rep(0,R),
                              Sigma = diag(R), empirical = TRUE)
    } else {
      V0 <- initRight
    }
  } else if (init == "svd") {
    U0 <- svdx$u
    V0 <- svdx$v
  } else if ( init=="rand") {
    U0 <- 1/(I-1) * mvrnorm(n = I, mu = rep(0,R),
                            Sigma = diag(R), empirical = TRUE)
    V0 <- 1/(J-1) * mvrnorm(n = J, mu = rep(0,R),
                            Sigma = diag(R), empirical = TRUE)
  } else {
    stop("Unkown error, contact support!")
  }

  return(list(U0 = U0, V0 = V0))
}
