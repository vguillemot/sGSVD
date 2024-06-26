#' Constrained SVD of a matrix.
#'
#' @param X a (data) matrix ;
#' @param Y a second (data) matrix ; this is optional and is only used for two-table methods ;
#' @param k the desired rank of the singular decomposition ;
#' @param init how the pseudo-singular vectors are initialized ; default to "svd", meaning that they are initialized with regular singular vectors ; any other value will initialize them randomly ;
#' @param initLeft how the left pseudo-singular vectors are initialized ; default to the same value as "init" ;
#' @param initRight how the right pseudo-singular vectors are initialized ; default to the same value as "init" ;
#' @param seed a random seed for result reproducibility ; if NULL (the default), no random seed will be used ;
#' @param rdsLeft a vector of radiuses
#' (>0) of the $L_1$ or $L_G$ balls for each of the k left vectors ;
#' @param rdsRight a vector of radiuses
#' (>0) of the $L_1$ or $L_G$ balls for each of the k right vectors ;
#' @param grpLeft vector describing the groups for the left vectors ; default to one group per row ;
#' @param grpRight vector describing the groups for the right vectors ; default to one group per column ;
#' @param orthogonality wh the orthogonality constraint will
#' @param OrthSpaceLeft
#' @param OrthSpaceRight
#' @param projPriority
#' @param projPriorityLeft
#' @param projPriorityRight
#' @param itermaxALS The maximum number of ALS iterations
#' @param itermaxPOCS The maximum number of POCS iterations
#' @param epsALS Precision for ALS
#' @param epsPOCS Precision for POCS
#'
#' @return Pseudo-singular vectors and values
#' @examples
#' X <- matrix(rnorm(20), 5, 4)
#' sparseSVD(X)
#' @author Vincent Guillemot
#' @export

sparseSVD <- function(X, Y = NULL, k = 2L,
                 init = "svd", initLeft = NULL, initRight = NULL, seed = NULL,
                 rdsLeft = rep(1, k), rdsRight = rep(1, k),
                 grpLeft = NULL, grpRight = NULL,
                 orthogonality = "loadings",
                 OrthSpaceLeft = NULL, OrthSpaceRight = NULL,
                 projPriority = "orth",
                 projPriorityLeft = projPriority,
                 projPriorityRight = projPriority,
                 compute_sparsity_index = TRUE,
                 correction4SI = "gsvd",
                 itermaxALS = 1000, itermaxPOCS = 1000,
                 epsALS = 1e-10, epsPOCS = 1e-10,
                 tol.si = .Machine$double.eps) {

  if (is.null(Y)) {
    Data <- X
  }else{
    Data <- t(X) %*% Y
    if (nrow(X) != nrow(Y))
      stop ("The two data tables should have matching numbers of rows.")
    N <- nrow(X)
  }

  # Test that the arguments are valid
  garb <- runTestsSVD(Data, k, init, initLeft, initRight, seed,
                   rdsLeft, rdsRight,
                   grpLeft, grpRight,
                   orthogonality, OrthSpaceLeft, OrthSpaceRight,
                   projPriority,
                   projPriorityLeft,
                   projPriorityRight)

  I <- nrow(Data)
  J <- ncol(Data)

  # Build initialization matrices either with SVD (prefered method)
  # or randomly
  res.init <- initializeSVD(Data, I, J, k, init, initLeft, initRight, seed)
  U0 <- res.init$U0
  V0 <- res.init$V0
  # Build projections based on the arguments
  projLeft <- makeComposedProjection(projPriority = projPriorityLeft, grp = grpLeft)
  projRight <- makeComposedProjection(projPriority = projPriorityRight, grp = grpRight)

  if (is.null(OrthSpaceLeft)) OrthSpaceLeft <- matrix(0, I, 1)
  if (is.null(OrthSpaceRight)) OrthSpaceRight <- matrix(0, J, 1)
  U <- matrix(NA, I, k)
  V <- matrix(NA, J, k)

  if (!is.null(Y)) {
    U.Rv <- matrix(NA, I, k)
    V.Ru <- matrix(NA, J, k)
    Lx <- Ly <- matrix(NA, N, k)
  }

  iter <- matrix(NA, k, 2, dimnames = list(paste0("Dim. ", 1:k), c("Total", "ALS")))
  d <- rep(NA, k)

  for (r in 1:k) {
    ## Power Iteration with orth projection
    res.als <- als(
      X = Data,                 # original matrix
      initLeft = U0[,r], initRight = V0[,r], # initialization vectors
      projLeft = projLeft, projRight = projRight,
      rdsLeft = rdsLeft[r], rdsRight = rdsRight[r],
      grpLeft = grpLeft, grpRight = grpRight,
      OrthSpaceLeft = OrthSpaceLeft,
      OrthSpaceRight = OrthSpaceRight,
      itermaxALS = itermaxALS, itermaxPOCS = itermaxPOCS,
      epsALS = epsALS, epsPOCS = epsPOCS)

    if (all(res.als$u == 0) | all(res.als$v == 0)) {
      stop ("Too many components are estimated. Try extracting fewer components.")
    }

    U[, r] <- res.als$u
    V[, r] <- res.als$v


    if (!is.null(Y)) {
      U.Rv[, r] <- projL2(Data %*% V[, r])$x
      V.Ru[, r] <- projL2(t(Data) %*% U[, r])$x
    }

    if (orthogonality == "loadings") {
      OrthSpaceLeft <- U[,(1:r),drop=FALSE]
      OrthSpaceRight <- V[,(1:r),drop=FALSE]
    } else if (orthogonality == "scores") {
      if (is.null(Y))
        stop ("Y is missing! The `scores` orthogonality option is for two-table methods.")
      OrthSpaceLeft <- U.Rv[,(1:r),drop=FALSE]
      OrthSpaceRight <- V.Ru[,(1:r),drop=FALSE]
    } else if (orthogonality == "both") {
      if (is.null(Y))
        stop ("Y is missing! The `both` orthogonality option is for two-table methods.")
      ULx.bind <- cbind(U[, (1:r), drop = FALSE],
                        U.Rv[, (1:r), drop = FALSE])
      VLy.bind <- cbind(V[, (1:r), drop = FALSE],
                        V.Ru[, (1:r), drop = FALSE])
      ULx <- unique.column(ULx.bind, n.round = 10)
      VLy <- unique.column(VLy.bind, n.round = 10)

      ## Replace QR by SVD to avoid strange behaviors
      # OrthSpaceLeft <- qr.Q(qr(ULx))
      # OrthSpaceRight <- qr.Q(qr(VLy))

      res.svd.left <- svd(ULx, nv = 0)
      res.svd.right <- svd(VLy, nv = 0)

      OrthSpaceLeft <- res.svd.left$u[, res.svd.left$d > 1e-16, drop = FALSE]
      OrthSpaceRight <- res.svd.right$u[, res.svd.right$d > 1e-16, drop = FALSE]
    } else {
      stop ("Check what you entered for orthogonality. Please use eiter loadings (default), scores, or both.")
    }

    iter[r,] <- c(res.als$iterTOTAL, res.als$iterALS)
    d[r] <- res.als$d

  }
  # print(d)
  oD <- order(d, decreasing = TRUE)
  # oD <- 1:k
  res <- list(d = d[oD], u = U[, oD], v = V[, oD], order = oD, iter = iter)

  res$rdsLeft <- rdsLeft
  res$rdsRight <- rdsRight
  res$grpLeft <- grpLeft
  res$grpRight <- grpRight

  # Compute the Sparsity Index
  if (compute_sparsity_index) {
    res.SI <- sparseIndex(
        res.sgsvd = res,
        singularValues = svd(X, 0, 0)$d,
        correction = correction4SI,
        tol = tol.si)
    res$SI <- res.SI
  }

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

runTestsSVD <- function(X, k, init, initLeft, initRight,
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
  if (!is.integer(k)) stop("R should be an integer.")
  if (k <= 1) stop("K should be > 1.")

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

initializeSVD <- function(X, I, J, k, init, initLeft, initRight, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  if (any(c(init, initLeft, initRight) == "svd")) {
    svdx <- svd(X, nu=k, nv=k)
  }

  if (is.null(init)) {
    if (initLeft == "svd") {
      U0 <- svdx$u
    } else if (initLeft == "rand") {
      U0 <- 1/(I-1) * MASS::mvrnorm(n = I, mu = rep(0,k),
                              Sigma = diag(k), empirical = TRUE)
    } else {
      U0 <- initLeft
    }

    if (initRight == "svd") {
      V0 <- svdx$u
    } else if (initRight == "rand") {
      V0 <- 1/(I-1) * MASS::mvrnorm(n = I, mu = rep(0,k),
                              Sigma = diag(k), empirical = TRUE)
    } else {
      V0 <- initRight
    }
  } else if (init == "svd") {
    U0 <- svdx$u
    V0 <- svdx$v
  } else if ( init=="rand") {
    U0 <- 1/(I-1) * MASS::mvrnorm(n = I, mu = rep(0,k),
                            Sigma = diag(k), empirical = TRUE)
    V0 <- 1/(J-1) * MASS::mvrnorm(n = J, mu = rep(0,k),
                            Sigma = diag(k), empirical = TRUE)
  } else {
    stop("Unkown error, contact support!")
  }

  return(list(U0 = U0, V0 = V0))
}
