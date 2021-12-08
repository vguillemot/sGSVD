#' @export
#'
#' @title Sparse generalized eigenvalue decomposition
#'
#' @description
#' \code{sparse_sparse_geigen} takes in a positive definite matrix of constraints (\code{W}), and a positive real number setting the required amount of sparsity in the resulting pseudo-eigen vectors.
#'
#' @param X a square, semi-positive symmetric data matrix to decompose
#' @param W \bold{W}eights -- the constraints applied to the matrix and thus the eigen vectors.
#' @param R
#' @param init
#' @param seed
#' @param rds
#' @param grp
#' @param orthogonality
#' @param OrthSpace
#' @param projPriority
#' @param itermaxALS
#' @param itermaxPOCS
#' @param epsALS
#' @param epsPOCS
#' @param k total number of components to return though the full variance will still be returned (see \code{d_full}). If 0, the full set of components are returned.
#'
#' @return A list with eight elements:
#' \item{d_full}{A vector containing the singular values of X above the tolerance threshold (based on eigenvalues).}
#' \item{l_full}{A vector containing the eigen values of X above the tolerance threshold (\code{tol}).}
#' \item{d}{A vector of length \code{min(length(d_full), k)} containing the retained singular values of X}
#' \item{l}{A vector of length \code{min(length(l_full), k)} containing the retained eigen values of X}
#' \item{v}{Eigenvectors. Dimensions are \code{ncol(X)} by k.}
#' \item{q}{Generalized eigenvectors. Dimensions are \code{ncol(X)} by k.}
#' \item{fj}{Component scores. Dimensions are \code{ncol(X)} by k.}
#'
#' @seealso \code{\link{tolerance_eigen}}, \code{\link{gsvd}} and \code{\link{gplssvd}}
#'
#' @examples
#'
#' ## (Metric) Multidimensional Scaling
#' data(wine, package="GSVD")
#' D <- as.matrix(dist(wine$objective))
#' masses <- rep(1/nrow(D), nrow(D))
#' Xi <- diag(nrow(D)) - ( rep(1,nrow(D)) %o% masses )
#' S <- Xi %*% (-(D^2) / 2) %*% t(Xi)
#' mds.res_sparse_geigen <- sparse_geigen(S)
#'
#' ## Principal components analysis: "covariance"
#' cov_X <- as.matrix(cov(wine$objective))
#' cov_pca.res_sparse_geigen <- sparse_geigen(cov_X)
#'
#' ## Principal components analysis: "correlation"
#' cor_X <- as.matrix(cor(wine$objective))
#' cor_pca.res_sparse_geigen <- sparse_geigen(cor_X)
#'
#' @author Derek Beaton
#' @keywords multivariate

sparseGEIGEN <- function(X, W, k = 0, R = 2L,
                          init = NULL, seed = NULL,
                          rds = rep(1, R),
                          grp = NULL,
                          orthogonality = "loadings",
                          OrthSpace = NULL,
                          projPriority = "orth",
                          itermaxALS = 1000, itermaxPOCS = 1000,
                          epsALS = 1e-10, epsPOCS = 1e-10){

  # preliminaries
  X_dimensions <- dim(X)
  ## stolen from MASS::ginv()
  if (length(X_dimensions) > 2 || !(is.numeric(X) || is.complex(X))){
    stop("sparse_geigen: 'X' must be a numeric or complex matrix")
  }
  if (!is.matrix(X)){
    X <- as.matrix(X)
  }

  # check square-ness here.
  if(X_dimensions[1] != X_dimensions[2]){
    stop("sparse_geigen: X must be square (i.e., have the same number of rows and columns)")
  }

  # a few things about W for stopping conditions
  W_is_missing <- missing(W)
  if(!W_is_missing){

    W_is_vector <- is.vector(W)

    if(!W_is_vector){

      if( nrow(W) != ncol(W) | nrow(W) != X_dimensions[2] ){
        stop("sparse_geigen: nrow(W) does not equal ncol(W) or ncol(X)")
      }

      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(W)){
        stop("sparse_geigen: W is empty (i.e., all 0s")
      }
    }

    if(W_is_vector){
      if(length(W)!=X_dimensions[1]){
        stop("sparse_geigen: length(W) does not equal nrow(X)")
      }

      # if you gave me all zeros, I'm stopping.
      #if(all(abs(W)<=tol)){
      if(!are_all_values_positive(W)){
        stop("sparse_geigen: W is not strictly positive values.")
      }
    }
  }

  ## convenience checks *could* be removed* if problematic
  # convenience checks & conversions; these are meant to minimize W's memory footprint
  if(!W_is_missing){
    if( !W_is_vector) {
      # if( is_identity_matrix(W) ){
      #   W_is_missing <- T
      #   W <- substitute() # neat! this makes it go missing
      # }

      if( !W_is_vector & is_diagonal_matrix(W) ){
        W <- diag(W)
        W_is_vector <- T  # now it's a vector
      }
    }

    if( W_is_vector & all(W==1) ){
      W_is_missing <- T
      W <- substitute() # neat! this makes it go missing
    }
  }

  # this manipulates X as needed
  if(!W_is_missing){
    if( W_is_vector ){

      sqrt_W <- sqrt(W)
      X <- t(t(X * sqrt_W) * sqrt_W)

    }else{

      W <- as.matrix(W)
      # sqrt_W <- W %^% (1/2)
      # sqrt_W <- sqrt_psd_matrix(W)

      ## woopsies before. my assumption previously was that W is symmetric, but I doesn't have to be. it probably should be, but that's the user's problem.
      if(isSymmetric(W)){
          ## ok so now I also check. if it's symmetric we can just go straight for it
        sqrt_W <- sqrt_psd_matrix(W)
        X <- sqrt_W %*% X %*% sqrt_W
      }else{
        X <- sqrt_psd_matrix(W) %*% X %*% sqrt_psd_matrix(t(W))
      }


    }
  }
  # sparsity parameter

  if (!is.vector(rds)) {
    warning("'rds' is not a vector, converting it automatically")
    rds <- as.vector(rds)
  }

  if (!is.numeric(rds)) {
    stop("rds must be a numeric vector")
  }

  if (length(rds) != k) {
    warning("rds is not length k, converting")
    if (length(x) == 1) rds <- rep(rds, k)
    else stop("rds could not be converted automatically to the correct format")
  }

  # all the decomposition things
  if (k <= 0) {
    k <- min(X_dimensions)
  }

  # if (missing(symmetric)) {
  #   symmetric <- isSymmetric(X)
  # }

  res <- sparseEIGEN(X, R = R,
                     init = init, seed = seed,
                     rds = rds,
                     grp = grp,
                     orthogonality = orthogonality,
                     OrthSpace = OrthSpace,
                     projPriority = projPriority,
                     itermaxALS = itermaxALS, itermaxPOCS = itermaxPOCS,
                     epsALS = epsALS, epsPOCS = epsPOCS)

  res$l_full <- res$values
    res$values <- NULL
  res$d_full <- sqrt(res$l_full)
  # res$tau <- (res$l_full/sum(res$l_full)) * 100

  components_to_return <- min(length(res$d_full),k)

  res$d <- res$d_full[1:components_to_return]
  res$l <- res$l_full[1:components_to_return]
  res$v <- res$vectors[,1:components_to_return, drop = FALSE]
    res$vectors <- NULL




  # make scores according to weights
  if(!W_is_missing){
    if(W_is_vector){

      res$q <- res$v / sqrt_W
      res$fj <- t(t(res$q * W) * res$d)

    }else{

      # res$q <- (W %^% (-1/2)) %*% res$v
      res$q <- invsqrt_psd_matrix(W) %*% res$v
      res$fj <- t(t(W %*% res$q) * res$d)

    }
  }else{

    res$q <- res$v
    res$fj <- t(t(res$q) * res$d)

  }

  rownames(res$fj) <- rownames(res$v) <- rownames(res$q) <- colnames(X)

  class(res) <- c("sparse_geigen", "GSVD", "list")
  return(res)

}
