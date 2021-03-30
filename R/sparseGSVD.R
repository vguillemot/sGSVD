#' @title Sparse Generalized Singular Value Decomposition
#' @description Constrained SVD of a matrix (wrapper of c++ functions).
#' @param X a (data) matrix;
#' @param Y a second (data) matrix; this is optional and is only used for a two-table method, Default: NULL
#' @param LW PARAM_DESCRIPTION
#' @param RW PARAM_DESCRIPTION
#' @param k the desired rank of the singular decomposition, Default: 0
#' @param tol PARAM_DESCRIPTION, Default: .Machine$double.eps
#' @param init How to initialize the algorithm, Default: 'svd'
#' @param initLeft PARAM_DESCRIPTION, Default: NULL
#' @param initRight PARAM_DESCRIPTION, Default: NULL
#' @param seed PARAM_DESCRIPTION, Default: NULL
#' @param rdsLeft The radius (>0) of the
#' $L_1$ ball for each left vector, Default: rep(1, k)
#' @param rdsRight The radius (>0) of the $L_1$ balls for each right vector, Default: rep(1, k)
#' @param grpLeft PARAM_DESCRIPTION, Default: NULL
#' @param grpRight PARAM_DESCRIPTION, Default: NULL
#' @param orthogonality PARAM_DESCRIPTION, Default: 'loadings'
#' @param OrthSpaceLeft PARAM_DESCRIPTION, Default: NULL
#' @param OrthSpaceRight PARAM_DESCRIPTION, Default: NULL
#' @param projPriority PARAM_DESCRIPTION, Default: 'orth'
#' @param projPriorityLeft PARAM_DESCRIPTION, Default: projPriority
#' @param projPriorityRight PARAM_DESCRIPTION, Default: projPriority
#' @param itermaxALS The maximum number of ALS iterations, Default: 1000
#' @param itermaxPOCS The maximum number of the POCs iterations, Default: 1000
#' @param epsALS Precision in ALS, Default: 1e-10
#' @param epsPOCS Precision in POCs, Default: 1e-10
#' @return Pseudo-singular vectors and values
#' @details DETAILS
#' @examples
#' X <- matrix(rnorm(20), 5, 4)
#' sparseSVD(X)
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname sparseGSVD
#' @author Vincent Guillemot, Ju-Chi Yu
#' @export


sparseGSVD <- function(X, Y = NULL, LW, RW, LM, RM, k = 0, tol = .Machine$double.eps,
                       init = "svd", initLeft = NULL, initRight = NULL, seed = NULL,
                       rdsLeft = rep(1, k), rdsRight = rep(1, k),
                       grpLeft = NULL, grpRight = NULL,
                       orthogonality = "loadings",
                       OrthSpaceLeft = NULL, OrthSpaceRight = NULL,
                       projPriority = "orth",
                       projPriorityLeft = projPriority,
                       projPriorityRight = projPriority,
                       itermaxALS = 1000, itermaxPOCS = 1000,
                       epsALS = 1e-10, epsPOCS = 1e-10){

  # preliminaries
  Y_is_missing <- missing(Y)
  if( !Y_is_missing ){
    if ( !is.matrix(Y) ){
      Y <- as.matrix(Y)
    }
  }
  if ( !Y_is_missing ) Y_dimensions <- dim(Y)

  X_dimensions <- dim(X)

  ## stolen from MASS::ginv()
  if (length(X_dimensions) > 2 || !(is.numeric(X) || is.complex(X))){
    stop("gsvd: 'X' must be a numeric or complex matrix")
  }
  if ( !is.matrix(X) ){
    X <- as.matrix(X)
  }


  # a few things about LW for stopping conditions
  LW_is_missing <- missing(LW)

  if ( !LW_is_missing ){

    LW_is_vector <- is.vector(LW)

    if ( !LW_is_vector ){
      if ( Y_is_missing ){
        if ( nrow(LW) != ncol(LW) | nrow(LW) != X_dimensions[1] ){
          stop("gsvd: nrow(LW) does not equal ncol(LW) or nrow(X)")
        }
      }else{
        if ( nrow(LW) != ncol(LW) | nrow(LW) != X_dimensions[2] ){
          stop("gsvd: nrow(LW) does not equal ncol(LW) or ncol(X)")
        }
      }
      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(LW)){
        stop("gsvd: LW is empty (i.e., all 0s")
      }
    }

    if(LW_is_vector){
      if ( Y_is_missing ){
        if(length(LW)!=X_dimensions[1]){
          stop("gsvd: length(LW) does not equal nrow(X)")
        }
      }else{
        if(length(LW)!=X_dimensions[2]){
          stop("gsvd: length(LW) does not equal ncol(X)")
        }
      }


      # if you gave me all zeros, I'm stopping.
      # if(all(abs(LW)<=tol)){
      if(!are_all_values_positive(LW)){
        stop("gsvd: LW is not strictly positive values")
      }
    }
  }

  # a few things about RW for stopping conditions
  RW_is_missing <- missing(RW)
  if ( !RW_is_missing ){

    RW_is_vector <- is.vector(RW)

    if ( !RW_is_vector ){
      if ( Y_is_missing ){
        if( nrow(RW) != ncol(RW) | nrow(RW) != X_dimensions[2] ){
          stop("gsvd: nrow(RW) does not equal ncol(RW) or ncol(X)")
        }
      }else{
        if( nrow(RW) != ncol(RW) | nrow(RW) != Y_dimensions[2] ){
          stop("gsvd: nrow(RW) does not equal ncol(RW) or ncol(Y)")
        }
      }
      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(RW)){
        stop("gsvd: RW is empty (i.e., all 0s")
      }
    }

    if(RW_is_vector){
      if ( Y_is_missing ){
        if(length(RW)!=X_dimensions[2]){
          stop("gsvd: length(RW) does not equal ncol(X)")
        }
      }else{
        if(length(RW)!=Y_dimensions[2]){
          stop("gsvd: length(RW) does not equal ncol(Y)")
        }
      }

      # if you gave me all zeros, I'm stopping.
      # if(all(abs(RW)<=tol)){
      if(!are_all_values_positive(RW)){
        stop("gsvd: RW is not strictly positive values")
      }
    }
  }

  # a few things about LM for stopping conditions
  LM_is_missing <- missing(LM)

  if ( !LM_is_missing ){

    LM_is_vector <- is.vector(LM)

    if ( !LM_is_vector ){
      if ( Y_is_missing ){
        stop("gsvd: RM is only used when there are two data tables")
      }else{
        if ( nrow(LM) != ncol(LM) | nrow(LM) != X_dimensions[1] ){
          stop("gsvd: nrow(LM) does not equal ncol(LM) or nrow(X)")
        }
      }
      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(LM)){
        stop("gsvd: LM is empty (i.e., all 0s")
      }
    }

    if(LM_is_vector){
      if ( Y_is_missing ){
        stop("gsvd: RM is only used when there are two data tables")
      }else{
        if(length(LM)!=X_dimensions[1]){
          stop("gsvd: length(LM) does not equal nrow(X)")
        }
      }
      # if you gave me all zeros, I'm stopping.
      if(!are_all_values_positive(LM)){
        stop("gsvd: LM is not strictly positive values")
      }
    }
  }

  # a few things about RM for stopping conditions
  RM_is_missing <- missing(RM)
  if ( !RM_is_missing ){

    RM_is_vector <- is.vector(RM)

    if ( !RM_is_vector ){
      if ( Y_is_missing ){
          stop("gsvd: RM is only used when there are two data tables")
      }else{
        if( nrow(RM) != ncol(RM) | nrow(RM) != Y_dimensions[1] ){
          stop("gsvd: nrow(RM) does not equal ncol(RM) or nrow(Y)")
        }
      }
      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(RM)){
        stop("gsvd: RM is empty (i.e., all 0s")
      }
    }

    if(RM_is_vector){
      if ( Y_is_missing ){
        stop("gsvd: RM is only used when there are two data tables")
      }else{
        if(length(RM)!=Y_dimensions[1]){
          stop("gsvd: length(RM) does not equal nrow(Y)")
        }
      }
      # if you gave me all zeros, I'm stopping.
      if(!are_all_values_positive(RM)){
        stop("gsvd: RM is not strictly positive values")
      }
    }
  }



  ## convenience checks *could* be removed* if problematic
  # convenience checks & conversions; these are meant to minimize LW's memory footprint
  if(!LW_is_missing){
    if( !LW_is_vector){

      if( is_diagonal_matrix(LW) ){
        LW <- diag(LW)
        LW_is_vector <- T  # now it's a vector
      }
    }

    if( LW_is_vector & all(LW==1) ){
      LW_is_missing <- T
      LW <- substitute() # neat! this makes it go missing
    }
  }

  if(!LM_is_missing){
    if( !LM_is_vector){

      if( is_diagonal_matrix(LM) ){
        LM <- diag(LM)
        LM_is_vector <- T  # now it's a vector
      }
    }

    if( LM_is_vector & all(LM==1) ){
      LM_is_missing <- T
      LM <- substitute() # neat! this makes it go missing
    }
  }

  # convenience checks & conversions; these are meant to minimize RW's memory footprint
  if(!RW_is_missing){
    if( !RW_is_vector ){

      if( !RW_is_vector & is_diagonal_matrix(RW) ) {
        RW <- diag(RW)
        RW_is_vector <- TRUE  # now it's a vector
      }
    }

    if( RW_is_vector & all(RW == 1) ) {
      RW_is_missing <- TRUE
      RW <- substitute() # neat! this makes it go missing
    }
  }

  if(!RM_is_missing){
    if( !RM_is_vector ){

      if( !RM_is_vector & is_diagonal_matrix(RM) ) {
        RM <- diag(RM)
        RM_is_vector <- TRUE  # now it's a vector
      }
    }

    if( RM_is_vector & all(RM == 1) ) {
      RM_is_missing <- TRUE
      RM <- substitute() # neat! this makes it go missing
    }
  }


  # this manipulates X as needed based on XLW
  if ( !LW_is_missing ) { ## plain SVD
    if ( LW_is_vector ) {
      sqrt_LW <- sqrt(LW)
      if (Y_is_missing) {   ## one-table
        X <- X * sqrt_LW
      } else {              ## two-table
        X <- t(t(X) * sqrt_LW)
      }
    } else {              ## GSVD
      LW <- as.matrix(LW)
      if (Y_is_missing) {   ## one-table
        X <- sqrt_psd_matrix(LW) %*% X
      } else {              ## two-table
      X <- X %*% sqrt_psd_matrix(LW)
    }
    }
  }

  if ( !Y_is_missing ) { ## two-tables with masses for X
    if ( !LM_is_missing ) {
      if (LM_is_vector) {
        X <- X * sqrt(LM)
      }
      else {
        LM <- as.matrix(LM)
        X <- sqrt_psd_matrix(LM) %*% X
      }
    }
  }

  # this manipulates X (or Y) as needed based on XRW
  if ( !RW_is_missing ) { ## plain SVD
    if( RW_is_vector ) {
      sqrt_RW <- sqrt(RW)
      # X <- sweep(X,2, sqrt_RW,"*") ## replace the sweep with * & t()
      if (Y_is_missing) {   ## one-table
        X <- t(t(X) * sqrt_RW)
      } else {              ## two-table
        Y <-  t(t(Y) * sqrt_RW)
      }
    } else {              ## GSVD
      RW <- as.matrix(RW)
      if (Y_is_missing) {   ## one-table
        X <- X %*% sqrt_psd_matrix(RW)
      } else {              ## two-table
        Y <- Y %*% sqrt_psd_matrix(RW)
      }
    }
  }

  if ( !Y_is_missing ) { ## two-tables with masses for Y
    if ( !RM_is_missing ) {
      if (RM_is_vector) {
        Y <- Y * sqrt(RM)
      }
      else {
        RM <- as.matrix(RM)
        Y <- sqrt_psd_matrix(RM) %*% Y
      }
    }
  }

  # all the decomposition things
  if (k <= 0) {
    k <- min(X_dimensions)
  }

  # res <- tolerance_svd(X, nu = k, nv = k, tol = tol)
  res <- sparseSVD(X = X, Y = Y, k = k,
                   init=init, initLeft = initLeft, initRight = initRight, seed = seed,
                   rdsLeft = rdsLeft, rdsRight = rdsRight,
                   grpLeft = grpLeft, grpRight = grpRight,
                   orthogonality = orthogonality,
                   OrthSpaceLeft = OrthSpaceLeft, OrthSpaceRight = OrthSpaceRight,
                   projPriority = projPriority,
                   projPriorityLeft = projPriorityLeft,
                   projPriorityRight = projPriorityRight,
                   itermaxALS = itermaxALS, itermaxPOCS = itermaxPOCS,
                   epsALS = epsALS, epsPOCS = epsPOCS)

  res$d_full <- res$d
  res$l_full <- res$d_full^2
  # res$tau <- (res$l_full/sum(res$l_full)) * 100
  components.to.return <- min(length(res$d_full), k) #a safety check
  res$d <- res$d_full[1:components.to.return]
  res$l <- res$d^2
  res$U <- res$U[,1:components.to.return, drop = FALSE]
  res$V <- res$V[,1:components.to.return, drop = FALSE]

  # make scores according to weights
  if(!LW_is_missing){
    if(LW_is_vector){

      # res$p <- sweep(res$U,1,1/sqrt_LW,"*")
      res$p <- res$U / sqrt_LW
      # res$fi <- sweep(sweep(res$p,1,LW,"*"),2,res$d,"*")
      res$fi <- t(t(res$p * LW) * res$d)

    }else{

      # res$p <- (LW %^% (-1/2)) %*% res$U
      res$p <- invsqrt_psd_matrix(LW) %*% res$U
      # res$fi <- sweep((LW %*% res$p),2,res$d,"*")
      res$fi <- t(t(LW %*% res$p) * res$d)

    }
  }else{

    res$p <- res$U
    # res$fi <- sweep(res$p,2,res$d,"*")
    res$fi <- t(t(res$p) * res$d)

  }

  if(!RW_is_missing){
    if(RW_is_vector){

      # res$q <- sweep(res$V,1,1/sqrt_RW,"*")
      res$q <- res$V / sqrt_RW
      # res$fj <- sweep(sweep(res$q,1,RW,"*"),2,res$d,"*")
      res$fj <- t(t(res$q * RW) * res$d)

    }else{

      res$q <- invsqrt_psd_matrix(RW) %*% res$V
      # res$fj <- sweep((RW %*% res$q),2,res$d,"*")
      res$fj <- t(t(RW %*% res$q) * res$d)

    }
  }else{

    res$q <- res$V
    # res$fj <- sweep(res$q,2,res$d,"*")
    res$fj <- t(t(res$q)  * res$d)

  }

  if (is.null(Y)){
    rownames(res$fi) <- rownames(res$U) <- rownames(res$p) <- rownames(X)
    rownames(res$fj) <- rownames(res$V) <- rownames(res$q) <- colnames(X)
  }else{
    rownames(res$fi) <- rownames(res$U) <- rownames(res$p) <- colnames(X)
    rownames(res$fj) <- rownames(res$V) <- rownames(res$q) <- colnames(Y)
  }

  class(res) <- c("sGSVD", "sSVD", "list")
  return(res)

}
