#' L1-norm, L2-norm or LG-norm of a vector of numerics
#'
#' @param vec vector of numeric values
#' @param grp vector describing the groups
#'
#' @return the L1-norm, L2-norm or LG-norm of vec
#' @name norm
#'
#' @examples
#' x <- c(-0.1, 1, 0.5)
#' g <- c(1, 1, 2)
#' normL1(x) # = 1.6
#' normL2(x) # ~= 1.12
#' normLG(x, g) # ~= 1.5
NULL
#' @rdname norm
#' @export
normL1 <- function(vec) {
  return(sum(abs(vec)))
}
#' @rdname norm
#' @export
normL2 <- function(vec) {
  return(sqrt(sum(vec**2)))
}
#' @rdname norm
#' @export
normLG <- function(vec, grp) {
  return(sum(sqrt(rowsum(vec^2, grp))))
}


