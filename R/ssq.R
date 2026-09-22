#' Compute the sum of squares.
#'
#' @param u A vector of numerics
#' @return The sum of the squared coefficients of vector u.
#' @examples 
#' ssq(1:10)
#' @export
ssq <- function(u) sum(u**2)
