#' Exclude duplicated columns with a tolerance
#'
#' @param x a marix with columns to screen
#' @param n.round to allow tolerance of the difference between values, specify the rounding digit
#' @return The unique columns of \eqn{x}.
#' @examples
#' x <- matrix(c(1,2,3, 1.001, 2.002, 3.003), nrow = 3, ncol = 2, byrow = FALSE)
#' unique.column(x, round = 2
#' @export
unique.column <- function(x, n.round) {
  suppressMessages(require("data.table"))
  xt.round <- t(round(x, n.round))
  return(as.matrix(x[,!duplicated(xt.round)]))
}
