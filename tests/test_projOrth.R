rm(list = ls())

library(sGSVD)
library(microbenchmark)
library(ggplot2)

projOrth2 <- function(vec, OrthSpace) {
  Mtx <- crossprod(OrthSpace, vec)
  MMtx <- OrthSpace %*% Mtx
  res <- vec - MMtx
  return(list(
    x = res,
    lambda = NA,
    k = NaN
  ))
}
hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, `+`) }
somedata <- hilbert(100)[, 1:50]

Vs <- svd(somedata)$v[, 1:3]
x <- projL2(rnorm(50))$x
# We want to project x orthogonally to Vs

o1 <- sGSVD:::projOrth(x, Vs)$x
# o2 <- lsfit(Vs, x, intercept = FALSE)$res
o2 <- projOrth2(x, Vs)$x

crossprod(o1, Vs)
crossprod(o2, Vs)

sum((o1 - o2)**2)

mb <- microbenchmark(
  projOrth = sGSVD:::projOrth(x, Vs)$x,
  projOrth2 = projOrth2(x, Vs)$x,
  times = 1000L
)

autoplot(mb)

