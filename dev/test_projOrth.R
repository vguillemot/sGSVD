# Dev-only benchmark/regression check for the crossprod-based projOrth()
# optimization: compares sGSVD:::projOrth() against a reference
# reimplementation (projOrth2) for correctness and speed.
# Requires the package to be installed (uses :::), e.g. via `R CMD INSTALL .`.

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

set.seed(42)
Vs <- svd(somedata)$v[, 1:3]
x <- projL2(rnorm(50))$x
# We want to project x orthogonally to Vs
# (x is never near-null here, so projOrth's near-null random-vector
# fallback is intentionally not exercised by this comparison)

o1 <- sGSVD:::projOrth(x, Vs)$x
o2 <- projOrth2(x, Vs)$x

stopifnot(sum((o1 - o2)^2) < 1e-12)

mb <- microbenchmark(
  projOrth = sGSVD:::projOrth(x, Vs)$x,
  projOrth2 = projOrth2(x, Vs)$x,
  times = 1000L
)
print(mb)

if (interactive()) autoplot(mb)


