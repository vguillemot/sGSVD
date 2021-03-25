library(ExPosition)
library(GSVD)
library(GPLS)

tol  = 1e-13

# test SVD ==================================
data(words)
pca.words.res <- epPCA(words$data, graphs = FALSE)

X <- scale(as.matrix(words$data))
sgsvd.words.res <- sparseGSVD(X, init = "svd", rdsLeft = rep(sqrt(nrow(words$data)), 2), rdsRight = rep(sqrt(ncol(words$data)), 2))

test_that("sparseGSVD gives back plain SVD", {
  expect_equal(sgsvd.words.res$U, pca.words.res$ExPosition.Data$pdq$p)
  expect_equal(sgsvd.words.res$V, pca.words.res$ExPosition.Data$pdq$q)
  expect_equal(sgsvd.words.res$d, pca.words.res$ExPosition.Data$pdq$Dv)
  expect_equal(sgsvd.words.res$fi, pca.words.res$ExPosition.Data$fi)
  expect_equal(sgsvd.words.res$fj, pca.words.res$ExPosition.Data$fj)
})



# test SVD with X and Y ========================
data("wine", package = "GSVD")

X <- as.matrix(wine$objective, rownames.force = TRUE)
Y <- as.matrix(wine$subjective, rownames.force = TRUE)
X4svd <- scale(X)
Y4svd <- scale(Y)

# run PLSC
plscor_results <- pls_cor(wine$objective, wine$subjective, components = 0)
# run sPLSC with no sparsification
spls.res.both <- sparseGSVD(X4svd, Y4svd, orthogonality = "both", k = 5L, rdsLeft = rep(sqrt(ncol(wine$objective)), 5), rdsRight = rep(sqrt(ncol(wine$subjective)), 5))
spls.res.loadings <- sparseGSVD(X4svd, Y4svd, orthogonality = "loadings", k = 5L, rdsLeft = rep(sqrt(ncol(wine$objective)), 5), rdsRight = rep(sqrt(ncol(wine$subjective)), 5))
spls.res.score <- sparseGSVD(X4svd, Y4svd, orthogonality = "score", k = 5L, rdsLeft = rep(sqrt(ncol(wine$objective)), 5), rdsRight = rep(sqrt(ncol(wine$subjective)), 5))

test_that("sparsePLSC gives back plain PLSC (orthogonality = both)", {
  expect_equal(abs(spls.res.both$U), abs(plscor_results$u), tolerance = tol)
  expect_equal(abs(spls.res.both$V), abs(plscor_results$v), tolerance = tol)
  expect_equal(abs(spls.res.both$d), abs(plscor_results$d), tolerance = tol)
  expect_equal(abs(spls.res.both$fi), abs(plscor_results$fi), tolerance = tol)
  expect_equal(abs(spls.res.both$fj), abs(plscor_results$fj), tolerance = tol)
})

test_that("sparsePLSC gives back plain PLSC (orthogonality = loadings)", {
  expect_equal(abs(spls.res.loadings$U), abs(plscor_results$u), tolerance = tol)
  expect_equal(abs(spls.res.loadings$V), abs(plscor_results$v), tolerance = tol)
  expect_equal(abs(spls.res.loadings$d), abs(plscor_results$d), tolerance = tol)
  expect_equal(abs(spls.res.loadings$fi), abs(plscor_results$fi), tolerance = tol)
  expect_equal(abs(spls.res.loadings$fj), abs(plscor_results$fj), tolerance = tol)
})

test_that("sparsePLSC gives back plain PLSC (orthogonality = score)", {
  expect_equal(abs(spls.res.score$U), abs(plscor_results$u), tolerance = tol)
  expect_equal(abs(spls.res.score$V), abs(plscor_results$v), tolerance = tol)
  expect_equal(abs(spls.res.score$d), abs(plscor_results$d), tolerance = tol)
  expect_equal(abs(spls.res.score$fi), abs(plscor_results$fi), tolerance = tol)
  expect_equal(abs(spls.res.score$fj), abs(plscor_results$fj), tolerance = tol)
})


# test GSVD ================================
data(authors)
X.ca <- authors
ca.authors.res <- epCA(X.ca, graphs = FALSE)

N <- sum(X.ca)
X.ca.proc <- 1/N * X.ca
Lv <- rowSums(X.ca.proc)
Rv <- colSums(X.ca.proc)
X.ca.proc <- X.ca.proc - Lv %*% t(Rv)
LW <- 1/Lv
RW <- 1/Rv

gsvd.words.res <- gsvd(as.matrix(X.ca.proc), LW = diag(LW), RW = diag(RW))
sgsvd.words.res <- sparseGSVD(as.matrix(X.ca.proc), LW = diag(LW), RW = diag(RW), k = 2L, init = "svd", rdsLeft = rep(sqrt(nrow(X.ca)), 2), rdsRight = rep(sqrt(ncol(X.ca)), 2))

test_that("sparseGSVD gives back plain GSVD", {
  expect_equal(sgsvd.words.res$p, ca.authors.res$ExPosition.Data$M * ca.authors.res$ExPosition.Data$pdq$p, tolerance = tol)
  expect_equal(sgsvd.words.res$p, gsvd.words.res$p, tolerance = tol)
  expect_equal(sgsvd.words.res$q, ca.authors.res$ExPosition.Data$pdq$q, tolerance = tol)
  expect_equal(sgsvd.words.res$d, ca.authors.res$ExPosition.Data$pdq$Dv, tolerance = tol)
  expect_equal(sgsvd.words.res$fi, ca.authors.res$ExPosition.Data$fi, tolerance = tol)
  expect_equal(sgsvd.words.res$fj, ca.authors.res$ExPosition.Data$fj, tolerance = tol)
})

