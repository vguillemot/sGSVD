library(ExPosition)
library(GSVD)
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
  expect_equal(round(sgsvd.words.res$p,10), round(ca.authors.res$ExPosition.Data$M * ca.authors.res$ExPosition.Data$pdq$p,10))
  expect_equal(round(sgsvd.words.res$p,10), round(gsvd.words.res$p,10))
  expect_equal(round(sgsvd.words.res$q,10), round(ca.authors.res$ExPosition.Data$pdq$q,10))
  expect_equal(round(sgsvd.words.res$d,10), round(ca.authors.res$ExPosition.Data$pdq$Dv,10))
  expect_equal(round(sgsvd.words.res$fi,10), round(ca.authors.res$ExPosition.Data$fi,10))
  expect_equal(round(sgsvd.words.res$fj,10), round(ca.authors.res$ExPosition.Data$fj,10))
})

