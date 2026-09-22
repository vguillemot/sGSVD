#' L1-, L2-, L1L2- or LG- projection operators applied to a vector of numerics
#'
#' @param vec vector of numeric values
#' @param rds radius
#' @param grp vector describing the groups
#' @param OrthSpace matrix defining the orthogonal space
#' @param method which algorithm to use, either "regular" (default) or "fast"
#' @param itermax the maximum number of POCS iterations
#' @param eps precision for the POCS iterations
#'
#' @return the L1-, L2-, L1L2- or LG- projection of a vector of numerics
#' @name proj
#'
#' @examples
#' x <- c(-0.1, 1, 0.5)
#' M <- rbind(c(-1, 1, 0) / sqrt(2))
#' r <- 1
#' g <- c(1, 1, 2)
#' projL1(vec = x, rds = r)$x # = (0, 0.75, 0.25)
#' projL2(vec = x)$x # ~= (-0.1, 0.9, 0.4)
#' projL1L2(vec = x, rds = r)$x # = (0, 1, 0)
#' projLG(vec = x, rds = r, grp = g)$x # ~= (-0.1, 0.7, 0.2)
#' projLGL2(vec = x, rds = r, grp = g)$x # ~= (-0.1, 1.0, 0)
#' projOrth(vec = x, OrthSpace = M)$x # = (-0.55, 0.55, 0.50)
NULL

#' @rdname proj
#' @export

projL1 <- function(vec, rds) {
  if (normL1(vec) <= rds) {
    return(list(x = vec, lambda = 0, k = NaN))
  }
  u <- sort.int(abs(vec), decreasing = TRUE, method = "quick")
  n <- length(vec)
  ukmaok <- (cumsum(u) - rds)/(1:n)
  K <- max(which(ukmaok < u))
  lambda <- ukmaok[K]
  res <- proxL1(vec, lambda)
  return(list(x = res, lambda = lambda, k = NaN))
}

#' @rdname proj
#' @export

projL2 <- function(vec, rds = 1) {
  norm2x <- normL2(vec)
  if (norm2x < .Machine$double.eps)
    return(list(x = 0*vec, lambda = NA, k = NaN))
  res <- vec / norm2x
  return(list(x = res, lambda = NA, k = NaN))
}

#' @rdname proj
#' @export

projL1L2 <- function(vec, rds, method = "regular") {
  if (method == "fast") {
    res <- projL1L2fast(vec, rds)
  } else {
    res <- projL1L2regular(vec, rds)
  }
  return(res)
}

#' Fast (randomized pivot search) L1L2 projection
#'
#' @param vec vector of numeric values
#' @param rds radius
#'
#' @return the L1L2 projection of vec
#' @noRd
projL1L2fast <- function(vec, rds) {
  norm2_x <- normL2(vec)
  if (norm2_x < .Machine$double.eps) {
    warning("Projecting a vector with a very small L2-norm!")
    return(list(x = 0 * vec, lambda = NA, k = NaN))
  }
  if (sum(abs(vec / norm2_x)) <= rds)
    return(list(
      x = vec / norm2_x,
      lambda = max(abs(vec)),
      k = NaN
    ))
  uneq <- vec != 0
  L <- sum(!uneq)
  p <- abs(vec[uneq])
  # Check if multiple maximum
  MAX <- max(p)
  bMAX <- p == MAX
  nMAX <- sum(bMAX)
  if (rds < sqrt(nMAX)) {
    warning("Minimum radius is: ", sqrt(nMAX))
    x_soft        <- rep(0, length(vec))
    x_soft[abs(vec) == MAX]  <- 1 / sqrt(nMAX)
    return(list(x = x_soft, lambda = MAX, k = NaN))
  } else if (rds == sqrt(nMAX)) {
    warning("radius is equal to sqrt(nMAX)")
    x_soft        <- rep(0, length(vec))
    x_soft[abs(vec) == MAX]  <- 1 / sqrt(nMAX)
    # x_soft[bMAX]  <- 1 / sqrt(nMAX)
    return(list(x = x_soft, lambda = MAX, k = NaN))
  }
  #Initialize parameters
  s_1 <- s_2 <- nb <- 0
  while (T) {
    N       <- length(p)
    if (N == 0) {
      warning("length(p) = 0")
      break
    }
    # Choose next a_k
    a_k     <- p[round(runif(1, 1, N), 0)]
    while (a_k == MAX) {
      a_k     <- p[round(runif(1, 1, N), 0)]
    }
    # print(a_k)
    # Make a partition of list p
    p_inf_ak <- p < a_k
    p_sup_ak <- p > a_k
    p_high  <- p[p_inf_ak]
    p_low   <- p[p_sup_ak]
    # Evaluation decreasing rank of a_k
    nb_a_k  <- sum(p == a_k)
    k       <- nb + sum(p_sup_ak) + nb_a_k
    # Compute constraint value
    aksq <- a_k ^ 2
    s_low_1 <- sum(p_low) + nb_a_k * a_k
    s_low_2 <- ssq(p_low) + nb_a_k * aksq
    if (s_2 + s_low_2 - 2 * a_k * (s_1 + s_low_1) + k * aksq <= 1e-20) {print("help, i'm stuck !") ; next}
    psi_a_k <- (s_1 + s_low_1 - k * a_k) /
      sqrt(s_2 + s_low_2 - 2 * a_k * (s_1 + s_low_1) + k * aksq)
    ## Minor tweak: put an abs
    # psi_a_k <- (s_1 + s_low_1 - k * a_k) /
    #   sqrt(abs(s_2 + s_low_2 - 2 * a_k * (s_1 + s_low_1) + k * aksq))
    # print(s_2 + s_low_2 - 2 * a_k * (s_1 + s_low_1) + k * aksq)
    #Choose partition depending on the constraint
    # print(psi_a_k)
    # print(rds)
    if (psi_a_k > rds) {
      if (length(p_low) == 0)
        break
      p         <- p_low
    } else {
      if (length(p_high) == 0) {
        break
      } else {
        a_k_1     <- max(p_high)
        psi_a_k_1 <- (s_1 + s_low_1 - k * a_k_1) /
          sqrt(s_2 + s_low_2 - 2 * a_k_1 * (s_1 + s_low_1) + k * a_k_1 ^ 2)
        if (psi_a_k_1 > rds) {
          break
        }
        p   <- p_high
        nb  <- k
        s_1 <- s_1 + s_low_1
        s_2 <- s_2 + s_low_2
      }
    }
  }
  # Compute lambda and the soft tresholded vector to return
  lambda <- a_k -
    (rds * sqrt((k - psi_a_k ^ 2) / (k - rds ^ 2)) - psi_a_k) *
    (s_1 + s_low_1 - k * a_k) / (psi_a_k * (k))
  x_soft <- sign(vec) * pmax(0, abs(vec) - lambda)
  return(list(
    x = x_soft / normL2(x_soft) ,
    lambda = lambda,
    k = NaN
  ))
}

#' L1/L2-norm ratio (and L2-norm) of a sorted vector after soft-thresholding at each of its own entries
#'
#' Equivalent to \code{psi(xtilde, xtilde)} but computed in O(n) via cumulative
#' sums instead of O(n^2) (one \code{proxL1} pass per entry), exploiting the
#' fact that \code{xtilde} is sorted in decreasing order. Also returns the
#' L2-norm at each threshold, since \code{projL1L2regular} needs it and it is
#' already computed as an intermediate here.
#'
#' @param xtilde a vector of numeric values sorted in decreasing order
#'
#' @return a list with \code{psi} (the L1/L2-norm ratio of \code{xtilde}
#' soft-thresholded at each \code{xtilde[j]}) and \code{L2} (the corresponding
#' L2-norms)
#' @noRd
psi_sorted_self <- function(xtilde) {
  n <- length(xtilde)
  csum <- cumsum(xtilde)
  csumsq <- cumsum(xtilde^2)
  csum_prev <- c(0, csum[-n])
  csumsq_prev <- c(0, csumsq[-n])
  jm1 <- seq_len(n) - 1
  L1 <- csum_prev - jm1 * xtilde
  L2sq <- csumsq_prev - 2 * xtilde * csum_prev + jm1 * xtilde^2
  list(psi = L1 / sqrt(L2sq), L2 = sqrt(L2sq))
}

#' L1L2 projection via a sorted search on \code{psi}
#'
#' @param vec vector of numeric values
#' @param rds radius
#'
#' @return the L1L2 projection of vec
#' @noRd
projL1L2regular <- function(vec, rds) {
  norm2_x <- normL2(vec)
  if (norm2_x < .Machine$double.eps) {
    warning("Projecting a vector with a very small L2-norm!")
    return(list(x = 0 * vec, lambda = NA, k = NaN))
  }
  if (sum(abs(vec / norm2_x)) <= rds)
    return(list(
      x = vec / norm2_x,
      lambda = max(abs(vec)),
      k = NaN
    ))
  uneq <- vec != 0
  L <- sum(!uneq)
  p <- abs(vec[uneq])
  # Check if multiple maximum
  MAX <- max(p)
  bMAX <- p == MAX
  nMAX <- sum(bMAX)
  if (rds < sqrt(nMAX)) {
    warning("Minimum radius is: ", sqrt(nMAX))
    x_soft        <- rep(0, length(vec))
    x_soft[abs(vec) == MAX]  <- 1 / sqrt(nMAX)
    return(list(x = x_soft, lambda = MAX, k = NaN))
  } else if (rds == sqrt(nMAX)) {
    warning("radius is equal to sqrt(nMAX)")
    x_soft        <- rep(0, length(vec))
    x_soft[abs(vec) == MAX]  <- 1 / sqrt(nMAX)
    # x_soft[bMAX]  <- 1 / sqrt(nMAX)
    return(list(x = x_soft, lambda = MAX, k = NaN))
  }
  # 1. Take the absolute value of $\x$
  #   and sort its elements in decreasing order
  # to get $\widetilde{\x}$\;
  xtilde <- sort.int(abs(vec), decreasing = TRUE, method = "quick")
  psi_res <- psi_sorted_self(xtilde)
  psi_xtilde <- psi_res$psi
  # 2. Find i such that
  # $\psi(\widetilde{x}_{i+1})\leq c<\psi(\widetilde x_{i})$\;
  i <- max(which(psi_xtilde <= rds))
  # 3. Let $\displaystyle \delta = \frac{\normTwo{S(\widetilde{\x}	,
  #                                                 \widetilde x_i)}}{i}\left( c\sqrt{\frac{i-\psi(\widetilde x_i)^2}{i-c^2}}
  #                                                                            - \psi(\widetilde x_i)\right)$\;
  t1 <- psi_res$L2[i] / i
  t2 <- (i - psi_xtilde[i]^2) / (i - rds^2)
  t3 <- psi_xtilde[i]
  delta <- t1 * (rds * sqrt(t2) - t3)
  # 4. Compute $S(\x, \lambda)$ with $\lambda = \widetilde x_i - \delta$
  lambda <- xtilde[i] - delta
  x_soft <- sign(vec) * pmax(0, abs(vec) - lambda)
  return(list(
    x = x_soft / normL2(x_soft) ,
    lambda = lambda,
    k = NaN
  ))
}

#' L1-norm of \code{x} after soft-thresholding at \code{lambda}
#'
#' @param x vector of numeric values
#' @param lambda one or more soft-thresholding parameters
#'
#' @return the L1-norm of the soft-thresholded \code{x}, for each \code{lambda}
#' @noRd
phi <- function(x, lambda) {
  sapply(lambda, function(l) normL1(proxL1(vec = x, lambda = l)))
}

#' L1/L2-norm ratio of \code{x} after soft-thresholding at \code{lambda}
#'
#' @param x vector of numeric values
#' @param lambda one or more soft-thresholding parameters
#'
#' @return the L1/L2-norm ratio of the soft-thresholded \code{x}, for each \code{lambda}
#' @noRd
psi <- function(x, lambda) {
  sapply(lambda, function(l) {
    x_soft <- proxL1(vec = x, lambda = l)
    normL1(x_soft) / normL2(x_soft)
  })
}


#' @rdname proj
#' @export

projLG <- function(vec, rds, grp) {
  if (normLG(vec, grp) <= rds)
    return(list(x = vec, lambda = NA, k = NaN))
  vecnorm <- tapply(vec, grp, normL2)
  resproj <- projL1(vecnorm, rds)
  lambda <- resproj$lambda
  res <- proxLG(vec, lambda, grp)
  return(list(x = res, lambda = lambda, k = NaN))
}

#' @rdname proj
#' @export

projLGL2 <- function(vec, rds, grp, method = "regular") {
  if (method == "fast") {
    res <- projLGL2fast(vec, rds, grp)
  } else {
    res <- projLGL2regular(vec, rds, grp)
  }
  return(res)
}

#' Fast (randomized pivot search) group-L2 projection
#'
#' @param vec vector of numeric values
#' @param rds radius
#' @param grp vector describing the groups
#'
#' @return the group-L2 projection of vec
#' @noRd
projLGL2fast <- function(vec, rds, grp) {
  grpvec <- tapply(vec, grp, normL2)
  norm2_x <- normL2(grpvec)
  if (norm2_x < .Machine$double.eps) {
    warning("Projecting a vector with a very small L2-norm!")
    return(list(x = 0 * vec, lambda = NA, k = NaN))
  }
  if (normLG(projL2(vec)$x, grp) <= rds)
    return(list(
      x = projL2(vec)$x,
      lambda = max(abs(grpvec)),
      k = NaN
    ))
  uneq <- grpvec != 0
  L <- sum(!uneq)
  p <- abs(grpvec[uneq])
  # Check if multiple maximum
  MAX <- max(p)
  bMAX <- p == MAX
  nMAX <- sum(bMAX)
  if (rds < sqrt(nMAX)) {
    warning("Minimum radius is: ", sqrt(nMAX))
    x_soft        <- rep(0, length(vec))
    x_soft[abs(vec) == MAX]  <- 1 / sqrt(nMAX)
    return(list(x = x_soft, lambda = MAX, k = NaN))
  } else if (rds == sqrt(nMAX)) {
    warning("radius is equal to sqrt(nMAX)")
    projvec <- rep(0, length(vec))
    for (g in unique(grp)[bMAX]) {
      ig <- grp == g
      projvec[ig] <- projL2(vec[ig])$x
    }
    return(list(x = projL2(projvec)$x, lambda = MAX, k = NaN))
  }
  #Initialize parameters
  s_1 <- s_2 <- nb <- 0
  while (T) {
    N       <- length(p)
    if (N == 0) {
      warning("length(p) = 0")
      break
    }
    # Choose next a_k
    a_k     <- p[round(runif(1, 1, N), 0)]
    while (a_k == MAX) {
      a_k     <- p[round(runif(1, 1, N), 0)]
    }
    # print(a_k)
    # Make a partition of list p
    p_inf_ak <- p < a_k
    p_sup_ak <- p > a_k
    p_high  <- p[p_inf_ak]
    p_low   <- p[p_sup_ak]
    # Evaluation decreasing rank of a_k
    nb_a_k  <- sum(p == a_k)
    k       <- nb + sum(p_sup_ak) + nb_a_k
    # Compute constraint value
    aksq <- a_k ^ 2
    s_low_1 <- sum(p_low) + nb_a_k * a_k
    s_low_2 <- ssq(p_low) + nb_a_k * aksq
    psi_a_k <- (s_1 + s_low_1 - k * a_k) /
      sqrt(s_2 + s_low_2 - 2 * a_k * (s_1 + s_low_1) + k * aksq)
    #Choose partition depending on the constraint
    if (psi_a_k > rds) {
      if (length(p_low) == 0)
        break
      p         <- p_low
    } else {
      if (length(p_high) == 0) {
        break
      } else {
        a_k_1     <- max(p_high)
        psi_a_k_1 <- (s_1 + s_low_1 - k * a_k_1) /
          sqrt(s_2 + s_low_2 - 2 * a_k_1 * (s_1 + s_low_1) + k * a_k_1 ^ 2)
        if (psi_a_k_1 > rds) {
          break
        }
        p   <- p_high
        nb  <- k
        s_1 <- s_1 + s_low_1
        s_2 <- s_2 + s_low_2
      }
    }
  }
  # Compute lambda and the soft tresholded vector to return
  lambda <- a_k -
    (rds * sqrt((k - psi_a_k ^ 2) / (k - rds ^ 2)) - psi_a_k) *
    (s_1 + s_low_1 - k * a_k) / (psi_a_k * (k))
  projvec <- projL2(proxLG(vec, lambda, grp))$x
  return(list(
    x = projvec,
    lambda = lambda,
    k = NaN
  ))
}


#' Group-L2 projection via a sorted search (\code{projL1L2regular} on group norms)
#'
#' @param vec vector of numeric values
#' @param rds radius
#' @param grp vector describing the groups
#'
#' @return the group-L2 projection of vec
#' @noRd
projLGL2regular <- function(vec, rds, grp) {
  grpvec <- tapply(vec, grp, normL2)
  norm2_x <- normL2(grpvec)
  if (norm2_x < .Machine$double.eps) {
    warning("Projecting a vector with a very small L2-norm!")
    return(list(x = 0 * vec, lambda = NA, k = NaN))
  }
  if (normLG(projL2(vec)$x, grp) <= rds)
    return(list(
      x = projL2(vec)$x,
      lambda = max(abs(grpvec)),
      k = NaN
    ))
  uneq <- grpvec != 0
  L <- sum(!uneq)
  p <- abs(grpvec[uneq])
  # Check if multiple maximum
  MAX <- max(p)
  bMAX <- p == MAX
  nMAX <- sum(bMAX)
  if (rds < sqrt(nMAX)) {
    warning("Minimum radius is: ", sqrt(nMAX))
    x_soft        <- rep(0, length(vec))
    x_soft[abs(vec) == MAX]  <- 1 / sqrt(nMAX)
    return(list(x = x_soft, lambda = MAX, k = NaN))
  } else if (rds == sqrt(nMAX)) {
    warning("radius is equal to sqrt(nMAX)")
    projvec <- rep(0, length(vec))
    for (g in unique(grp)[bMAX]) {
      ig <- grp == g
      projvec[ig] <- projL2(vec[ig])$x
    }
    return(list(x = projL2(projvec)$x, lambda = MAX, k = NaN))
  }
  # Compute lambda with projL1L2regular
  lambda <- projL1L2regular(grpvec, rds)$lambda
  projvec <- projL2(proxLG(vec, lambda, grp))$x
  return(list(
    x = projvec,
    lambda = lambda,
    k = NaN
  ))
}


#' @rdname proj
#' @export

projOrth <- function(vec, OrthSpace) {
  # Mtx <- t(OrthSpace) %*% vec -> switch to crossprod
  Mtx <- crossprod(OrthSpace, vec)
  MMtx <- OrthSpace %*% Mtx
  res <- vec - MMtx
  # If the vector is almost null, then replace
  # it by a random vector
  if (normL2(res) < 1e-16) res <- projL2(rnorm(length(res)))$x
  return(list(
    x = res,
    lambda = NA,
    k = NaN
  ))
}

#' @rdname proj
#' @export

projPos <- function(vec) {
  res <- pmax(vec, 0)
  return(list(x = res, lambda = NA, k = NaN))
}


#' @rdname proj
#' @export

projL1L2_then_projOrth <- function(vec, rds, grp = NULL, OrthSpace, itermax, eps) {
  vecnew <- vecold <- vec
  for (iter in 1:itermax) {
    res.projL1L2 <- projL1L2(vecold, rds)
    vecnew <- projOrth(res.projL1L2$x, OrthSpace)$x
    if (normL2(vecnew - vecold) < eps) break
    vecold <- vecnew
  }
  return(list(x = vecnew, lambda = res.projL1L2$lambda, k = iter))
}

#' @rdname proj
#' @export
projLGL2_then_projOrth <- function(vec, rds, grp, OrthSpace, itermax, eps)  {
  vecnew <- vecold <- vec
  for (iter in 1:itermax) {
    res.projLGL2 <- projLGL2(vecold, rds, grp)
    vecnew <- projOrth(res.projLGL2$x, OrthSpace)$x
    if (normL2(vecnew - vecold) < eps) break
    vecold <- vecnew
  }
  return(list(x = vecnew, lambda = res.projLGL2$lambda, k = iter))
}

#' @rdname proj
#' @export
projOrth_then_projL1L2 <- function(vec, rds, grp = NULL, OrthSpace, itermax, eps)  {
  vecnew <- vecold <- vec
  for (iter in 1:itermax) {
    res.projL1L2 <- projL1L2(projOrth(vecold, OrthSpace)$x, rds)
    vecnew <- res.projL1L2$x
    if (normL2(vecnew - vecold) < eps) break
    vecold <- vecnew
  }
  return(list(x = vecnew, lambda = res.projL1L2$lambda, k = iter))
}

#' @rdname proj
#' @export
projOrth_then_projLGL2 <- function(vec, rds, grp, OrthSpace, itermax, eps)  {
  vecnew <- vecold <- vec
  for (iter in 1:itermax) {
    res.projLGL2 <- projLGL2(projOrth(vecold, OrthSpace)$x, rds, grp)
    vecnew <- res.projLGL2$x
    if (normL2(vecnew - vecold) < eps) break
    vecold <- vecnew
  }
  return(list(x = vecnew, lambda = res.projLGL2$lambda, k = iter))
}


#' @rdname proj
#' @export
projOrth_then_projPos_then_projL1L2 <- function(vec, rds, grp = NULL, OrthSpace, itermax, eps)  {
  vecnew <- vecold <- vec
  for (iter in 1:itermax) {
    res.projL1L2 <- projL1L2(projPos(projOrth(vecold, OrthSpace)$x)$x, rds)
    vecnew <- res.projL1L2$x
    if (normL2(vecnew - vecold) < eps) break
    vecold <- vecnew
  }
  return(list(x = vecnew, lambda = res.projL1L2$lambda, k = iter))
}

#' @rdname proj
#' @export
projOrth_then_projPos_then_projLGL2 <- function(vec, rds, grp, OrthSpace, itermax, eps)  {
  vecnew <- vecold <- vec
  for (iter in 1:itermax) {
    res.projLGL2 <- projLGL2(projPos(projOrth(vecold, OrthSpace)$x)$x, rds, grp)
    vecnew <- res.projLGL2$x
    if (normL2(vecnew - vecold) < eps) break
    vecold <- vecnew
  }
  return(list(x = vecnew, lambda = res.projLGL2$lambda, k = iter))
}


