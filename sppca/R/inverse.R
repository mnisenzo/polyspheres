
#' @title Nadaraya--Watson inversion from the hypersphere to the polysphere
#'
#' @description Kernel interpolator \eqn{\hat m_h: \mathcal{S}^p \to
#' \mathcal{S}^{d_1, \ldots, d_r}} built from the pairs
#' \eqn{(\boldsymbol{y}_i, \boldsymbol{x}_i)}, \eqn{i = 1, \ldots, n}:
#' \deqn{\hat m_h(\boldsymbol{y}) := \mathrm{proj}\Big(\sum_{i = 1}^n
#' W_{i, h}(\boldsymbol{y}) \boldsymbol{x}_i\Big), \quad
#' W_{i, h}(\boldsymbol{y}) \propto K\{\rho \arccos(\boldsymbol{y}'
#' \boldsymbol{y}_i) / h\},}
#' that is, a weighted extrinsic mean with intrinsic weights, with Gaussian
#' kernel \eqn{K}.
#'
#' @param y matrix of size \code{c(m, p + 1)} with evaluation points.
#' @param Y matrix of size \code{c(n, p + 1)} with the hyperspherical
#' images of the sample.
#' @param X matrix of size \code{c(n, sum(d) + r)} with the sample.
#' @param d vector of dimensions of the polysphere.
#' @param h bandwidth.
#' @param rho radius of \eqn{\mathcal{S}^p} (distances are scaled by it).
#' @return A matrix of size \code{c(m, sum(d) + r)}.
#' @export
nw_sph_polysph <- function(y, Y, X, d, h, rho = 1) {

  Th <- rho * dist_sph_cross(y, Y)
  # Log-sum-exp stabilization of the Gaussian weights
  L <- -0.5 * (Th / h)^2
  L <- L - apply(L, 1, max)
  W <- exp(L)
  W <- W / rowSums(W)
  return(proj_polysph(W %*% X, d = d))

}


#' @title Leave-one-out cross-validation bandwidth for
#' \code{nw_sph_polysph}
#'
#' @description Minimizes
#' \deqn{\mathrm{CV}(h) = \sum_{i = 1}^n \delta^2\{\hat m_{-i, h}
#' (\boldsymbol{y}_i), \boldsymbol{x}_i\}.}
#' The leave-one-out fits are obtained in closed form by zeroing the
#' diagonal of the kernel matrix before normalizing the weights (no refits),
#' so each evaluation costs \eqn{O(n^2 (d_1 + \cdots + d_r + r))}.
#'
#' @inheritParams nw_sph_polysph
#' @param h_range range of bandwidths (in the units of
#' \eqn{\rho\arccos}) for the search. Defaults to
#' \code{c(q01, max)} of the pairwise distances.
#' @param plot_cv plot the CV curve?
#' @return A list with \code{h} (CV bandwidth) and \code{cv} (the CV curve
#' on a log-grid).
#' @export
bw_cv_nw_sph_polysph <- function(Y, X, d, rho = 1, h_range = NULL,
                                 plot_cv = FALSE) {

  Th <- rho * dist_sph_cross(Y, Y)
  if (is.null(h_range)) {

    h_range <- c(quantile(Th[upper.tri(Th)], 0.01), max(Th))

  }
  cv <- function(log_h) {
    L <- -0.5 * (Th / exp(log_h))^2
    diag(L) <- -Inf
    L <- L - apply(L, 1, max)
    W <- exp(L)
    W <- W / rowSums(W)
    fit <- proj_polysph(W %*% X, d = d)
    sum(dist_polysph_pair(fit, X, d = d, squared = TRUE))
  }
  grid <- seq(log(h_range[1]), log(h_range[2]), l = 40)
  vals <- sapply(grid, cv)
  i <- which.min(vals)
  lo <- grid[max(1, i - 1)]
  hi <- grid[min(length(grid), i + 1)]
  opt <- optimize(cv, interval = c(lo, hi))
  if (plot_cv) {

    plot(exp(grid), vals, type = "o", log = "x", xlab = "h", ylab = "CV(h)")
    abline(v = exp(opt$minimum), col = 2)

  }
  return(list(h = exp(opt$minimum), cv = cbind(h = exp(grid), cv = vals)))

}


#' @title Local polynomial mode of variation on the polysphere
#'
#' @description Estimates the conditional extrinsic mean
#' \eqn{t \mapsto \mathrm{proj}\, E[\boldsymbol{X} \mid \xi_j = t]} of the
#' polyspherical sample given its \eqn{j}-th principal score
#' \eqn{\xi_j}, by local constant (\code{p = 0}) or local linear
#' (\code{p = 1}) regression (\code{\link{ll_scores_polysph}}; for
#' \code{p = 0} it coincides with \code{\link[polykde]{kre_polysph}}).
#' The bandwidth, in units of the standard deviation of the scores, is
#' chosen by leave-one-out cross-validation.
#'
#' @inheritParams nw_sph_polysph
#' @param t evaluation points for the score.
#' @param scores vector of size \code{n} with the scores.
#' @param p degree of the local fit.
#' @param h bandwidth multiplier. If \code{NULL}, chosen by leave-one-out
#' cross-validation.
#' @return A matrix of size \code{c(length(t), sum(d) + r)} with attribute
#' \code{"h"}.
#' @export
mode_polysph <- function(t, scores, X, d, p = 1, h = NULL) {

  # Bandwidths are expressed in units of sd(scores), as in ll_scores_polysph
  if (is.null(h)) {

    h_grid <- 10^seq(-1, 0.5, l = 31)
    cv <- sapply(h_grid, function(hh) {
      fit <- ll_scores_polysph(cbind(scores), cbind(scores), X, d = d, h = hh,
                               deg = p, loo = TRUE)
      sum(dist_polysph_pair(fit, X, d = d, squared = TRUE))
    })
    h <- h_grid[which.min(cv)]

  }
  out <- ll_scores_polysph(cbind(t), cbind(scores), X, d = d, h = h, deg = p)
  attr(out, "h") <- h
  return(out)

}


#' @title Local linear regression of polyspherical data on scores
#'
#' @description Estimates \eqn{\boldsymbol{s} \mapsto \mathrm{proj}\,
#' E[\boldsymbol{X} \mid \boldsymbol{\xi}^{(q)} = \boldsymbol{s}]}, the
#' conditional extrinsic mean of the sample given its first \eqn{q}
#' principal scores, by local linear (\code{deg = 1}) or local constant
#' (\code{deg = 0}) regression with a Gaussian product kernel whose
#' bandwidths are \eqn{h} times the standard deviations of the scores,
#' followed by projection onto the polysphere. This provides a common,
#' method-agnostic inversion of the scores of any dimension-reduction
#' method; with \code{loo = TRUE} the fits at the sample scores are
#' leave-one-out, giving honest reconstruction errors.
#'
#' @param s matrix of size \code{c(m, q)} with evaluation scores.
#' @param S matrix of size \code{c(n, q)} with the sample scores.
#' @param X matrix of size \code{c(n, sum(d) + r)} with the sample.
#' @param d vector of dimensions.
#' @param h bandwidth multiplier.
#' @param deg degree of the local fit, \code{0} or \code{1}.
#' @param loo leave-one-out fits (\code{s} must then be \code{S})?
#' @param ridge ridge regularization of the local design.
#' @return A matrix of size \code{c(m, sum(d) + r)}.
#' @export
ll_scores_polysph <- function(s, S, X, d, h, deg = 1, loo = FALSE,
                              ridge = 1e-6) {

  s <- rbind(s)
  S <- cbind(S)
  q <- ncol(S)
  sds <- apply(S, 2, sd)
  hs <- h * sds
  out <- matrix(0, nrow = nrow(s), ncol = ncol(X))
  for (i in seq_len(nrow(s))) {

    Z <- sweep(S, 2, s[i, ]) %*% diag(1 / hs, q)
    w <- exp(-0.5 * rowSums(Z^2))
    if (loo) w[i] <- 0
    w <- w / max(sum(w), .Machine$double.xmin)
    if (deg == 0) {

      out[i, ] <- colSums(w * X)

    } else {

      D <- cbind(1, Z)
      A <- crossprod(D, w * D) + ridge * diag(c(0, rep(1, q)))
      B <- crossprod(D, w * X)
      out[i, ] <- tryCatch(solve(A, B)[1, ], error = function(e) colSums(w * X))

    }

  }
  return(proj_polysph(out, d = d))

}


#' @title Cross-validated reconstruction errors from scores
#'
#' @description For \eqn{q = 0, \ldots, Q}, computes the leave-one-out
#' residual sum of squared polyspherical distances
#' \eqn{\sum_i \delta^2\{\boldsymbol{x}_i, \hat m_{-i}(\boldsymbol{\xi}_i^{(q)})\}}
#' of \code{\link{ll_scores_polysph}}, with the bandwidth multiplier chosen
#' to minimize it over \code{h_grid}. Rank \eqn{q = 0} uses the extrinsic
#' mean.
#'
#' @inheritParams ll_scores_polysph
#' @param Q maximum rank.
#' @param h_grid grid of bandwidth multipliers.
#' @return A vector of size \code{Q + 1} with attribute \code{"h"}.
#' @export
residual_scores <- function(S, X, d, Q = 5, deg = 1,
                            h_grid = 10^seq(-0.7, 0.5, l = 13)) {

  mu <- mean_polysph(X, d = d)
  r0 <- sum(dist_polysph_pair(X, matrix(mu, nrow(X), length(mu), byrow = TRUE),
                              d = d, squared = TRUE))
  hs <- numeric(Q)
  res <- sapply(seq_len(Q), function(q) {
    cv <- sapply(h_grid, function(h) {
      fit <- ll_scores_polysph(S[, 1:q, drop = FALSE], S[, 1:q, drop = FALSE],
                               X, d = d, h = h, deg = deg, loo = TRUE)
      sum(dist_polysph_pair(fit, X, d = d, squared = TRUE))
    })
    hs[q] <<- h_grid[which.min(cv)]
    min(cv)
  })
  out <- c(r0, res)
  names(out) <- 0:Q
  attr(out, "h") <- hs
  return(out)

}
