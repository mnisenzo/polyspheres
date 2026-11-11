
#' @title Nadaraya--Watson inversion from the hypersphere to the polysphere
#'
#' @description Kernel interpolator \eqn{\hat m_h: \mathcal{S}^p \to
#' \mathcal{S}^{d_1, \ldots, d_r}} built from the pairs
#' \eqn{(\boldsymbol{y}_i, \boldsymbol{x}_i)}, \eqn{i = 1, \ldots, n}:
#' \deqn{\hat m_h(\boldsymbol{y}) := \mathrm{proj}\Big(\sum_{i = 1}^n
#' W_{i, h}(\boldsymbol{y}) \boldsymbol{x}_i\Big), \quad
#' W_{i, h}(\boldsymbol{y}) \propto K\{\rho \arccos(\boldsymbol{y}'
#' \boldsymbol{y}_i) / h\},}
#' that is, a weighted extrinsic mean with intrinsic weights, for the
#' Gaussian kernel \eqn{K(t) = e^{-t^2 / 2}}.
#'
#' @param y matrix of size \code{c(m, p + 1)} with evaluation points.
#' @param Y matrix of size \code{c(n, p + 1)} with the hyperspherical
#' images of the sample.
#' @param X matrix of size \code{c(n, sum(d) + r)} with the sample.
#' @param d vector of size \code{r} with the dimensions of the polysphere.
#' @param h bandwidth, in the units of \eqn{\rho \arccos}.
#' @param rho radius of \eqn{\mathcal{S}^p}. Defaults to \code{1}.
#' @return A matrix of size \code{c(m, sum(d) + r)} with points on the
#' polysphere.
#' @details The weights are computed in the log scale and shifted by their
#' rowwise maxima before exponentiation, so the estimator is well defined
#' for arbitrarily small bandwidths (it then tends to the nearest-neighbor
#' interpolator).
#' @examples
#' set.seed(1)
#' X <- r_circles_polysph(n = 50, r = 2)
#' th <- attr(X, "theta")
#' Y <- cbind(cos(th), sin(th))
#' y <- cbind(cos(c(0, pi)), sin(c(0, pi)))
#' nw_sph_polysph(y = y, Y = Y, X = X, d = c(2, 2), h = 0.2)
#' @export
nw_sph_polysph <- function(y, Y, X, d, h, rho = 1) {

  # Check inputs
  y <- rbind(y)
  if (ncol(y) != ncol(Y)) {

    stop("y and Y must have the same number of columns.")

  }
  if (nrow(Y) != nrow(X)) {

    stop("Y and X must have the same number of rows.")

  }

  # Squared scaled distances, shifted by their rowwise minima (log-sum-exp)
  Th2 <- (rho * dist_sph_cross(y, Y))^2
  Th2 <- Th2 - row_min(Th2)

  # Normalized Gaussian weights and projected weighted means
  W <- exp(-0.5 * Th2 / h^2)
  W <- W / rowSums(W)
  return(proj_polysph(W %*% X, d = d))

}


#' @title Leave-one-out cross-validation bandwidth for
#' \code{nw_sph_polysph}
#'
#' @description Minimizes
#' \deqn{\mathrm{CV}(h) = \sum_{i = 1}^n \delta^2\{\hat m_{-i, h}
#' (\boldsymbol{y}_i), \boldsymbol{x}_i\},}
#' where \eqn{\hat m_{-i, h}} is \code{\link{nw_sph_polysph}} computed
#' without the \eqn{i}-th pair.
#'
#' @inheritParams nw_sph_polysph
#' @param h_range range of bandwidths for the search. Defaults to the
#' \code{0.01} quantile and the maximum of the scaled pairwise distances.
#' @param n_grid size of the logarithmic grid used to bracket the minimum,
#' which is then refined by \code{\link[stats]{optimize}}. Defaults to
#' \code{15}.
#' @param n_cv maximum number of observations at which the leave-one-out
#' loss is evaluated. If \code{n > n_cv}, a random subsample of size
#' \code{n_cv} is used (each fit still uses all the other observations),
#' which gives an unbiased estimate of the CV criterion at cost
#' \eqn{O(n_{\mathrm{cv}} n)}. Defaults to \code{1000}.
#' @param plot_cv plot the CV curve? Defaults to \code{FALSE}.
#' @return A list with \code{h} (the CV bandwidth) and \code{cv} (a matrix
#' with the CV curve on the grid).
#' @details Since the weights are normalized, the leave-one-out weights are
#' obtained by setting the diagonal of the kernel matrix to zero before the
#' normalization, \eqn{W_{-i, j, h} = W_{j, h}(\boldsymbol{y}_i) /
#' \{1 - W_{i, h}(\boldsymbol{y}_i)\}}, so no refitting is needed. The
#' rowwise shifts that stabilize the weights do not depend on \eqn{h} and are
#' computed once. Each evaluation of \eqn{\mathrm{CV}(h)} costs one
#' \eqn{n \times n} by \eqn{n \times (d_1 + \cdots + d_r + r)} matrix product.
#' @examples
#' set.seed(1)
#' X <- r_circles_polysph(n = 60, r = 2)
#' th <- attr(X, "theta")
#' Y <- cbind(cos(th), sin(th))
#' bw_cv_nw_sph_polysph(Y = Y, X = X, d = c(2, 2))$h
#' @export
bw_cv_nw_sph_polysph <- function(Y, X, d, rho = 1, h_range = NULL,
                                 n_grid = 15, n_cv = 1000, plot_cv = FALSE) {

  # Check inputs
  if (nrow(Y) != nrow(X)) {

    stop("Y and X must have the same number of rows.")

  }

  # Rows at which the leave-one-out loss is evaluated
  n <- nrow(Y)
  ind <- if (n > n_cv) sort(sample(n, n_cv)) else seq_len(n)

  # Scaled squared distances without the diagonal, shifted rowwise once
  Th <- rho * dist_sph_cross(Y[ind, , drop = FALSE], Y)
  if (is.null(h_range)) {

    h_range <- c(quantile(Th[Th > 0], 0.01), max(Th))

  }
  Th2 <- Th^2
  Th2[cbind(seq_along(ind), ind)] <- Inf
  Th2 <- Th2 - row_min(Th2)
  X_ind <- X[ind, , drop = FALSE]

  # Leave-one-out CV loss
  cv <- function(log_h) {

    W <- exp(-0.5 * Th2 / exp(2 * log_h))
    W <- W / rowSums(W)
    fit <- proj_polysph(W %*% X, d = d)
    sum(dist_polysph_pair(fit, X_ind, d = d, squared = TRUE))

  }

  # Grid search followed by a local refinement
  grid <- seq(log(h_range[1]), log(h_range[2]), l = n_grid)
  vals <- vapply(grid, cv, numeric(1))
  i <- which.min(vals)
  opt <- optimize(cv, interval = grid[c(max(1, i - 1), min(n_grid, i + 1))])
  if (plot_cv) {

    plot(exp(grid), vals, type = "o", log = "x", xlab = "h", ylab = "CV(h)")
    abline(v = exp(opt$minimum), col = 2)

  }
  return(list(h = exp(opt$minimum), cv = cbind(h = exp(grid), cv = vals)))

}


#' @title Local linear regression of polyspherical data on scores
#'
#' @description Estimates \eqn{\boldsymbol{s} \mapsto \mathrm{proj}\,
#' E[\boldsymbol{X} \mid \boldsymbol{\xi}^{(q)} = \boldsymbol{s}]}, the
#' conditional extrinsic mean of a polyspherical sample given its first
#' \eqn{q} principal scores, by local linear (\code{deg = 1}) or local
#' constant (\code{deg = 0}) regression with a Gaussian product kernel,
#' followed by projection onto the polysphere. This is a method-agnostic
#' inversion of the scores of any dimension-reduction method. With
#' \code{loo = TRUE}, fits at the sample scores are leave-one-out, giving
#' honest reconstruction errors.
#'
#' @param s matrix of size \code{c(m, q)} with evaluation scores.
#' @param S matrix of size \code{c(n, q)} with the sample scores.
#' @param X matrix of size \code{c(n, sum(d) + r)} with the sample.
#' @param d vector of size \code{r} with the dimensions of the polysphere.
#' @param h bandwidth multiplier: the bandwidth of the \eqn{k}-th score is
#' \code{h * sd(S[, k])}.
#' @param deg degree of the local fit, \code{0} or \code{1}. Defaults to
#' \code{1}.
#' @param loo leave-one-out fits? If \code{TRUE}, \code{s} must be \code{S}.
#' Defaults to \code{FALSE}.
#' @param ridge ridge regularization of the local slopes, relative to the
#' local sample size. Defaults to \code{1e-6}.
#' @return A matrix of size \code{c(m, sum(d) + r)} with points on the
#' polysphere.
#' @details The local linear estimator is linear in the responses,
#' \eqn{\hat m(\boldsymbol{s}) = \sum_i \ell_i(\boldsymbol{s})
#' \boldsymbol{x}_i}, with equivalent-kernel weights
#' \eqn{\ell_i(\boldsymbol{s}) = w_i\, \boldsymbol{e}_1'\boldsymbol{M}^{-1}
#' (1, \boldsymbol{z}_i')'}, where \eqn{\boldsymbol{z}_i} are the scaled
#' differences, \eqn{w_i} the kernel weights, and \eqn{\boldsymbol{M}} the
#' \eqn{(q + 1) \times (q + 1)} matrix of local weighted moments. The
#' moments are computed for all evaluation points at once with matrix
#' operations, only the small systems \eqn{\boldsymbol{M}\boldsymbol{a} =
#' \boldsymbol{e}_1} are solved pointwise, and all fits follow from a single
#' product \eqn{\boldsymbol{L}\boldsymbol{X}}. The cost is
#' \eqn{O(m n q^2 + m n D)} instead of \eqn{O(m n q D)} for
#' \eqn{D = d_1 + \cdots + d_r + r}.
#' @examples
#' set.seed(1)
#' X <- r_circles_polysph(n = 80, r = 3)
#' S <- cbind(attr(X, "theta"))
#' fit <- ll_scores_polysph(s = S, S = S, X = X, d = rep(2, 3), h = 0.3,
#'                          loo = TRUE)
#' mean(dist_polysph_pair(fit, X, d = rep(2, 3)))
#' @export
ll_scores_polysph <- function(s, S, X, d, h, deg = 1, loo = FALSE,
                              ridge = 1e-6) {

  # Check inputs
  s <- rbind(s)
  S <- cbind(S)
  q <- ncol(S)
  if (ncol(s) != q) {

    stop("s and S must have the same number of columns.")

  }
  if (nrow(S) != nrow(X)) {

    stop("S and X must have the same number of rows.")

  }
  if (loo && nrow(s) != nrow(S)) {

    stop("loo = TRUE requires s = S.")

  }

  # Scaled differences Z[[k]][i, j] = (S[j, k] - s[i, k]) / h_k
  hs <- h * apply(S, 2, sd)
  Z <- lapply(seq_len(q), function(k) outer(-s[, k], S[, k], "+") / hs[k])

  # Gaussian product kernel weights (m x n). Each row is rescaled by its
  # largest weight, which leaves the (local constant or local linear)
  # estimator unchanged but avoids underflow for isolated evaluation points
  E <- Reduce(`+`, lapply(Z, `^`, 2))
  if (loo) {

    diag(E) <- Inf

  }
  W <- exp(-0.5 * (E - row_min(E)))

  # Equivalent-kernel weights
  if (deg == 0) {

    L <- W

  } else {

    # Local moments: m0 = sum w, m1_k = sum w z_k, m2_kl = sum w z_k z_l
    m0 <- rowSums(W)
    WZ <- lapply(Z, `*`, W)
    m1 <- vapply(WZ, rowSums, numeric(nrow(s)))
    m1 <- matrix(m1, nrow = nrow(s))
    a <- matrix(0, nrow = nrow(s), ncol = q + 1)
    M <- diag(q + 1)
    for (i in seq_len(nrow(s))) {

      M[1, 1] <- m0[i]
      M[1, -1] <- M[-1, 1] <- m1[i, ]
      for (k in seq_len(q)) {

        for (l in k:q) {

          M[k + 1, l + 1] <- M[l + 1, k + 1] <- sum(WZ[[k]][i, ] * Z[[l]][i, ])

        }
        M[k + 1, k + 1] <- M[k + 1, k + 1] + ridge * m0[i]

      }
      a[i, ] <- tryCatch(solve(M, c(1, rep(0, q))),
                         error = function(e) c(1 / m0[i], rep(0, q)))

    }

    # l_i(s) = w_i (a_0 + sum_k a_k z_ik)
    L <- W * a[, 1]
    for (k in seq_len(q)) {

      L <- L + WZ[[k]] * a[, k + 1]

    }

  }

  # Normalize (exact for deg = 0; the local linear weights already sum to 1)
  if (deg == 0) {

    L <- L / pmax(rowSums(L), .Machine$double.xmin)

  }
  return(proj_polysph(L %*% X, d = d))

}


#' @title Local linear mode of variation on the polysphere
#'
#' @description Estimates the conditional extrinsic mean
#' \eqn{t \mapsto \mathrm{proj}\, E[\boldsymbol{X} \mid \xi = t]} of a
#' polyspherical sample given a scalar score \eqn{\xi}, by local constant
#' (\code{p = 0}) or local linear (\code{p = 1}) regression. This is the
#' estimator \code{\link[polykde]{kre_polysph}} with a bandwidth chosen by
#' leave-one-out cross-validation, computed with
#' \code{\link{ll_scores_polysph}}.
#'
#' @inheritParams nw_sph_polysph
#' @param t vector with the evaluation points of the score.
#' @param scores vector of size \code{n} with the scores.
#' @param p degree of the local fit. Defaults to \code{1}.
#' @param h bandwidth multiplier (in units of \code{sd(scores)}). If
#' \code{NULL} (default), chosen by leave-one-out cross-validation.
#' @param h_grid grid of bandwidth multipliers for the cross-validation.
#' @return A matrix of size \code{c(length(t), sum(d) + r)} with attribute
#' \code{"h"}, the bandwidth multiplier used.
#' @examples
#' set.seed(1)
#' X <- r_circles_polysph(n = 80, r = 3, theta_range = c(0, pi))
#' m <- mode_polysph(t = c(0.5, 1.5, 2.5), scores = attr(X, "theta"), X = X,
#'                   d = rep(2, 3))
#' attr(m, "h")
#' @export
mode_polysph <- function(t, scores, X, d, p = 1, h = NULL,
                         h_grid = 10^seq(-1, 0.5, l = 31)) {

  # Cross-validated bandwidth
  if (is.null(h)) {

    cv <- vapply(h_grid, function(hh) {
      fit <- ll_scores_polysph(cbind(scores), cbind(scores), X, d = d,
                               h = hh, deg = p, loo = TRUE)
      sum(dist_polysph_pair(fit, X, d = d, squared = TRUE))
    }, numeric(1))
    h <- h_grid[which.min(cv)]

  }

  # Fit
  out <- ll_scores_polysph(cbind(t), cbind(scores), X, d = d, h = h, deg = p)
  attr(out, "h") <- h
  return(out)

}


#' @title Cross-validated reconstruction errors from scores
#'
#' @description For \eqn{q = 0, \ldots, Q}, computes the leave-one-out
#' residual sum of squared polyspherical distances
#' \deqn{\sum_{i = 1}^n \delta^2\{\boldsymbol{x}_i,
#' \hat m_{-i}(\boldsymbol{\xi}_i^{(q)})\}}
#' of \code{\link{ll_scores_polysph}}, with the bandwidth multiplier chosen
#' to minimize it over \code{h_grid}. Rank \eqn{q = 0} uses the extrinsic
#' mean. One minus the ratio to the \eqn{q = 0} value is the leave-one-out
#' explained variance \eqn{\mathrm{EV}(q)}, which compares the scores of
#' different methods with a common inversion.
#'
#' @inheritParams ll_scores_polysph
#' @param Q maximum rank. Defaults to \code{5}.
#' @param h_grid grid of bandwidth multipliers.
#' @return A vector of size \code{Q + 1} with attribute \code{"h"}, the
#' selected bandwidth multipliers.
#' @examples
#' set.seed(1)
#' X <- r_circles_polysph(n = 80, r = 3)
#' tp <- tangent_pca(X, d = rep(2, 3))
#' res <- residual_scores(tp$scores, X, d = rep(2, 3), Q = 2)
#' 1 - res / res[1]
#' @export
residual_scores <- function(S, X, d, Q = 5, deg = 1,
                            h_grid = 10^seq(-0.7, 0.5, l = 13)) {

  # Rank zero: extrinsic mean
  S <- cbind(S)
  Q <- min(Q, ncol(S))
  mu <- mean_polysph(X, d = d)
  r0 <- sum(dist_polysph_pair(X, rep_rows(mu, nrow(X)), d = d, squared = TRUE))

  # Ranks 1, ..., Q, each with its own CV bandwidth
  hs <- numeric(Q)
  res <- numeric(Q)
  for (q in seq_len(Q)) {

    Sq <- S[, seq_len(q), drop = FALSE]
    cv <- vapply(h_grid, function(h) {
      fit <- ll_scores_polysph(Sq, Sq, X, d = d, h = h, deg = deg, loo = TRUE)
      sum(dist_polysph_pair(fit, X, d = d, squared = TRUE))
    }, numeric(1))
    hs[q] <- h_grid[which.min(cv)]
    res[q] <- min(cv)

  }
  out <- c(r0, res)
  names(out) <- 0:Q
  attr(out, "h") <- hs
  return(out)

}
