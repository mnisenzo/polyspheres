
#' @title Competing dimension-reduction methods on the polysphere
#'
#' @description Benchmarks for \code{\link{sppca}}, all returning scores
#' and rank-\eqn{q} reconstructions on the polysphere:
#' \itemize{
#'   \item \code{tangent_pca}: PCA of the logarithmic maps at the
#'   (blockwise) Fréchet mean, mapped back by the exponential map
#'   (principal geodesic analysis in its linearized form).
#'   \item \code{extrinsic_pca}: PCA of the embedded coordinates in
#'   \eqn{\mathbb{R}^{d_1 + \cdots + d_r + r}}, reconstructions projected
#'   back to the polysphere.
#'   \item \code{cpns}: composite principal nested spheres (Pizer et al.,
#'   2013): PNS on each \eqn{\mathcal{S}^{d_j}}, followed by PCA of the
#'   concatenated Euclideanized PNS residuals; reconstructions are mapped
#'   back through each marginal PNS.
#' }
#'
#' @inheritParams polysph_geometry
#' @param X matrix of size \code{c(n, sum(d) + r)} with the sample.
#' @param q rank of the reconstruction.
#' @param fit output of the corresponding fitting function.
#' @param iter number of iterations for the blockwise Fréchet mean.
#' @return Fitting functions return a list with \code{scores} and the
#' ingredients for \code{reconstruct_*}; the latter return a matrix of
#' size \code{c(n, sum(d) + r)}.
#' @references
#' Pizer, S. M., Jung, S., Goswami, D., Vicory, J., Zhao, X., Chaudhuri,
#' R., Damon, J. N., Huckemann, S. and Marron, J. S. (2013). Nested sphere
#' statistics of skeletal models. In \emph{Innovations for Shape Analysis},
#' pp. 93--115. Springer, Berlin.
#' @examples
#' set.seed(1)
#' X <- r_circles_polysph(n = 100, r = 3)
#' d <- rep(2, 3)
#' tp <- tangent_pca(X, d = d)
#' cp <- cpns(X, d = d)
#' # Residual sum of squares of the rank-1 reconstructions
#' sum(dist_polysph_pair(reconstruct_tangent_pca(tp, 1, d = d), X, d = d)^2)
#' sum(dist_polysph_pair(reconstruct_cpns(cp, 1, d = d), X, d = d)^2)
#' @name competitors
NULL


#' @rdname competitors
#' @export
frechet_mean_polysph <- function(X, d, iter = 50) {

  mu <- mean_polysph(X, d = d)
  for (it in seq_len(iter)) {

    v <- colMeans(log_polysph(X, mu = mu, d = d))
    mu_new <- drop(exp_polysph(rbind(v), mu = mu, d = d))
    if (max(abs(mu_new - mu)) < 1e-10) {

      break

    }
    mu <- mu_new

  }
  return(mu)

}


#' @rdname competitors
#' @export
tangent_pca <- function(X, d) {

  mu <- frechet_mean_polysph(X, d = d)
  V <- log_polysph(X, mu = mu, d = d)
  vbar <- colMeans(V)
  s <- svd(sweep(V, 2, vbar))
  scores <- s$u %*% diag(s$d)
  return(list(mu = mu, vbar = vbar, rot = s$v, scores = scores,
              var = 100 * s$d^2 / sum(s$d^2)))

}


#' @rdname competitors
#' @export
reconstruct_tangent_pca <- function(fit, q, d) {

  V <- matrix(fit$vbar, nrow = nrow(fit$scores), ncol = length(fit$vbar),
              byrow = TRUE)
  if (q > 0) {

    V <- V + fit$scores[, 1:q, drop = FALSE] %*% t(fit$rot[, 1:q, drop = FALSE])

  }
  return(exp_polysph(V, mu = fit$mu, d = d))

}


#' @rdname competitors
#' @export
extrinsic_pca <- function(X, d) {

  xbar <- colMeans(X)
  s <- svd(sweep(X, 2, xbar))
  return(list(xbar = xbar, rot = s$v, scores = s$u %*% diag(s$d),
              var = 100 * s$d^2 / sum(s$d^2)))

}


#' @rdname competitors
#' @export
reconstruct_extrinsic_pca <- function(fit, q, d) {

  Z <- matrix(fit$xbar, nrow = nrow(fit$scores), ncol = length(fit$xbar),
              byrow = TRUE)
  if (q > 0) {

    Z <- Z + fit$scores[, 1:q, drop = FALSE] %*% t(fit$rot[, 1:q, drop = FALSE])

  }
  return(proj_polysph(Z, d = d))

}


#' @rdname competitors
#' @param sphere_type type of marginal PNS fits, see \code{\link[shapes]{pns}}.
#' @export
cpns <- function(X, d, sphere_type = "small") {

  ind <- block_ind(d)
  marg <- vector("list", length(d))
  E <- NULL
  for (j in seq_along(d)) {

    Xj <- X[, ind[[j]], drop = FALSE]
    if (d[j] == 1) {

      # Circle: centered angles about the circular mean
      th <- atan2(Xj[, 2], Xj[, 1])
      m <- atan2(mean(sin(th)), mean(cos(th)))
      res <- rbind(((th - m + pi) %% (2 * pi)) - pi)
      marg[[j]] <- list(type = "circle", mean = m)

    } else {

      fj <- quiet(shapes::pns(t(Xj), sphere.type = sphere_type,
                              output = FALSE))
      res <- fj$resmat
      marg[[j]] <- list(type = "pns", PNS = fj$PNS)

    }
    E <- cbind(E, t(res))

  }
  ebar <- colMeans(E)
  s <- svd(sweep(E, 2, ebar))
  return(list(marg = marg, ebar = ebar, rot = s$v, scores = s$u %*% diag(s$d),
              var = 100 * s$d^2 / sum(s$d^2)))

}


#' @rdname competitors
#' @export
reconstruct_cpns <- function(fit, q, d) {

  n <- nrow(fit$scores)
  E <- matrix(fit$ebar, nrow = n, ncol = length(fit$ebar), byrow = TRUE)
  if (q > 0) {

    E <- E + fit$scores[, 1:q, drop = FALSE] %*% t(fit$rot[, 1:q, drop = FALSE])

  }
  out <- matrix(0, nrow = n, ncol = sum(d) + length(d))
  ind <- block_ind(d)
  col <- 0
  for (j in seq_along(d)) {

    if (fit$marg[[j]]$type == "circle") {

      th <- E[, col + 1] + fit$marg[[j]]$mean
      out[, ind[[j]]] <- cbind(cos(th), sin(th))
      col <- col + 1

    } else {

      out[, ind[[j]]] <- t(shapes::PNSe2s(t(E[, col + seq_len(d[j]),
                                               drop = FALSE]),
                                          fit$marg[[j]]$PNS))
      col <- col + d[j]

    }

  }
  return(proj_polysph(out, d = d))

}


#' @title Residual curves of competing methods
#'
#' @description Sums of squared polyspherical distances between the sample
#' and its rank-\eqn{q} reconstructions, \eqn{q = 0, \ldots, Q}, for
#' \code{\link{sppca}} and the methods in \code{\link{competitors}}.
#'
#' @inheritParams competitors
#' @param Q maximum rank.
#' @param fit_sppca optional \code{"sppca"} fit (otherwise computed with
#' default arguments).
#' @param methods methods to compare.
#' @param ... further arguments passed to \code{\link{sppca}}.
#' @return A matrix of size \code{c(Q + 1, length(methods))} with the
#' residuals, with attribute \code{"total"} (the residual of the blockwise
#' Fréchet mean) and \code{"time"} (fitting times).
#' @examples
#' set.seed(1)
#' X <- r_circles_polysph(n = 100, r = 3)
#' res <- compare_residuals(X, d = rep(2, 3), Q = 2, p = 4)
#' 1 - res / attr(res, "total")
#' @export
compare_residuals <- function(X, d, Q = 5, fit_sppca = NULL,
                              methods = c("SP-PCA", "Tangent PCA",
                                          "Extrinsic PCA", "CPNS"), ...) {

  res <- matrix(NA, nrow = Q + 1, ncol = length(methods),
                dimnames = list(0:Q, methods))
  times <- setNames(numeric(length(methods)), methods)
  rss <- function(Xh) sum(dist_polysph_pair(Xh, X, d = d, squared = TRUE))
  for (m in methods) {

    t0 <- proc.time()[[3]]
    res[, m] <- switch(m,
      "SP-PCA" = {
        f <- if (is.null(fit_sppca)) sppca(X, d = d, ...) else fit_sppca
        sapply(0:Q, function(q) {
          rss(reconstruct_sppca(f, min(q, ncol(f$scores))))
        })
      },
      "Tangent PCA" = {
        f <- tangent_pca(X, d = d)
        sapply(0:Q, function(q) rss(reconstruct_tangent_pca(f, q, d = d)))
      },
      "Extrinsic PCA" = {
        f <- extrinsic_pca(X, d = d)
        sapply(0:Q, function(q) rss(reconstruct_extrinsic_pca(f, q, d = d)))
      },
      "CPNS" = {
        f <- cpns(X, d = d)
        sapply(0:Q, function(q) rss(reconstruct_cpns(f, q, d = d)))
      })
    times[m] <- proc.time()[[3]] - t0

  }
  mu <- frechet_mean_polysph(X, d = d)
  attr(res, "total") <- rss(matrix(mu, nrow = nrow(X), ncol = length(mu),
                                   byrow = TRUE))
  attr(res, "time") <- times
  return(res)

}
