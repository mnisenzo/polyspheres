
#' @title Scaled Polyspherical Principal Component Analysis
#'
#' @description Dimension reduction for data on the polysphere
#' \eqn{\mathcal{S}^{d_1, \ldots, d_r}} in three steps:
#' \enumerate{
#'   \item \emph{Scaling}: landmark scaled spherical multidimensional
#'   scaling to a hypersphere \eqn{\mathcal{S}^p_\rho} with data-driven
#'   radius \eqn{\rho} (\code{\link{smds_polysph}}).
#'   \item \emph{Nesting}: principal nested spheres (PNS) of the scaled
#'   sample (\code{\link[shapes]{pns}}), giving scores and a backwards
#'   sequence of nested subspheres.
#'   \item \emph{Inversion}: Nadaraya--Watson interpolation from
#'   \eqn{\mathcal{S}^p} to the polysphere (\code{\link{nw_sph_polysph}})
#'   with leave-one-out cross-validated bandwidth, giving modes of variation
#'   and low-rank reconstructions on the original space.
#' }
#'
#' @inheritParams smds_polysph
#' @param sphere_type type of subspheres fitted by PNS, see
#' \code{\link[shapes]{pns}}: \code{"small"} (default), \code{"great"}, or
#' \code{"seq.test"} (small or great subspheres decided by sequential
#' likelihood ratio tests; considerably slower).
#' @param h bandwidth for the inversion. If \code{NULL}, chosen by
#' \code{\link{bw_cv_nw_sph_polysph}}.
#' @return An object of class \code{"sppca"}: a list with
#' \itemize{
#'   \item \code{scores}: \code{c(n, p)} matrix of PNS scores (in the
#'   units of the polyspherical distance, i.e., multiplied by \eqn{\rho}).
#'   \item \code{var_sphere}: percentages of variance explained on
#'   \eqn{\mathcal{S}^p} by each component.
#'   \item \code{smds}, \code{pns}: outputs of the scaling and nesting
#'   steps.
#'   \item \code{h}: inversion bandwidth.
#'   \item \code{X}, \code{d}: the data.
#'   \item \code{time}: elapsed times of each step.
#' }
#' @references
#' Jung, S., Dryden, I. L. and Marron, J. S. (2012). Analysis of principal
#' nested spheres. \emph{Biometrika}, 99(3):551--568.
#'
#' Zoubouloglou, P., García-Portugués, E. and Marron, J. S. (2023). Scaled
#' torus principal component analysis. \emph{Journal of Computational and
#' Graphical Statistics}, 32(3):1024--1035.
#' @examples
#' set.seed(1)
#' X <- r_circles_polysph(n = 150, r = 4, sigma = 0.05)
#' fit <- sppca(X, d = rep(2, 4), p = 4)
#' fit
#' plot(fit$scores[, 1:2], col = rainbow(150), pch = 16)
#' @export
sppca <- function(X, d, p = NULL, rho = NULL, n_landmarks = NULL,
                  landmarks = c("auto", "ward", "kmeans", "random", "all"),
                  sphere_type = "small", h = NULL, maxit = 500) {

  landmarks <- match.arg(landmarks)
  X <- proj_polysph(rbind(X), d = d)

  # Scaling
  smds <- smds_polysph(X = X, d = d, p = p, rho = rho,
                       n_landmarks = n_landmarks, landmarks = landmarks,
                       maxit = maxit)

  # Nesting
  t0 <- proc.time()[[3]]
  pns <- if (sphere_type %in% c("small", "great")) {

    pns_fit(smds$Y, sphere_type = sphere_type)

  } else {

    quiet(shapes::pns(t(smds$Y), sphere.type = sphere_type, output = FALSE))

  }
  t_pns <- proc.time()[[3]] - t0
  scores <- smds$rho * t(pns$resmat)
  colnames(scores) <- paste0("PC", seq_len(ncol(scores)))

  # Inversion bandwidth
  t0 <- proc.time()[[3]]
  if (is.null(h)) {

    h <- bw_cv_nw_sph_polysph(Y = smds$Y, X = X, d = d, rho = smds$rho)$h

  }
  t_inv <- proc.time()[[3]] - t0
  out <- list(scores = scores, var_sphere = pns$percent, smds = smds,
              pns = pns, h = h, X = X, d = d,
              time = c(smds$time, pns = t_pns, bandwidth = t_inv))
  class(out) <- "sppca"
  return(out)

}


#' @export
print.sppca <- function(x, ...) {

  cat("Scaled polyspherical PCA\n")
  cat(sprintf("  n = %d, polysphere dims = (%s), target sphere S^%d\n",
              nrow(x$X), paste(head(x$d, 6), collapse = ", "), x$smds$p))
  cat(sprintf("  radius = %.4f, landmarks = %d, stress = %.4g / %.4g\n",
              x$smds$rho, nrow(x$smds$landmarks$Z), x$smds$stress[1],
              x$smds$stress[2]))
  v <- head(x$var_sphere, 5)
  cat("  % variance on S^p:", paste(sprintf("%.1f", v), collapse = ", "),
      "\n")
  cat(sprintf("  bandwidth = %.4g, time = %.2fs\n", x$h, sum(x$time)))
  invisible(x)

}


#' @title Hyperspherical points from PNS scores
#'
#' @description Maps scores (in polyspherical units) back to
#' \eqn{\mathcal{S}^p} through the fitted nested spheres.
#'
#' @param fit an \code{"sppca"} object.
#' @param scores matrix of size \code{c(m, q)} with \code{q <= p} scores;
#' remaining scores are set to zero.
#' @return A matrix of size \code{c(m, p + 1)}.
#' @examples
#' set.seed(1)
#' X <- r_circles_polysph(n = 100, r = 3)
#' fit <- sppca(X, d = rep(2, 3), p = 3)
#' # The PNS mean on the scaled hypersphere
#' scores_to_sph(fit, scores = 0)
#' @export
scores_to_sph <- function(fit, scores) {

  scores <- rbind(scores)
  p <- ncol(fit$scores)
  full <- matrix(0, nrow = nrow(scores), ncol = p)
  full[, seq_len(ncol(scores))] <- scores
  return(t(shapes::PNSe2s(t(full / fit$smds$rho), fit$pns$PNS)))

}


#' @title Modes of variation and reconstructions of SP-PCA
#'
#' @description \code{mode_sppca} returns the \eqn{j}-th mode of variation
#' evaluated at scores \code{t}, either by inverting the PNS principal
#' curve through \code{\link{nw_sph_polysph}} (\code{method = "curve"}) or
#' by local linear regression of the sample on its scores
#' (\code{method = "regression"}): the estimated conditional extrinsic mean
#' \eqn{E[\boldsymbol{X} \mid \xi_1 = \cdots = \xi_{j-1} = 0, \xi_j = t]},
#' that is, the \eqn{j}-th mode through the center of the earlier ones (see
#' \code{\link{mode_polysph}} and \code{\link{ll_scores_polysph}}).
#' \code{reconstruct_sppca} returns the rank-\eqn{q} reconstruction of the
#' sample, obtained by zeroing the scores beyond the \eqn{q}-th and
#' inverting. \code{residual_sppca} gives the residual sum of squared
#' polyspherical distances for ranks \eqn{q = 0, \ldots, Q}.
#'
#' @param fit an \code{"sppca"} object.
#' @param j component.
#' @param t vector of scores at which to evaluate the mode. Defaults to a
#' grid between the \code{0.02} and \code{0.98} quantiles of the scores.
#' @param method \code{"curve"} or \code{"regression"}.
#' @param q rank of the reconstruction.
#' @param Q maximum rank in \code{residual_sppca}.
#' @return Matrices with points on the polysphere (\code{mode_sppca},
#' \code{reconstruct_sppca}) or a vector of residuals.
#' @examples
#' set.seed(1)
#' X <- r_circles_polysph(n = 150, r = 3, theta_range = c(0, 1.5 * pi))
#' d <- rep(2, 3)
#' fit <- sppca(X, d = d, p = 4)
#' m1 <- mode_sppca(fit, j = 1, method = "regression")
#' dim(m1)
#' # Residuals of the rank-q reconstructions, q = 0, 1, 2
#' residual_sppca(fit, Q = 2)
#' @name sppca_modes
NULL


#' @rdname sppca_modes
#' @export
mode_sppca <- function(fit, j = 1, t = NULL,
                       method = c("curve", "regression")) {

  method <- match.arg(method)
  if (is.null(t)) {

    t <- seq(quantile(fit$scores[, j], 0.02), quantile(fit$scores[, j], 0.98),
             l = 50)

  }
  if (method == "curve") {

    sc <- matrix(0, nrow = length(t), ncol = j)
    sc[, j] <- t
    y <- scores_to_sph(fit, sc)
    out <- nw_sph_polysph(y = y, Y = fit$smds$Y, X = fit$X, d = fit$d,
                          h = fit$h, rho = fit$smds$rho)

  } else {

    if (j == 1) {

      out <- mode_polysph(t = t, scores = fit$scores[, 1], X = fit$X,
                          d = fit$d)

    } else {

      # Condition on the first j scores, with the first j - 1 at zero:
      # E[X | xi_1 = ... = xi_{j-1} = 0, xi_j = t]
      S <- fit$scores[, seq_len(j), drop = FALSE]
      h_grid <- 10^seq(-0.7, 0.5, l = 13)
      cv <- vapply(h_grid, function(h) {
        f <- ll_scores_polysph(S, S, fit$X, d = fit$d, h = h, loo = TRUE)
        sum(dist_polysph_pair(f, fit$X, d = fit$d, squared = TRUE))
      }, numeric(1))
      s_eval <- cbind(matrix(0, nrow = length(t), ncol = j - 1), t)
      out <- ll_scores_polysph(s_eval, S, fit$X, d = fit$d,
                               h = h_grid[which.min(cv)])

    }

  }
  attr(out, "t") <- t
  return(out)

}


#' @rdname sppca_modes
#' @export
reconstruct_sppca <- function(fit, q) {

  if (q == 0) {

    y <- scores_to_sph(fit, matrix(0, nrow = nrow(fit$X), ncol = 1))

  } else {

    y <- scores_to_sph(fit, fit$scores[, seq_len(q), drop = FALSE])

  }
  return(nw_sph_polysph(y = y, Y = fit$smds$Y, X = fit$X, d = fit$d,
                        h = fit$h, rho = fit$smds$rho))

}


#' @rdname sppca_modes
#' @export
residual_sppca <- function(fit, Q = min(10, ncol(fit$scores))) {

  sapply(0:Q, function(q) {
    sum(dist_polysph_pair(reconstruct_sppca(fit, q), fit$X, d = fit$d,
                          squared = TRUE))
  })

}


#' @keywords internal
quiet <- function(expr) {

  sink(tempfile())
  on.exit(sink())
  invisible(force(expr))

}


#' @title Scores of new observations
#'
#' @description Projects new polyspherical observations onto a fitted
#' SP-PCA: out-of-sample scaled MDS through the landmarks (spherical Gower
#' interpolation refined by L-BFGS), followed by projection onto the
#' fitted nested spheres (\code{\link{pns_s2e}}).
#'
#' @param object an \code{"sppca"} object.
#' @param newdata matrix of size \code{c(m, sum(d) + r)}.
#' @param ... ignored.
#' @return A matrix of size \code{c(m, p)} with scores.
#' @examples
#' set.seed(1)
#' X <- r_circles_polysph(n = 120, r = 3)
#' fit <- sppca(X[1:100, ], d = rep(2, 3), p = 3)
#' predict(fit, newdata = X[101:120, ])[, 1:2]
#' @export
predict.sppca <- function(object, newdata, ...) {

  d <- object$d
  Xn <- proj_polysph(rbind(newdata), d = d)
  lm <- object$smds$landmarks
  rho <- object$smds$rho
  p <- object$smds$p
  # Out-of-sample scaled MDS through the landmarks, then PNS projection
  Delta <- dist_polysph_cross(Xn, lm$Z, d = d)
  Y <- oos_smds(Delta = Delta, T = lm$T, rho = rho)
  sc <- rho * t(pns_s2e(Y, object$pns$PNS))
  colnames(sc) <- colnames(object$scores)
  return(sc)

}
