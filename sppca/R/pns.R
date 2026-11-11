
#' @title Fast principal nested spheres
#'
#' @description Principal nested spheres (PNS) of Jung et al. (2012) with
#' small or great subspheres, computed with a vectorized intrinsic
#' Levenberg--Marquardt fit of each subsphere. The output has the same
#' structure as \code{\link[shapes]{pns}}, so that
#' \code{\link[shapes]{PNSe2s}} and \code{shapes:::PNSs2e} can be used on it.
#'
#' @param Y matrix of size \code{c(n, p + 1)} with the data on
#' \eqn{\mathcal{S}^p} (unit rows).
#' @param sphere_type \code{"small"} or \code{"great"} subspheres.
#' @return A list with entries \code{resmat} (\code{c(p, n)} matrix of scores,
#' first component first), \code{PNS} (fitted nested spheres, as in
#' \code{\link[shapes]{pns}}), and \code{percent} (percentages of variance
#' explained by each component).
#' @details At the level of \eqn{\mathcal{S}^k}, the subsphere
#' \eqn{\{\boldsymbol{y}: \arccos(\boldsymbol{v}'\boldsymbol{y}) = r\}}
#' minimizes \eqn{\sum_i \{\arccos(\boldsymbol{v}'\boldsymbol{y}_i) - r\}^2}.
#' The radius is profiled as the mean angle (small subspheres) or fixed at
#' \eqn{\pi / 2} (great subspheres), and \eqn{\boldsymbol{v}} is fitted by
#' Levenberg--Marquardt on \eqn{\mathcal{S}^k}. The analytic Jacobian
#' \eqn{-(\boldsymbol{y}_i - c_i \boldsymbol{v}) / s_i},
#' \eqn{c_i = \boldsymbol{v}'\boldsymbol{y}_i}, \eqn{s_i = (1 - c_i^2)^{1/2}},
#' is centered when the radius is profiled. As in
#' \code{\link[shapes]{pns}}, two initializations are used (the last
#' singular vector of the data and the last principal component), and the
#' best fit is kept. The data are then projected onto the subsphere, rotated
#' so that \eqn{\boldsymbol{v}} is the north pole, and rescaled to
#' \eqn{\mathcal{S}^{k - 1}}. On the final circle, the exact Fréchet mean is
#' found among the \eqn{n} candidate points of Hotz and Huckemann (2015).
#' @references
#' Hotz, T. and Huckemann, S. (2015). Intrinsic means on the circle:
#' uniqueness, locus and asymptotics. \emph{Annals of the Institute of
#' Statistical Mathematics}, 67(1):177--193.
#'
#' Jung, S., Dryden, I. L. and Marron, J. S. (2012). Analysis of principal
#' nested spheres. \emph{Biometrika}, 99(3):551--568.
#' @examples
#' set.seed(1)
#' th <- runif(100, 0, 2 * pi)
#' Y <- cbind(0.6 * cos(th), 0.6 * sin(th), 0.8) +
#'   0.05 * matrix(rnorm(300), 100, 3)
#' Y <- Y / sqrt(rowSums(Y^2))
#' fit <- pns_fit(Y, sphere_type = "small")
#' fit$percent
#' fit$PNS$dist  # radius of the fitted small circle, about acos(0.8)
#' @export
pns_fit <- function(Y, sphere_type = c("small", "great")) {

  sphere_type <- match.arg(sphere_type)
  x <- t(Y)
  k <- nrow(x)
  n <- ncol(x)
  d <- k - 1
  if (d < 2) {

    stop("PNS requires data on a sphere of dimension at least 2.")

  }

  # Trivial reduction of a null space (data spanning a proper subspace)
  sv <- svd(x, nu = k)
  maxd <- which(sv$d < 1e-15)[1]
  if (is.na(maxd) || k > n) {

    maxd <- min(k, n) + 1

  }
  nullspdim <- k - maxd + 1
  uu <- sv$u
  resmat <- matrix(NA, d, n)
  orthaxis <- vector("list", d)
  dist <- rep(NA, d - 1)
  cur <- x

  # Helper: project to the subsphere orthogonal complement of v and rescale
  descend <- function(cur, v, kk) {

    ns <- rot_mat(v) %*% cur
    ns[seq_len(kk - 1), , drop = FALSE] /
      rep_rows(sqrt(pmax(1 - ns[kk, ]^2, 1e-300)), kk - 1)

  }
  if (nullspdim > 0) {

    for (i in seq_len(nullspdim)) {

      v <- uu[, ncol(uu) - i + 1]
      resmat[i, ] <- acos(clamp1(drop(v %*% cur))) - pi / 2
      orthaxis[[i]] <- v
      dist[i] <- pi / 2
      cur <- descend(cur, v, k - i + 1)
      uu <- descend(uu, v, k - i + 1)

    }

  }

  # Backwards sequence of best-fitting subspheres
  if (nullspdim + 1 <= d - 1) {

    for (i in (nullspdim + 1):(d - 1)) {

      sp <- fit_subsphere(t(cur), small = sphere_type == "small")
      resmat[i, ] <- acos(clamp1(drop(sp$center %*% cur))) - sp$r
      orthaxis[[i]] <- sp$center
      dist[i] <- sp$r
      cur <- descend(cur, sp$center, k - i + 1)

    }

  }

  # Final circle: exact Frechet mean and centered angles
  th <- atan2(cur[2, ], cur[1, ])
  m <- frechet_mean_circle(th)
  orthaxis[[d]] <- m
  resmat[d, ] <- ((th - m + pi) %% (2 * pi)) - pi

  # Scale the residuals by the radii of the nested spheres (as in shapes)
  radii <- c(1, cumprod(sin(dist[seq_len(d - 1)])))
  resmat <- (radii * resmat)[d:1, , drop = FALSE]
  PNS <- list(radii = radii, orthaxis = orthaxis, dist = dist,
              pvalues = NaN, ratio = rep(NA, d - 1), basisu = NULL,
              sphere.type = sphere_type)
  PNS$scores <- t(resmat)
  PNS$mean <- c(shapes::PNSe2s(matrix(0, d, 1), PNS))
  var_pns <- rowSums(resmat^2) / n
  return(list(resmat = resmat, PNS = PNS,
              percent = 100 * var_pns / sum(var_pns)))

}


#' @title Least-squares subsphere of a sample on the sphere
#'
#' @description Fits \eqn{\{\boldsymbol{y}: \arccos(\boldsymbol{v}'
#' \boldsymbol{y}) = r\}} to a sample by Levenberg--Marquardt on the sphere;
#' see \code{\link{pns_fit}}.
#'
#' @param Y matrix of size \code{c(n, k)} with unit rows.
#' @param small fit a small subsphere (\code{TRUE}) or a great one?
#' @param maxit maximum number of iterations per initialization.
#' @return A list with \code{center} and \code{r} (\code{r <= pi / 2}).
#' @keywords internal
fit_subsphere <- function(Y, small = TRUE, maxit = 100) {

  k <- ncol(Y)
  obj <- function(v) {

    th <- acos(clamp1(drop(Y %*% v)))
    r <- if (small) mean(th) else pi / 2
    list(f = sum((th - r)^2), th = th, r = r)

  }
  lm_fit <- function(v) {

    v <- v / sqrt(sum(v^2))
    cur <- obj(v)
    lambda <- 1e-3
    for (it in seq_len(maxit)) {

      # Tangent Jacobian of the angles (centered if the radius is profiled)
      cc <- cos(cur$th)
      s <- sin(cur$th)
      s[s < 1e-8] <- 1e-8
      J <- -(Y - outer(cc, v)) / s
      if (small) J <- sweep(J, 2, colMeans(J))
      e <- cur$th - cur$r
      A <- crossprod(J)
      b <- crossprod(J, e)
      sc <- max(mean(diag(A)), 1e-12)
      delta <- -drop(solve(A + lambda * sc * diag(k) + tcrossprod(v), b))
      delta <- delta - sum(delta * v) * v
      v_new <- (v + delta) / sqrt(sum((v + delta)^2))
      new <- obj(v_new)
      if (new$f < cur$f) {

        gain <- (cur$f - new$f) / max(cur$f, .Machine$double.xmin)
        v <- v_new
        cur <- new
        lambda <- max(lambda / 3, 1e-12)
        if (gain < 1e-12) break

      } else {

        lambda <- lambda * 4
        if (lambda > 1e8) break

      }

    }
    list(center = v, r = cur$r, f = cur$f)

  }

  # Two initializations, as in shapes::pns: last singular vector of the data
  # and last principal component
  v1 <- svd(Y, nu = 0)$v[, k]
  v2 <- prcomp(Y)$rotation[, k]
  f1 <- lm_fit(v1)
  f2 <- lm_fit(v2)
  best <- if (f1$f <= f2$f) f1 else f2
  if (best$r > pi / 2) {

    best$center <- -best$center
    best$r <- pi - best$r

  }
  return(list(center = best$center, r = best$r))

}


#' @title Rotation taking a unit vector to the north pole
#'
#' @description Rotation matrix \eqn{\boldsymbol{R}} with \eqn{\boldsymbol{R}
#' \boldsymbol{b} = (0, \ldots, 0, 1)'} acting in the plane spanned by
#' \eqn{\boldsymbol{b}} and the north pole, as in \code{shapes:::rotMat}.
#'
#' @param b a vector.
#' @return A square matrix.
#' @keywords internal
rot_mat <- function(b) {

  dd <- length(b)
  b <- b / sqrt(sum(b^2))
  a <- c(rep(0, dd - 1), 1)
  ab <- sum(a * b)
  if (abs(ab - 1) < 1e-15) return(diag(dd))
  if (abs(ab + 1) < 1e-15) return(-diag(dd))
  alpha <- acos(ab)
  cv <- b - a * ab
  cv <- cv / sqrt(sum(cv^2))
  A <- tcrossprod(a, cv) - tcrossprod(cv, a)
  return(diag(dd) + sin(alpha) * A +
           (cos(alpha) - 1) * (tcrossprod(a) + tcrossprod(cv)))

}


#' @title Fréchet mean on the circle
#'
#' @description Exact intrinsic (Fréchet) mean of angles: the minimizer of
#' \eqn{\sum_i d^2(\theta_i, m)} is among the \eqn{n} points
#' \eqn{\bar\theta + 2\pi j / n}, \eqn{j = 0, \ldots, n - 1}, with
#' \eqn{\bar\theta} the arithmetic mean of the angles in \eqn{[0, 2\pi)}
#' (Hotz and Huckemann, 2015).
#'
#' @param th vector of angles.
#' @return An angle in \eqn{(-\pi, \pi]}.
#' @keywords internal
frechet_mean_circle <- function(th) {

  th <- th %% (2 * pi)
  n <- length(th)
  cand <- (mean(th) + 2 * pi * (0:(n - 1)) / n) %% (2 * pi)
  # Squared geodesic distances from all candidates (n x n, chunked if large)
  cost <- vapply(cand, function(m) {
    dd <- abs(th - m)
    sum(pmin(dd, 2 * pi - dd)^2)
  }, numeric(1))
  m <- cand[which.min(cost)]
  return(atan2(sin(m), cos(m)))

}


#' @title Scores of new points on fitted principal nested spheres
#'
#' @description Projects points of \eqn{\mathcal{S}^p} onto a fitted PNS
#' decomposition and returns their scores, scaled as in
#' \code{\link{pns_fit}} (first component first). It is the inverse of
#' \code{\link[shapes]{PNSe2s}} on the fitted nested spheres.
#'
#' @param Y matrix of size \code{c(m, p + 1)} with unit rows.
#' @param PNS the \code{PNS} element of the output of \code{\link{pns_fit}}
#' or \code{\link[shapes]{pns}}.
#' @return A matrix of size \code{c(p, m)} with the scores.
#' @examples
#' set.seed(1)
#' Y <- sweep(matrix(rnorm(150), 50, 3), 2, c(0, 0, 3), "+")
#' Y <- Y / sqrt(rowSums(Y^2))
#' fit <- pns_fit(Y)
#' max(abs(pns_s2e(Y, fit$PNS) - fit$resmat))
#' @export
pns_s2e <- function(Y, PNS) {

  cur <- t(rbind(Y))
  kk <- nrow(cur)
  if (kk != length(PNS$mean)) {

    stop("The dimension of Y does not match the PNS decomposition.")

  }
  if (!is.null(PNS$basisu)) {

    cur <- crossprod(PNS$basisu, cur)

  }
  res <- matrix(0, nrow = kk - 1, ncol = ncol(cur))

  # Signed distances to each nested subsphere, then descend
  for (i in seq_len(kk - 2)) {

    v <- PNS$orthaxis[[i]]
    res[i, ] <- acos(clamp1(drop(v %*% cur))) - PNS$dist[i]
    ns <- rot_mat(v) %*% cur
    cur <- ns[seq_len(kk - i), , drop = FALSE] /
      rep_rows(sqrt(pmax(1 - ns[kk - i + 1, ]^2, 1e-300)), kk - i)

  }

  # Angles on the final circle, centered at its Frechet mean
  th <- atan2(cur[2, ], cur[1, ])
  res[kk - 1, ] <- ((th - PNS$orthaxis[[kk - 1]] + pi) %% (2 * pi)) - pi
  return((PNS$radii * res)[(kk - 1):1, , drop = FALSE])

}
