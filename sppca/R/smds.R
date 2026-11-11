
#' @title Landmarks on the polysphere
#'
#' @description Selects \eqn{\ell} landmarks summarizing a polyspherical
#' sample by Ward clustering on the product geodesic distance (default),
#' \eqn{k}-means on the extrinsic coordinates (for large \eqn{n}), or random
#' subsampling. Cluster landmarks are the projected extrinsic means of the
#' clusters.
#'
#' @inheritParams polysph_geometry
#' @param X a matrix of size \code{c(n, sum(d) + r)} with the sample.
#' @param n_landmarks number of landmarks \eqn{\ell}. For
#' \code{smds_polysph}, defaults to \code{min(n, max(50, min(300, n / 5)))}.
#' @param method one of \code{"ward"}, \code{"kmeans"}, \code{"random"}, or
#' \code{"all"} (all the data, i.e., full scaled MDS).
#' @param D optional precomputed \code{n x n} distance matrix (only used by
#' \code{"ward"}).
#' @return A list with \code{Z} (\code{n_landmarks x (sum(d) + r)} matrix
#' of landmarks) and \code{cluster} (cluster labels of the sample).
#' @examples
#' X <- polykde::r_path_s2r(n = 200, r = 3, sigma = 0.1)
#' X <- matrix(aperm(X, c(1, 3, 2)), nrow = 200)
#' str(landmarks_polysph(X, d = rep(2, 3), n_landmarks = 20))
#' @export
landmarks_polysph <- function(X, d, n_landmarks,
                              method = c("ward", "kmeans", "random", "all"),
                              D = NULL) {

  method <- match.arg(method)
  n <- nrow(X)
  n_landmarks <- min(n_landmarks, n)
  if (method == "all" || n_landmarks == n) {

    return(list(Z = X, cluster = seq_len(n)))

  }
  cluster <- switch(method,
    ward = {
      if (is.null(D)) {

        D <- dist_polysph_cross(X, X, d = d)

      }
      cutree(hclust(as.dist(D), method = "ward.D2"), k = n_landmarks)
    },
    kmeans = {
      kmeans(X, centers = n_landmarks, nstart = 1, iter.max = 50,
             algorithm = "Hartigan-Wong")$cluster
    },
    random = {
      cl <- integer(n)
      cl[sample(n, n_landmarks)] <- seq_len(n_landmarks)
      cl
    })
  if (method == "random") {

    Z <- X[cluster > 0, , drop = FALSE][order(cluster[cluster > 0]), ,
                                        drop = FALSE]

  } else {

    Z <- rowsum(X, group = cluster, reorder = TRUE)
    Z <- proj_polysph(Z, d = d)

  }
  return(list(Z = Z, cluster = cluster))

}


#' @title Spherical stress, its gradient, and profiled radius
#'
#' @description Raw stress \eqn{\sum_{i, j} (\delta_{ij} - \rho
#' \theta_{ij})^2}, where \eqn{\theta_{ij}} is the angle between
#' \eqn{\boldsymbol{y}_i} and \eqn{\boldsymbol{t}_j}, and its gradient with
#' respect to the (unnormalized) \eqn{\boldsymbol{y}_i}'s. For a given
#' configuration the optimal radius is available in closed form,
#' \eqn{\hat\rho = \sum \delta_{ij}\theta_{ij} / \sum \theta_{ij}^2}, so the
#' radius is profiled out (by the envelope theorem the gradient with
#' respect to the configuration is unchanged).
#'
#' @param Y matrix \code{c(n, p + 1)} (unnormalized) to optimize.
#' @param U matrix \code{c(m, p + 1)} of unit vectors (fixed landmarks), or
#' \code{NULL} to use \code{Y} itself (landmark step).
#' @param Delta target distances, matrix \code{c(n, m)} (or \code{c(n, n)}).
#' @param rho radius; if \code{NULL} it is profiled.
#' @return A list with \code{value}, \code{gradient}, and \code{rho}.
#' @keywords internal
stress_sph <- function(Y, U = NULL, Delta, rho = NULL) {

  nY <- sqrt(rowSums(Y^2))
  Yh <- Y / nY
  self <- is.null(U)
  if (self) {

    U <- Yh

  }
  C <- clamp1(tcrossprod(Yh, U))
  Th <- acos(C)
  if (self) {

    diag(Th) <- 0

  }
  if (is.null(rho)) {

    rho <- sum(Delta * Th) / max(sum(Th^2), .Machine$double.eps)

  }
  R <- Delta - rho * Th
  if (self) {

    diag(R) <- 0

  }
  value <- sum(R^2)
  # Floor sin(theta) away from zero (coincident points have zero gradient)
  Sn <- sqrt(1 - C * C)
  Sn[Sn < 1e-10] <- 1e-10
  A <- 2 * rho * R / Sn
  if (self) {

    diag(A) <- 0
    # Each pair appears twice in the full matrix: both rows receive the
    # derivative of their own terms, which is the gradient of sum_{i != j}
    A <- 2 * A

  }
  G <- (A %*% U - rowSums(A * C) * Yh) / nY
  return(list(value = value, gradient = G, rho = rho))

}


#' @title Spherical classical MDS initialization
#'
#' @description Initial configuration on \eqn{\mathcal{S}^p} for a target
#' distance matrix \eqn{\Delta}: for a candidate radius \eqn{\rho}, the Gram
#' matrix of a configuration on \eqn{\mathcal{S}^p_\rho} is
#' \eqn{\cos(\Delta / \rho)}; its leading \eqn{p + 1} eigenvectors give the
#' configuration. The radius is chosen among a small grid by the resulting
#' (profiled) stress.
#'
#' @param Delta symmetric matrix of target distances.
#' @param p dimension of the target sphere.
#' @param rho_grid optional grid of candidate radii.
#' @return A list with \code{Y} (unit rows) and \code{rho}.
#' @keywords internal
init_smds <- function(Delta, p, rho_grid = NULL) {

  if (is.null(rho_grid)) {

    rho_grid <- max(Delta) / pi * c(1, 1.1, 1.25, 1.5, 2, 3, 5)

  }
  best <- list(value = Inf)
  for (rho in rho_grid) {

    G <- cos(pmin(Delta / rho, pi))
    e <- eigen(G, symmetric = TRUE)
    Y <- e$vectors[, 1:(p + 1), drop = FALSE] %*%
      diag(sqrt(pmax(e$values[1:(p + 1)], 1e-8)), p + 1)
    Y <- Y / sqrt(rowSums(Y^2))
    st <- stress_sph(Y = Y, Delta = Delta)
    if (st$value < best$value) {

      best <- list(value = st$value, Y = Y, rho = st$rho)

    }

  }
  return(best)

}


#' @title Scaled spherical multidimensional scaling with landmarks
#'
#' @description Maps a polyspherical sample \eqn{\boldsymbol{x}_1, \ldots,
#' \boldsymbol{x}_n \in \mathcal{S}^{d_1, \ldots, d_r}} to the hypersphere
#' \eqn{\mathcal{S}^p} with data-driven radius \eqn{\rho} by minimizing the
#' geodesic stress
#' \deqn{\sum_{i, j} \{\delta(\boldsymbol{x}_i, \boldsymbol{\zeta}_j) - \rho
#' \arccos(\boldsymbol{y}_i'\boldsymbol{t}_j)\}^2,}
#' where \eqn{\boldsymbol{\zeta}_1, \ldots, \boldsymbol{\zeta}_\ell} are
#' landmarks. Proceeds in two steps: (1) scaled spherical MDS of the
#' landmarks (radius and landmark images \eqn{\boldsymbol{t}_j} estimated
#' jointly); (2) out-of-sample placement of all observations given the
#' landmark images. Step (2) is initialized by a closed-form spherical
#' analogue of Gower's interpolation,
#' \eqn{\boldsymbol{y}_i \propto (\boldsymbol{T}'\boldsymbol{T})^{-1}
#' \boldsymbol{T}' \cos(\boldsymbol{\delta}_i / \rho)}, and refined by
#' pointwise Levenberg--Marquardt iterations vectorized over the
#' observations (see \code{\link{oos_smds}}). Step (1) uses L-BFGS-B with
#' analytic gradients and the radius profiled out in closed form, from a
#' spherical classical-MDS initialization.
#'
#' @inheritParams landmarks_polysph
#' @param p dimension of the target sphere \eqn{\mathcal{S}^p}. Defaults to
#' \code{min(sum(d), n_landmarks - 1)}.
#' @param rho radius of the target sphere. If \code{NULL} (default), it is
#' estimated.
#' @param landmarks landmark selection method, see
#' \code{\link{landmarks_polysph}}. The default, \code{"auto"}, uses Ward
#' clustering for \code{n <= 1000} and \eqn{k}-means otherwise.
#' @param maxit maximum number of L-BFGS iterations for each step. Defaults
#' to \code{500}.
#' @param tol relative reduction of the stress below which the L-BFGS
#' iterations of the landmark step stop (passed as
#' \code{factr = tol / .Machine$double.eps} to
#' \code{\link[stats]{optim}}). Defaults to \code{1e-8}.
#' @param refine refine the Gower initialization of step (2)? Defaults to
#' \code{TRUE}.
#' @return A list with:
#' \itemize{
#'   \item \code{Y}: \code{c(n, p + 1)} matrix of unit vectors.
#'   \item \code{rho}: estimated radius.
#'   \item \code{landmarks}: output of \code{\link{landmarks_polysph}} plus
#'   \code{T}, the landmark images on \eqn{\mathcal{S}^p}.
#'   \item \code{stress}: normalized stress (stress-1 squared) of the
#'   landmark and out-of-sample steps.
#'   \item \code{time}: elapsed times.
#' }
#' @examples
#' X <- polykde::r_path_s2r(n = 200, r = 4, sigma = 0.1)
#' X <- matrix(aperm(X, c(1, 3, 2)), nrow = 200)
#' fit <- smds_polysph(X, d = rep(2, 4), p = 3, n_landmarks = 40)
#' fit$rho
#' fit$stress
#' @export
smds_polysph <- function(X, d, p = NULL, rho = NULL, n_landmarks = NULL,
                         landmarks = c("auto", "ward", "kmeans", "random",
                                       "all"),
                         maxit = 500, tol = 1e-8, refine = TRUE) {

  # Check inputs
  X <- rbind(X)
  stopifnot(ncol(X) == sum(d) + length(d))
  landmarks <- match.arg(landmarks)
  n <- nrow(X)
  if (landmarks == "auto") {

    # Ward needs the n x n distance matrix; switch to k-means for large n
    landmarks <- if (n <= 1000) "ward" else "kmeans"

  }
  if (is.null(n_landmarks)) {

    n_landmarks <- min(n, max(50, min(300, ceiling(n / 5))))

  }
  if (is.null(p)) {

    p <- min(sum(d), n_landmarks - 1)

  }
  p <- min(p, n_landmarks - 1)
  times <- numeric(3)
  names(times) <- c("landmarks", "smds_landmarks", "out_of_sample")

  # Landmarks
  t0 <- proc.time()[[3]]
  lm <- landmarks_polysph(X = X, d = d, n_landmarks = n_landmarks,
                          method = landmarks)
  Z <- lm$Z
  ell <- nrow(Z)
  Delta_L <- dist_polysph_cross(Z, Z, d = d)
  times[1] <- proc.time()[[3]] - t0

  # Step 1: scaled spherical MDS of the landmarks
  t0 <- proc.time()[[3]]
  ini <- init_smds(Delta = Delta_L, p = p)
  fn_cache <- new.env()
  fn <- function(y) {
    st <- stress_sph(Y = matrix(y, nrow = ell), Delta = Delta_L, rho = rho)
    assign("gr", c(st$gradient), envir = fn_cache)
    st$value
  }
  gr <- function(y) {
    get("gr", envir = fn_cache)
  }
  opt <- optim(par = c(ini$Y), fn = fn, gr = gr, method = "L-BFGS-B",
               control = list(maxit = maxit, factr = tol / .Machine$double.eps))
  T_L <- matrix(opt$par, nrow = ell)
  T_L <- T_L / sqrt(rowSums(T_L^2))
  st_L <- stress_sph(Y = T_L, Delta = Delta_L, rho = rho)
  rho_hat <- st_L$rho
  stress_L <- st_L$value / sum(Delta_L[upper.tri(Delta_L)]^2) / 2
  times[2] <- proc.time()[[3]] - t0

  # Step 2: out-of-sample placement of all observations
  t0 <- proc.time()[[3]]
  if (landmarks == "all" || ell == n) {

    Y <- T_L
    stress_O <- stress_L

  } else {

    # Spherical Gower interpolation, refined by Levenberg--Marquardt
    Delta_O <- dist_polysph_cross(X, Z, d = d)
    Y <- oos_smds(Delta = Delta_O, T = T_L, rho = rho_hat,
                  maxit = if (refine) 30 else 0)
    stress_O <- stress_sph(Y = Y, U = T_L, Delta = Delta_O,
                           rho = rho_hat)$value / sum(Delta_O^2)

  }
  times[3] <- proc.time()[[3]] - t0
  lm$T <- T_L
  return(list(Y = Y, rho = rho_hat, p = p, landmarks = lm,
              stress = c(landmarks = stress_L, out_of_sample = stress_O),
              time = times))

}


#' @title Full-sample stress of a spherical configuration
#'
#' @description Normalized stress \eqn{\sum_{i < j} (\delta_{ij} - \rho
#' \theta_{ij})^2 / \sum_{i < j} \delta_{ij}^2} between the polyspherical
#' sample and its image on \eqn{\mathcal{S}^p_\rho}. Costs
#' \eqn{O(n^2)}; used for evaluation only.
#'
#' @inheritParams smds_polysph
#' @param Y \code{c(n, p + 1)} matrix of unit vectors.
#' @param rho radius. If \code{NULL}, the optimal radius is used.
#' @return The normalized stress (attribute \code{"rho"} gives the radius).
#' @examples
#' set.seed(1)
#' X <- r_circles_polysph(n = 100, r = 3)
#' fit <- smds_polysph(X, d = rep(2, 3), p = 3, n_landmarks = 30)
#' stress_polysph(X, fit$Y, d = rep(2, 3))
#' @export
stress_polysph <- function(X, Y, d, rho = NULL) {

  D <- dist_polysph_cross(X, X, d = d)
  Th <- dist_sph_cross(Y, Y)
  ut <- upper.tri(D)
  if (is.null(rho)) {

    rho <- sum(D[ut] * Th[ut]) / sum(Th[ut]^2)

  }
  out <- sum((D[ut] - rho * Th[ut])^2) / sum(D[ut]^2)
  attr(out, "rho") <- rho
  return(out)

}


#' @title Out-of-sample scaled spherical MDS
#'
#' @description Places \eqn{m} points on \eqn{\mathcal{S}^p} given their
#' polyspherical distances \eqn{\delta_{ij}} to \eqn{\ell} landmarks with
#' fixed images \eqn{\boldsymbol{t}_j} and radius \eqn{\rho}, by solving
#' the \eqn{m} decoupled problems
#' \deqn{\hat{\boldsymbol{y}}_i = \arg\min_{\boldsymbol{y} \in
#' \mathcal{S}^p} \sum_{j = 1}^\ell \{\delta_{ij} - \rho
#' \arccos(\boldsymbol{y}'\boldsymbol{t}_j)\}^2.}
#'
#' @param Delta matrix of size \code{c(m, l)} with the distances to the
#' landmarks.
#' @param T matrix of size \code{c(l, p + 1)} with the landmark images
#' (unit rows).
#' @param rho radius.
#' @param maxit maximum number of Levenberg--Marquardt iterations;
#' \code{0} returns the Gower initialization. Defaults to \code{30}.
#' @param tol relative decrease of each pointwise stress below which that
#' point is considered converged. Defaults to \code{1e-6}.
#' @return A matrix of size \code{c(m, p + 1)} with unit rows.
#' @details The initialization is the closed-form spherical Gower
#' interpolation: since \eqn{\boldsymbol{t}_j'\boldsymbol{y}_i \approx
#' \cos(\delta_{ij} / \rho)} at the optimum, \eqn{\boldsymbol{y}_i^{(0)}
#' \propto (\boldsymbol{T}'\boldsymbol{T})^{-1}\boldsymbol{T}'
#' \cos(\boldsymbol{\delta}_i / \rho)}.
#'
#' Each point is then refined by Levenberg--Marquardt iterations on the
#' sphere. With \eqn{c_j = \boldsymbol{y}'\boldsymbol{t}_j},
#' \eqn{s_j = (1 - c_j^2)^{1/2}}, residuals \eqn{r_j = \delta_{ij} - \rho
#' \arccos c_j}, and \eqn{\boldsymbol{P} = \boldsymbol{I} -
#' \boldsymbol{y}\boldsymbol{y}'}, the Gauss--Newton step in the tangent
#' space solves
#' \deqn{\rho\{\boldsymbol{P}\boldsymbol{A}\boldsymbol{P} +
#' \lambda\boldsymbol{I}\}\boldsymbol{v} = -\boldsymbol{P}\boldsymbol{T}'
#' (\boldsymbol{r} / \boldsymbol{s}), \quad \boldsymbol{A} =
#' \boldsymbol{T}'\mathrm{diag}(\boldsymbol{s}^{-2})\boldsymbol{T},}
#' followed by the retraction \eqn{\boldsymbol{y} \leftarrow
#' (\boldsymbol{y} + \boldsymbol{v}) / \|\boldsymbol{y} + \boldsymbol{v}\|}.
#' The matrices \eqn{\boldsymbol{A}} of all points are obtained at once as
#' the product of the \eqn{m \times \ell} matrix of \eqn{s_{ij}^{-2}} with
#' the \eqn{\ell \times (p + 1)^2} matrix of the \eqn{\boldsymbol{t}_j
#' \boldsymbol{t}_j'}, so only the small \eqn{(p + 1) \times (p + 1)} solves
#' are pointwise. The damping \eqn{\lambda} is adapted pointwise (steps that
#' increase the stress are rejected), so the iterations converge
#' quadratically near the solution.
#' @keywords internal
oos_smds <- function(Delta, T, rho, maxit = 30, tol = 1e-6) {

  m <- nrow(Delta)
  p1 <- ncol(T)

  # Spherical Gower interpolation
  Y <- cos(pmin(Delta / rho, pi)) %*% T %*%
    solve(crossprod(T) + 1e-10 * diag(p1))
  Y <- Y / sqrt(rowSums(Y * Y))
  if (maxit == 0) {

    return(Y)

  }

  # Pointwise stresses
  point_stress <- function(Y, D) {

    rowSums((D - rho * acos(clamp1(tcrossprod(Y, T))))^2)

  }

  # Outer products t_j t_j' as rows of an l x (p + 1)^2 matrix
  TT <- t(apply(T, 1, function(t) c(tcrossprod(t))))
  f <- point_stress(Y, Delta)
  lambda <- rep(1e-3, m)
  active <- rep(TRUE, m)
  for (it in seq_len(maxit)) {

    ia <- which(active)
    if (length(ia) == 0) {

      break

    }
    Ya <- Y[ia, , drop = FALSE]

    # Residuals, sines, and the Gauss--Newton ingredients for active points
    C <- clamp1(tcrossprod(Ya, T))
    S <- sqrt(1 - C * C)
    S[S < 1e-8] <- 1e-8
    R <- Delta[ia, , drop = FALSE] - rho * acos(C)
    A <- (1 / S^2) %*% TT
    b <- (R / S) %*% T
    b <- b - rowSums(b * Ya) * Ya

    # Damped steps, one small solve per point
    V <- matrix(0, nrow = length(ia), ncol = p1)
    for (k in seq_along(ia)) {

      y <- Ya[k, ]
      P <- diag(p1) - tcrossprod(y)
      Ak <- P %*% matrix(A[k, ], p1, p1) %*% P
      sc <- max(mean(diag(Ak)), 1e-12)
      V[k, ] <- -solve(Ak + (lambda[ia[k]] * sc) * diag(p1) + tcrossprod(y),
                       b[k, ]) / rho

    }

    # Retraction, acceptance, and damping update
    Yn <- Ya + V
    Yn <- Yn / sqrt(rowSums(Yn * Yn))
    fn <- point_stress(Yn, Delta[ia, , drop = FALSE])
    acc <- fn < f[ia]
    gain <- (f[ia] - fn) / pmax(f[ia], .Machine$double.xmin)
    Y[ia[acc], ] <- Yn[acc, , drop = FALSE]
    f[ia[acc]] <- fn[acc]
    lambda[ia] <- ifelse(acc, pmax(lambda[ia] / 3, 1e-12), lambda[ia] * 4)

    # A point stops when its relative gain or its step is negligible, or when
    # the damping has grown so much that steps are pure (tiny) gradient steps
    step <- sqrt(rowSums(V * V))
    active[ia] <- !((acc & gain < tol) | step < 1e-9 | lambda[ia] > 1e4)

  }
  return(Y)

}
