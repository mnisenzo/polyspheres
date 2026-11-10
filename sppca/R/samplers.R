
#' @title Samplers of structured polyspherical data
#'
#' @description \code{r_circles_polysph} simulates
#' \eqn{\boldsymbol{x}_i = (\boldsymbol{x}_{i1}, \ldots, \boldsymbol{x}_{ir})
#' \in (\mathcal{S}^2)^r} where the first \code{n_active} components move
#' along (great or small) circles at a common rate driven by a single latent
#' angle \eqn{\theta_i}, and the remaining components are concentrated von
#' Mises--Fisher noise about fixed directions. A one-dimensional mode of
#' variation therefore underlies the data. \code{r_clusters_polysph}
#' simulates a mixture of products of von Mises--Fisher distributions.
#'
#' @param n sample size.
#' @param r number of spheres \eqn{\mathcal{S}^2}.
#' @param sigma standard deviation of the tangent Gaussian noise added to
#' the circles.
#' @param n_active number of spheres following the common circular path.
#' @param kappa concentration of the inactive components.
#' @param height height of the circles' planes (\code{0} gives great
#' circles, values in \code{(0, 1)} small circles).
#' @param theta_range range of the latent angle.
#' @param centers list of \code{K} matrices of size \code{c(r, 3)} with the
#' cluster centers (for \code{r_clusters_polysph}); if \code{NULL}, random.
#' @param K number of clusters.
#' @return A matrix of size \code{c(n, 3 r)} with attributes \code{"theta"}
#' (latent angle) or \code{"cluster"} (labels).
#' @examples
#' X <- r_circles_polysph(n = 100, r = 3)
#' dim(X)
#' @name samplers
NULL


#' @rdname samplers
#' @export
r_circles_polysph <- function(n, r, sigma = 0.05, n_active = r,
                              kappa = 200, height = 0,
                              theta_range = c(0, 2 * pi)) {

  theta <- sort(runif(n, theta_range[1], theta_range[2]))
  X <- matrix(0, nrow = n, ncol = 3 * r)
  for (j in seq_len(r)) {

    # Random orthonormal frame for each sphere
    Q <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
    if (j <= n_active) {

      a <- sqrt(1 - height^2)
      P <- cbind(a * cos(theta), a * sin(theta), height) %*% t(Q)
      P <- P + sigma * matrix(rnorm(3 * n), n, 3)

    } else {

      P <- rotasym::r_vMF(n = n, mu = Q[, 3], kappa = kappa)

    }
    X[, 3 * (j - 1) + 1:3] <- P / sqrt(rowSums(P^2))

  }
  attr(X, "theta") <- theta
  return(X)

}


#' @rdname samplers
#' @export
r_clusters_polysph <- function(n, r, K = 3, kappa = 50, centers = NULL) {

  if (is.null(centers)) {

    centers <- replicate(K, rotasym::r_unif_sphere(n = r, p = 3),
                         simplify = FALSE)

  }
  cl <- sort(sample(K, size = n, replace = TRUE))
  X <- matrix(0, nrow = n, ncol = 3 * r)
  for (k in unique(cl)) {

    ik <- which(cl == k)
    for (j in seq_len(r)) {

      X[ik, 3 * (j - 1) + 1:3] <- rotasym::r_vMF(n = length(ik),
                                                 mu = centers[[k]][j, ],
                                                 kappa = kappa)

    }

  }
  attr(X, "cluster") <- cl
  return(X)

}
