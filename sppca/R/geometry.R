
#' @title Polyspherical geometry utilities
#'
#' @description Basic operations on the polysphere
#' \eqn{\mathcal{S}^{d_1,\ldots,d_r} := \mathcal{S}^{d_1} \times \cdots
#' \times \mathcal{S}^{d_r} \subset \mathbb{R}^{d_1 + \cdots + d_r + r}}:
#' block indexes, projection, geodesic product distances, extrinsic means,
#' and exponential/logarithmic maps.
#'
#' @param x,y matrices of size \code{c(nx, sum(d) + r)} and
#' \code{c(ny, sum(d) + r)} with points on the polysphere (rows).
#' @param d vector of size \code{r} with the dimensions of the spheres.
#' @param w optional vector of weights of length \code{nrow(x)}.
#' @param mu a vector of size \code{sum(d) + r} with a point on the
#' polysphere (base point of the tangent space).
#' @param v a matrix of size \code{c(nv, sum(d) + r)} with tangent vectors at
#' \code{mu}.
#' @param squared return squared distances? Defaults to \code{FALSE}.
#' @return
#' \itemize{
#'   \item \code{block_ind}: a list of length \code{r} with column indexes.
#'   \item \code{proj_polysph}: the projection of \code{x}.
#'   \item \code{dist_polysph_cross}: a matrix of size \code{c(nx, ny)}.
#'   \item \code{dist_polysph_pair}: a vector of size \code{nx} with
#'   rowwise distances.
#'   \item \code{mean_polysph}: the (weighted) extrinsic mean.
#'   \item \code{log_polysph}, \code{exp_polysph}: tangent vectors and
#'   points on the polysphere.
#' }
#' @details The product distance is
#' \eqn{\delta(\boldsymbol{x}, \boldsymbol{y}) = \{\sum_{j = 1}^r
#' \arccos^2(\boldsymbol{x}_j'\boldsymbol{y}_j)\}^{1/2}}. All operations are
#' vectorized over blocks with \code{tcrossprod}, so the cost of
#' \code{dist_polysph_cross} is \eqn{O(n_x n_y (d_1 + \cdots + d_r + r))}.
#' @examples
#' d <- c(2, 2)
#' x <- polykde::r_unif_polysph(n = 5, d = d)
#' D <- dist_polysph_cross(x, x, d = d)
#' max(abs(D - t(D)))
#' @name polysph_geometry
NULL


#' @rdname polysph_geometry
#' @export
block_ind <- function(d) {

  end <- cumsum(d + 1)
  ini <- c(1, head(end, -1) + 1)
  return(mapply(seq.int, ini, end, SIMPLIFY = FALSE))

}


#' @rdname polysph_geometry
#' @export
proj_polysph <- function(x, d) {

  x <- rbind(x)
  for (ind in block_ind(d)) {

    nrm <- sqrt(rowSums(x[, ind, drop = FALSE]^2))
    nrm[nrm == 0] <- 1
    x[, ind] <- x[, ind, drop = FALSE] / nrm

  }
  return(x)

}


#' @rdname polysph_geometry
#' @export
dist_polysph_cross <- function(x, y, d, squared = FALSE) {

  x <- rbind(x)
  y <- rbind(y)
  D2 <- matrix(0, nrow = nrow(x), ncol = nrow(y))
  for (ind in block_ind(d)) {

    cp <- tcrossprod(x[, ind, drop = FALSE], y[, ind, drop = FALSE])
    D2 <- D2 + acos(pmin(pmax(cp, -1), 1))^2

  }
  if (squared) {

    return(D2)

  } else {

    return(sqrt(D2))

  }

}


#' @rdname polysph_geometry
#' @export
dist_polysph_pair <- function(x, y, d, squared = FALSE) {

  x <- rbind(x)
  y <- rbind(y)
  stopifnot(nrow(x) == nrow(y))
  D2 <- 0
  for (ind in block_ind(d)) {

    cp <- rowSums(x[, ind, drop = FALSE] * y[, ind, drop = FALSE])
    D2 <- D2 + acos(pmin(pmax(cp, -1), 1))^2

  }
  if (squared) {

    return(D2)

  } else {

    return(sqrt(D2))

  }

}


#' @rdname polysph_geometry
#' @export
mean_polysph <- function(x, d, w = NULL) {

  x <- rbind(x)
  if (is.null(w)) {

    w <- rep(1 / nrow(x), nrow(x))

  }
  return(drop(proj_polysph(rbind(colSums(w * x)), d = d)))

}


#' @rdname polysph_geometry
#' @export
log_polysph <- function(x, mu, d) {

  x <- rbind(x)
  v <- x
  for (ind in block_ind(d)) {

    m <- mu[ind]
    cp <- pmin(pmax(drop(x[, ind, drop = FALSE] %*% m), -1), 1)
    th <- acos(cp)
    u <- x[, ind, drop = FALSE] - outer(cp, m)
    nu <- sqrt(rowSums(u^2))
    fac <- ifelse(nu > 1e-12, th / nu, 1)
    v[, ind] <- u * fac

  }
  return(v)

}


#' @rdname polysph_geometry
#' @export
exp_polysph <- function(v, mu, d) {

  v <- rbind(v)
  x <- v
  for (ind in block_ind(d)) {

    vv <- v[, ind, drop = FALSE]
    # Remove any normal component (numerical safety)
    vv <- vv - outer(drop(vv %*% mu[ind]), mu[ind])
    nv <- sqrt(rowSums(vv^2))
    fac <- ifelse(nv > 1e-12, sin(nv) / nv, 1)
    x[, ind] <- outer(cos(nv), mu[ind]) + vv * fac

  }
  return(proj_polysph(x, d = d))

}


#' @title Angular distances on the hypersphere
#'
#' @description Geodesic distances between the rows of two matrices,
#' after normalizing them to the unit sphere.
#'
#' @param x,y matrices with points (not necessarily normalized) as rows.
#' @return A matrix of size \code{c(nrow(x), nrow(y))}.
#' @keywords internal
dist_sph_cross <- function(x, y) {

  x <- rbind(x) / sqrt(rowSums(rbind(x)^2))
  y <- rbind(y) / sqrt(rowSums(rbind(y)^2))
  return(acos(pmin(pmax(tcrossprod(x, y), -1), 1)))

}
