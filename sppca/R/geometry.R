
#' @title Polyspherical geometry utilities
#'
#' @description Basic operations on the polysphere
#' \eqn{\mathcal{S}^{d_1,\ldots,d_r} := \mathcal{S}^{d_1} \times \cdots
#' \times \mathcal{S}^{d_r} \subset \mathbb{R}^{d_1 + \cdots + d_r + r}}:
#' block indexes, projection, product geodesic distances, extrinsic means,
#' and exponential and logarithmic maps.
#'
#' @param x,y matrices of size \code{c(nx, sum(d) + r)} and
#' \code{c(ny, sum(d) + r)} with points on the polysphere as rows. Vectors
#' are interpreted as a single point.
#' @param d vector of size \code{r} with the dimensions of the spheres.
#' @param w optional vector of weights of length \code{nrow(x)}. Defaults to
#' equal weights.
#' @param mu a vector of size \code{sum(d) + r} with a point on the
#' polysphere, the base point of the tangent space.
#' @param v a matrix of size \code{c(nv, sum(d) + r)} with tangent vectors at
#' \code{mu}.
#' @param squared return squared distances? Defaults to \code{FALSE}.
#' @return
#' \itemize{
#'   \item \code{block_ind}: a list of length \code{r} with the column
#'   indexes of each sphere.
#'   \item \code{proj_polysph}: the projection of \code{x} onto the
#'   polysphere.
#'   \item \code{dist_polysph_cross}: a matrix of size \code{c(nx, ny)} with
#'   the distances between the rows of \code{x} and \code{y}.
#'   \item \code{dist_polysph_pair}: a vector of size \code{nx} with the
#'   rowwise distances between \code{x} and \code{y}.
#'   \item \code{mean_polysph}: the (weighted) extrinsic mean.
#'   \item \code{log_polysph}, \code{exp_polysph}: tangent vectors at
#'   \code{mu} and points on the polysphere.
#' }
#' @details The product geodesic distance is
#' \deqn{\delta(\boldsymbol{x}, \boldsymbol{y}) = \Big\{\sum_{j = 1}^r
#' \arccos^2(\boldsymbol{x}_j'\boldsymbol{y}_j)\Big\}^{1/2}.}
#' Rowwise blockwise sums are computed with a single product against the
#' \eqn{(d_1 + \cdots + d_r + r) \times r} block-indicator matrix, so
#' \code{proj_polysph} and \code{dist_polysph_pair} have no loops over the
#' spheres. \code{dist_polysph_cross} costs \eqn{O(n_x n_y (d_1 + \cdots +
#' d_r + r))} through one \code{tcrossprod} per sphere.
#' @examples
#' d <- c(2, 2)
#' x <- polykde::r_unif_polysph(n = 5, d = d)
#' D <- dist_polysph_cross(x, x, d = d)
#' max(abs(D - t(D)))
#' all.equal(diag(dist_polysph_cross(x, x[5:1, ], d = d)),
#'           dist_polysph_pair(x, x[5:1, ], d = d))
#' @name polysph_geometry
NULL


#' @rdname polysph_geometry
#' @export
block_ind <- function(d) {

  end <- cumsum(d + 1)
  ini <- c(1, head(end, -1) + 1)
  return(mapply(seq.int, ini, end, SIMPLIFY = FALSE))

}


#' @title Block-indicator matrix
#'
#' @description Matrix \eqn{\boldsymbol{B}} of size
#' \code{c(sum(d) + r, r)} with \eqn{B_{kj} = 1} if column \eqn{k} belongs to
#' the \eqn{j}-th sphere, so that \code{(x * x) \%*\% B} gives the blockwise
#' squared norms of the rows of \code{x}.
#'
#' @inheritParams polysph_geometry
#' @return A numeric matrix.
#' @keywords internal
block_mat <- function(d) {

  B <- matrix(0, nrow = sum(d + 1), ncol = length(d))
  B[cbind(seq_len(sum(d + 1)), rep(seq_along(d), d + 1))] <- 1
  return(B)

}


#' @title Clamp to \eqn{[-1, 1]}
#'
#' @description In-place clamp used before \code{acos}, cheaper than
#' \code{pmin(pmax(x, -1), 1)}.
#'
#' @param x a numeric vector or matrix.
#' @return \code{x} with its entries clamped to \eqn{[-1, 1]}.
#' @keywords internal
clamp1 <- function(x) {

  x[x > 1] <- 1
  x[x < -1] <- -1
  return(x)

}


#' @rdname polysph_geometry
#' @export
proj_polysph <- function(x, d) {

  if (!is.matrix(x)) {

    x <- rbind(x)

  }

  # Blockwise norms, expanded back to the columns of x
  blk <- rep(seq_along(d), d + 1)
  nrm <- sqrt((x * x) %*% block_mat(d))
  nrm[nrm == 0] <- 1
  return(x / nrm[, blk, drop = FALSE])

}


#' @rdname polysph_geometry
#' @export
dist_polysph_cross <- function(x, y, d, squared = FALSE) {

  if (!is.matrix(x)) x <- rbind(x)
  if (!is.matrix(y)) y <- rbind(y)

  # Sum of squared angles over the spheres
  D2 <- 0
  for (ind in block_ind(d)) {

    D2 <- D2 + acos(clamp1(tcrossprod(x[, ind, drop = FALSE],
                                      y[, ind, drop = FALSE])))^2

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

  if (!is.matrix(x)) x <- rbind(x)
  if (!is.matrix(y)) y <- rbind(y)
  if (nrow(x) != nrow(y)) {

    stop("x and y must have the same number of rows.")

  }

  # Blockwise inner products in one product, then squared angles
  D2 <- rowSums(acos(clamp1((x * y) %*% block_mat(d)))^2)
  if (squared) {

    return(D2)

  } else {

    return(sqrt(D2))

  }

}


#' @rdname polysph_geometry
#' @export
mean_polysph <- function(x, d, w = NULL) {

  if (!is.matrix(x)) x <- rbind(x)
  if (is.null(w)) {

    w <- rep(1 / nrow(x), nrow(x))

  }
  return(drop(proj_polysph(rbind(colSums(w * x)), d = d)))

}


#' @rdname polysph_geometry
#' @export
log_polysph <- function(x, mu, d) {

  if (!is.matrix(x)) x <- rbind(x)
  B <- block_mat(d)
  blk <- rep(seq_along(d), d + 1)

  # Blockwise cosines with mu and the corresponding angles
  cp <- clamp1((x * rep_rows(mu, nrow(x))) %*% B)
  th <- acos(cp)

  # Component of x orthogonal to mu, rescaled to have norm equal to the angle
  u <- x - cp[, blk, drop = FALSE] * rep_rows(mu, nrow(x))
  nu <- sqrt((u * u) %*% B)
  fac <- ifelse(nu > 1e-12, th / nu, 1)
  return(u * fac[, blk, drop = FALSE])

}


#' @rdname polysph_geometry
#' @export
exp_polysph <- function(v, mu, d) {

  if (!is.matrix(v)) v <- rbind(v)
  B <- block_mat(d)
  blk <- rep(seq_along(d), d + 1)
  M <- rep_rows(mu, nrow(v))

  # Remove any normal component (numerical safety)
  v <- v - ((v * M) %*% B)[, blk, drop = FALSE] * M

  # Blockwise exponential map
  nv <- sqrt((v * v) %*% B)
  fac <- ifelse(nv > 1e-12, sin(nv) / nv, 1)
  x <- cos(nv)[, blk, drop = FALSE] * M + v * fac[, blk, drop = FALSE]
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

  if (!is.matrix(x)) x <- rbind(x)
  if (!is.matrix(y)) y <- rbind(y)
  x <- x / sqrt(rowSums(x * x))
  y <- y / sqrt(rowSums(y * y))
  return(acos(clamp1(tcrossprod(x, y))))

}


#' @title Internal helpers
#'
#' @description \code{row_min} returns the rowwise minima of a matrix using
#' the C-level \code{\link[base]{max.col}} (much faster than
#' \code{apply(M, 1, min)}); \code{rep_rows} stacks a vector \code{m} times.
#'
#' @param M a numeric matrix.
#' @param v a vector.
#' @param m number of rows.
#' @keywords internal
row_min <- function(M) {

  M[cbind(seq_len(nrow(M)), max.col(-M, ties.method = "first"))]

}


#' @rdname row_min
#' @keywords internal
rep_rows <- function(v, m) {

  matrix(v, nrow = m, ncol = length(v), byrow = TRUE)

}
