library(rotasym)
library(polykde)
library(waddR)
library(sphunif)
library(shapes)
library(smacof)
library(sdetorus)
library(rgl)
library(MASS)
library(sdetorus)
library(rgl)
library(plot3D)
library(testthat)
library(usethis)
source("Script1_Setup_and_Functions.R")

# Compare centroid and nearest neighbors on the same clustering
compare_centers <-
  function(x,
           angular = TRUE,
           dims = rep(1, ncol(x)),
           landmark_size = .1,
           method = "ward.D2") {
    n <- nrow(x)
    d <- ncol(x)
    d_sphere <- sum(dims)
      # Compute polysphere distances [s^d_1, s^d_2...]
    ind_dj <- comp_ind_dj(dims)

    # Convert from hyperspherical to cartesian coordinates if needed
    cart_x <- x
    if (angular) {
      cart_x <- polyAnglesToSphere(x, dims)
    }
    # Go straight to upper/lower part
    dist <- dist_polysph_matrix(cart_x, ind_dj = ind_dj, std = FALSE)
    # Make sure count of landmarks is integral
    ls <- floor(n * landmark_size)
    # Divide data to clusters
    hc <- hclust(as.dist(dist), method = method)
    ct <- cutree(hc, k = ls)
    # Pick the point closest to the center of each cluster as a landmark
    landmarks <- numeric(ls)
    for (i in 1:ls) {
      grpi <- x[ct == i, ]
      # Convert to data frame for the dist_to_centroids function
      if (is.vector(grpi)) {
        landmarks[i] <- which(ct == i)
      } else {
        dgrpi <- data.frame(grpi, Group = rep(1, nrow(grpi)))
        grp_dists <- dist_to_centroids(dist_subset(dist, which(ct == i)), dgrpi$Group)
        landmarks[i] <- strtoi(grp_dists$Item[which.min(grp_dists$CentroidDistance)])
      }
    }
    landmarks <- cart_x[landmarks,]
    landmarks_c <- matrix(0, ls, ncol(cart_x))
    for (i in 1:ls) {
      grpi <- cart_x[ct == i, ]
      # Convert to data frame for the dist_to_centroids function
      if (is.vector(grpi)) {
        landmarks_c[i, ] <- cart_x[which(ct == i), ]
      } else {
        landmarks_c[i, ] <- polysphere_center(grpi, dims)
      }
    }
    return(rbind(landmarks,landmarks_c))
    }
# Find subset of polysphere points to run interpolation on
get_center <-
  function(x,
           angular = TRUE,
           dims = rep(1, ncol(x)),
           landmark_size = .1,
           method = "ward.D2") {
    n <- nrow(x)
    d <- ncol(x)
    d_sphere <- sum(dims)

    # Compute polysphere distances [s^d_1, s^d_2...]
    ind_dj <- comp_ind_dj(dims)

    # Convert from hyperspherical to cartesian coordinates if needed
    cart_x <- x
    if (angular) {
      cart_x <- polyAnglesToSphere(x, dims)
    }
    # Go straight to upper/lower part
    dist <- dist_polysph_matrix(cart_x, ind_dj = ind_dj, std = FALSE)
    # Make sure count of landmarks is integral
    ls <- floor(n * landmark_size)
    # Divide data to clusters
    hc <- hclust(as.dist(dist), method = method)
    ct <- cutree(hc, k = ls)
    # Pick the point closest to the center of each cluster as a landmark
    landmarks <- numeric(ls)
    for (i in 1:ls) {
      grpi <- x[ct == i, ]
      # Convert to data frame for the dist_to_centroids function
      if (is.vector(grpi)) {
        landmarks[i] <- which(ct == i)
      } else {
        dgrpi <- data.frame(grpi, Group = rep(1, nrow(grpi)))
        grp_dists <- dist_to_centroids(dist_subset(dist, which(ct == i)), dgrpi$Group)
        landmarks[i] <- strtoi(grp_dists$Item[which.min(grp_dists$CentroidDistance)])
      }
    }
    return(landmarks)
  }
# Find centroids of polysphere points to run interpolation on
get_centroids <-
  function(x,
           angular = TRUE,
           dims = rep(1, ncol(x)),
           landmark_size = .1,
           method = "ward.D2") {
    n <- nrow(x)
    d <- ncol(x)
    d_sphere <- sum(dims)

    # Compute polysphere distances [s^d_1, s^d_2...]
    ind_dj <- comp_ind_dj(dims)

    # Convert from hyperspherical to cartesian coordinates if needed
    cart_x <- x
    if (angular) {
      cart_x <- polyAnglesToSphere(x, dims)
    }
    # Go straight to upper/lower part
    dist <- dist_polysph_matrix(cart_x, ind_dj = ind_dj, std = FALSE)
    # Make sure count of landmarks is integral
    ls <- floor(n * landmark_size)
    # Divide data to clusters
    hc <- hclust(as.dist(dist), method = method)
    ct <- cutree(hc, k = ls)
    # Pick the point closest to the center of each cluster as a landmark
    landmarks <- matrix(0, ls, ncol(cart_x))
    for (i in 1:ls) {
      grpi <- cart_x[ct == i, ]
      # Convert to data frame for the dist_to_centroids function
      if (is.vector(grpi)) {
        landmarks[i, ] <- cart_x[which(ct == i), ]
      } else {
        landmarks[i, ] <- polysphere_center(grpi, dims)
      }
    }
    return(landmarks)
  }
# # Find subset of polysphere points to run interpolation on.
# get_centers <-
#   function(x, dist, landmark_size, method = "ward.D2") {
#     n <- nrow(x)
#     d <- ncol(x)
#     d_sphere <- sum(dims)
#     # Initial guess for radius is number of spheres in polysphere
#     r_star <- sqrt(length(dims)) # Compare sum of squares
# 
#     # Compute polysphere distances [s^d_1, s^d_2...]
#     ind_dj <- comp_ind_dj(dims)
# 
#     # Convert from hyperspherical to cartesian coordinates if needed.
#     cart_x <- x
#     if (angular) {
#       cart_x <- polyAnglesToSphere(x, dims)
#     }
#     # Go straight to upper/lower part
#     dist <- dist_polysph_matrix(cart_x, ind_dj = ind_dj, std = FALSE)
#     # Make sure count of landmarks is integral
#     ls <- floor(n * landmark_size)
#     # Divide data to clusters
#     hc <- hclust(as.dist(dist), method = method)
#     ct <- cutree(hc, k = ls)
#     # Pick the point closest to the center of each cluster as a landmark
#     landmarks <- numeric(ls)
#     grpi <- x[ct == i, ]
# 
#     # Convert to data frame for the dist_to_centroids function
#     if (is.vector(grpi)) {
#       landmarks[i] <- which(ct == i)
#     } else {
#       dgrpi <- data.frame(grpi, Group = rep(1, nrow(grpi)))
#       grp_dists <- dist_to_centroids(dist[which(ct == i), which(ct == i)], dgrpi$Group)
# 
#       landmarks[i] <- strtoi(grp_dists$Item[which.min(grp_dists$CentroidDistance)])
#     }
# 
#     return(landmarks)
#   }

# get Frechet mean on polysphere
polysphere_center <- function(x, dims,tol=5e-4) {
  n <- nrow(x)
  d <- ncol(x)
  d_sphere <- sum(dims)

  # Compute polysphere distances [s^d_1, s^d_2...]
  ind_dj <- comp_ind_dj(dims)
  x <- matrix(x, nrow = n)
  x_initial <- x[1, ]
  polydist_initial <- polydist(x_initial, x, ind_dj)
  for (i in 1:n) {
    polydist_new <- polydist(x[i, ], x, ind_dj)
    if (polydist_initial > polydist_new) {
      polydist_initial <- polydist_new
      x_initial <- x[i, ]
    }
  }
  res <- optim(
    fn = polydist,
    par = c(t(matrix(x[1, ]))), group = x, ind_dj = ind_dj,
    method = "BFGS",
    control = list(reltol = tol)
  )
  ans <- res$par
  return(normalize_poly(res$par, dims))
}
# returns sum of squared distances from points on polysphere
polydist <- function(x1, group, ind_dj) {
  group_n <- nrow(group)
  x1 <- x1 / sqrt(sum(x1^2))
  squared_dists <- polykde::dist_polysph(matrix(group, nrow = group_n), t(matrix(x1)), norm_x = TRUE, norm_y = TRUE, ind_dj)**2
  return(sum(squared_dists))
}
# normalize vector to unit polysphere
normalize_poly <- function(x, dims) {
  results <- rep(0, length(dims) + sum(dims))
  start <- 1
  for (i in 1:length(dims)) {
    end <- start + dims[i]
    temp <- x[start:end]
    results[start:end] <- temp / sqrt(sum(temp^2))
    start <- end + 1
  }

  return(results)
}
