# Compare centroid and nearest neighbors on the same clustering
compare_centers <-
  function(x,
           angular = TRUE,
           dims = rep(1, ncol(x)),
           landmark_size = .1,
           method = "ward.D2") {
    # Divide data to clusters
    clusters <- get_clusters(x,
                             angular = angular,
                             dims = dims,
                             landmark_size = landmark_size,
                             method = method)
    hc <- clusters$hc
    ct <- clusters$ct
    ls <- clusters$ls
    cart_x <- clusters$cart_x

    # Pick the point closest to the center of each cluster as a landmark
    landmarks <- numeric(ls)
    for (i in 1:ls) {
      grpi <- x[ct == i, ]

      # Convert to data frame for the dist_to_centroids function
      if (is.vector(grpi)) {
        landmarks[i] <- which(ct == i)
        landmarks_c[i, ] <- cart_x[which(ct == i), ]
      } else {
        dgrpi <- data.frame(grpi, Group = rep(1, nrow(grpi)))
        grp_dists <- dist_to_centroids(dist_subset(dist, which(ct == i)), dgrpi$Group)
        landmarks[i] <- strtoi(grp_dists$Item[which.min(grp_dists$CentroidDistance)])
        landmarks_c[i, ] <- polysphere_center(grpi, dims)
      }
    }
    landmarks <- cart_x[landmarks, ]
    landmarks_c <- matrix(0, ls, ncol(cart_x))
    return(rbind(landmarks, landmarks_c))
  }

# Find subset of polysphere points to run interpolation on
get_center <-
  function(x,
           angular = TRUE,
           dims = rep(1, ncol(x)),
           landmark_size = .1,
           method = "ward.D2") {
    # Divide data to clusters
    clusters <- get_clusters(x,
                             angular = angular,
                             dims = dims,
                             landmark_size = landmark_size,
                             method = method)
    hc <- clusters$hc
    ct <- clusters$ct
    ls <- clusters$ls
    cart_x <- clusters$cart_x

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
    
    # Divide data to clusters
    clusters <- get_clusters(x,
                 angular = angular,
                 dims = dims,
                 landmark_size = landmark_size,
                 method = method)
    hc <- clusters$hc
    ct <- clusters$ct
    ls <- clusters$ls
    cart_x <- clusters$cart_x

    # Pick the centroid of each cluster as a landmark
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

# Return the clustering based on a dataset and the number of points
get_clusters <-
  function(x,
           angular = TRUE,
           dims = rep(1, ncol(x)),
           landmark_size = .1,
           method = "ward.D2") {
    n <- nrow(x)
    
    # Make sure count of landmarks is integer
    ls <- floor(n * landmark_size)
    
    # Compute polysphere distances [s^d_1, s^d_2...]
    ind_dj <- comp_ind_dj(dims)
    
    # Convert from hyperspherical to Cartesian coordinates if needed
    cart_x <- x
    if (angular) {
      cart_x <- angles_to_polysph(x, dims)
    }
    
    # Go straight to upper/lower part
    dist <- dist_polysph_matrix(cart_x, ind_dj = ind_dj, std = FALSE)
    
    # Divide data to clusters
    hc <- hclust(as.dist(dist), method = method)
    ct <- cutree(hc, k = ls)
    return(list(hc=hc,ct=ct, ls=ls, cart_x=cart_x))
  }

# Get Frechet mean on polysphere
polysphere_center <- function(x, dims) {
  y <- colSums(x)
  return(normalize_poly(y, dims))
}

# Normalize vector to unit polysphere
normalize_poly <- function(x, dims) {
  if (is.null(nrow(x))) {
    x <- matrix(x, nrow = 1)
  }
  results <- polykde::proj_polysph(x, comp_ind_dj(dims))

  return(results)
}
