### Load Libraries & Setup
# install.packages(c("rotasym", "sphunif", "shapes", "smacof" , "sdetorus", "rgl", "BiocManager"))
# BiocManager::install("waddR")

### Setup
library(rotasym)
library(polykde)
library(waddR)
library(sphunif)
library(shapes)
library(smacof)
library(sdetorus)
library(rgl)
library(plot3D)
library(usethis)

## MDS & PNS for polysphere data (steps 1 and 2 of PCA)
mds_pns <-
  function(x,
           init = "cmds",
           angular = TRUE,
           dims = rep(1, ncol(x)),
           col = NULL,
           method = "Main",
           tol = 5e-3,
           landmark_tol = 5e-4,
           pc = min(3, ncol(x)),
           landmark_size = .1,
           gradient = TRUE) {
    
    # Dimensions of dataset
    n <- nrow(x)
    d_sphere <- sum(dims)

    ls <- n * landmark_size

    if (is.null(col)) {
      
      col <- rep("black", n)
      
    }

    # Initial guess for radius is number of spheres in polysphere
    r_star <- sqrt(length(dims)) # Compare sum of squares

    # Compute polysphere indices [s^d_1, s^d_2...]
    ind_dj <- comp_ind_dj(dims)

    # Convert from hyperspherical to cartesian coordinates if needed
    cart_x <- x
    if (angular) {
      
      cart_x <- angles_to_polysph(x, dims)
      
    }
    
    # Go straight to upper/lower part
    dist <- dist_polysph_matrix(cart_x, ind_dj = ind_dj, std = FALSE)
    d_ij_polysphere <- dist[upper.tri(dist)]
    r_star_initial <- r_star
    stress <- function(y, r_star, n, d_sphere, d_ij_polysphere) {
      r_star <- y[n * (d_sphere + 1) + 1] # r_star in optimization
      y <-
        matrix(y[1:(n * (d_sphere + 1))], nrow = n, ncol = d_sphere + 1)
      y <- y / sqrt(rowSums(y^2)) # Normalization of y

      dim(y) <- c(dim(y), 1) # As an array
      
      # Spherical distance
      ind_tri <- sphunif::upper_tri_ind(n)
      d_ij_sphere <- sphunif::Psi_mat(
        data = y,
        use_ind_tri = TRUE,
        ind_tri = ind_tri
      )

      d_ij_sphere <- drop(d_ij_sphere)
      
      # Stress function (only sum the distinct distances and disregard diagonal)
      return(sum((d_ij_polysphere - r_star * d_ij_sphere)^2) / length(d_ij_polysphere))
    }
    
    # Take landmarks as the centroids of clusters
    landmarks <- get_centroids(x = x, dims = dims, angular = angular, landmark_size = landmark_size)
    
    # Identify landmarks that are contained in the dataset to avoid redundancy
    not_landmarks <- unique_indices(cart_x, landmarks)
    landmark_subset <- subset_indices(cart_x, landmarks)
    landmark_subset_direct <- subset_indices(landmarks, cart_x)

    # Project the landmarks onto the sphere
    dist_ls <- dist_polysph_matrix(landmarks, ind_dj = ind_dj, std = FALSE)
    d_ij_polysphere_ls <- dist_ls[upper.tri(dist_ls)]
    
    # Initialization by cmdscale
    y_ls <-
      cmdscale(
        d = dist_ls,
        k = d_sphere + 1,
        eig = FALSE,
        x.ret = FALSE
      )
    y <- matrix(0, nrow(cart_x), d_sphere + 1)
    
    # Define stress function on landmarks
    landmark_stress <- function(y, d_sphere, r_star, d_ij_polysphere_ls) {
      
      r_star <- y[ls * (d_sphere + 1) + 1] # r_star in optimization
      y <-
        matrix(y[1:(ls * (d_sphere + 1))], nrow = ls, ncol = d_sphere + 1)
      y <- y / sqrt(rowSums(y^2)) # Normalization of y

      dim(y) <- c(dim(y), 1) # As an array
     
      # Spherical distance
      ind_tri <- sphunif::upper_tri_ind(ls)
      d_ij_sphere <- sphunif::Psi_mat(
        data = y,
        use_ind_tri = TRUE,
        ind_tri = ind_tri
      )
      d_ij_sphere <- drop(d_ij_sphere)
      
      # Stress function (only sum the distinct distances and disregard diagonal)
      return(sum((d_ij_polysphere_ls - r_star * d_ij_sphere)^2) / length(d_ij_polysphere_ls))
    
    }

    # Solves MDS for a single point directly given landmark centroids ls
    projection_stress <- function(yp, d_sphere, dims, n, yl, r_star, dist) {
      
      # Landmark points defined outside the function
      y1 <-
        matrix(yp, ncol = d_sphere + 1)
      y1 <- y1 / sqrt(rowSums(y1^2)) # Normalization of y
      y2 <- y1
      dim(y2) <- c(dim(y2), 1) # As an array

      ind_dj <- comp_ind_dj(d_sphere)
      sph_mat <- dist_polysph_cross(x = y1, y = yl, ind_dj = ind_dj, norm_x = TRUE, norm_y = TRUE, std = FALSE)

      # Stress function (only sum the distinct distances and disregard diagonal)
      out <- sum((dist - r_star * sph_mat)^2) / length(dist)
      attr(out, "gradient") <- grad_mds(dist = dist, rho = r_star, x = c(y1), landmarks = yl)
      return(out)
      
    }
    
    # Solve landmark MDS first
    res_landmark <- optim(
      fn = landmark_stress,
      par = c(y_ls, r_star), d_sphere = d_sphere, d_ij_polysphere_ls = d_ij_polysphere_ls,
      method = "BFGS",
      control = list(reltol = landmark_tol)
    )
    r_star <- res_landmark$par[ls * (d_sphere + 1) + 1] # r_star in optimization
    yl <-
      matrix(res_landmark$par[1:(ls * (d_sphere + 1))], nrow = ls, ncol = d_sphere + 1)
    
    # Normalization of y
    yl <- yl / sqrt(rowSums(yl^2))
    y[landmark_subset, ] <- yl[landmark_subset_direct, ]

    ind_dj_poly <- comp_ind_dj(dims)
    if (length(not_landmarks) > 0) {

      # Distances between original points and landmarks on the polysphere
      dist_initial <- dist_polysph_cross(x = cart_x[not_landmarks, ], y = landmarks, norm_x = TRUE, norm_y = TRUE, ind_dj = ind_dj_poly, std = FALSE)

      # Initial guess for non landmark points is weighted average of landmarks
      yp <- weighted_sphere(rows = length(not_landmarks), landmarks = yl, dist_initial = dist_initial, d_sphere = d_sphere)
      y_initial <- y
      y_initial[not_landmarks, ] <- yp
      res <- matrix(0, nrow = nrow(yp), ncol = d_sphere + 1)
      for (i in 1:length(not_landmarks)) {
        
        # Compare landmarks to non-landmarks
        if (gradient) {
          
          resi <- nlm(projection_stress, p = c(yp[i, ]), r_star = r_star, yl = yl, d_sphere = d_sphere, n = n, dist = dist_initial[i, ])
          res[i, ] <- resi$estimate
          
        } else {
          
          resi <- optim(
            fn = projection_stress,
            par = c(yp[i, ]), r_star = r_star, yl = yl, d_sphere = d_sphere, n = n, dist = dist_initial[i, ],
            method = "BFGS",
            control = list(reltol = tol)
          )
          res[i, ] <- resi$par
        
        }
        
      }
      
      y[not_landmarks, ] <- res
      
    } else {
      
      y <- yl
      
    }
    
    y <- y / sqrt(rowSums(y^2))
    
    # Visualize output
    stress_initial <- stress(
      y = c(y_initial, r_star_initial),
      r_star = r_star_initial, n = n, d_sphere = d_sphere, d_ij_polysphere = d_ij_polysphere
    )
    fun <- stress(c(y, r_star), r_star, n = n, d_sphere = d_sphere, d_ij_polysphere = d_ij_polysphere)
    mds_info <- list(
      stress_init = stress_initial,
      stress_final = fun,
      rstar = r_star_initial,
      r.hat = r_star
    )
    out <-
      matrix(y, nrow = n, ncol = d_sphere + 1)
    out <- out / sqrt(rowSums(out^2))
    general <-
      pns(t(out), sphere.type = "seq.test", R = 1) # PNS with sequential testing for type of fitting sphere
    
    # Plotting
    if (ncol(x) < pc && d_sphere > 2) {
      
      pairs(
        t(general$resmat[1:ncol, ]),
        col = col,
        xlim = c(-pi, pi), ylim = c(-pi, pi),
        main = paste("First", ncol, "components, r* =", r_star, collapse = " ")
      )
      
    } else if (d_sphere > 2) {
      
      pairs(
        t(general$resmat[1:pc, ]),
        col = col,
        xlim = c(-pi, pi), ylim = c(-pi, pi),
        main = paste("First", pc, "components, r* =", r_star, collapse = " ")
      )
      
    }

    return(
      list(
        mds_data = out,
        mds_info = mds_info,
        pns_data = t(general$resmat),
        pns_info = general$PNS,
        pns_var = general$percent,
        radius = r_star
      )
    )
    
  }

## Inverting principal curve onto polysphere by MDS
inverse_mds_opt <-
  function(x,
           mds,
           init = "cmds",
           angular = TRUE,
           dims = rep(1, ncol(x)),
           col = NULL,
           method = "Main",
           tol = 5e-3,
           landmark_tol = 5e-4,
           pc = 1,
           grid_size = 100,
           plot_pc = TRUE,
           gradient = TRUE) {
    
    # Convert to angular coordinates
    if (angular) {
      
      x <- angles_to_polysph(x, dims)
      
    }
    
    n <- nrow(x)
    d_sphere <- sum(dims)
    comp_x <- mds$mds_data
    rect_x <- angles_to_polysph(x, dims)

    # Evenly spaced points along first principal curve
    sc <- march_pns1_shapes(seq(0, mds$pns_info$radii[d_sphere] * 2 * pi, l = grid_size),
      PNS = mds$pns_info
    )
    dist_initial <- dist_polysph_cross(sc, comp_x, ind_dj = comp_ind_dj(d_sphere), std = FALSE)

    range <- seq(1, n, floor(n / n))
    r_star <- mds$mds_info$rstar
    poly_pc_initial <- weighted_polysphere(rows = nrow(sc), landmarks = rect_x[range, ], dist_initial = dist_initial[, range], dims = dims)
    res <- matrix(0, nrow = nrow(poly_pc_initial), ncol = sum(dims) + length(dims))

    for (i in 1:grid_size) {
      
      resi <- optim(
        fn = poly_projection_stress,
        par = c(poly_pc_initial[i, ]), r_star = r_star, yl = rect_x[range, ], dims = dims, d_sphere = d_sphere, n = n, dist = dist_initial[i, range],
        method = "BFGS",
        control = list(reltol = tol)
      )
      res[i, ] <- normalize_poly(resi$par, dims)
      
    }
    
    poly_pc <- res
    return(poly_pc)
    
  }

# Inverts PNS onto the polysphere using either Nadaraya-Watson or inverse MDS
inverse_mds_nw <-function(x, result, dims = rep(1, ncol(x)), col = rainbow(nrow(x)), grid_size = 100, method = "NW", angular = TRUE, plot = TRUE) {
  
  d_sphere <- sum(dims)
  sc1 <- march_pns1_shapes(seq(0, result$pns_info$radii[d_sphere] * 2 * pi, l = grid_size),
                           PNS = result$pns_info
  )
  opt_h <- opt_nw_sph_polysph(x = sc1, X = comp_x, Y = x, radius = radius, dims = dims, K = dnorm)
  dist_x <- dist_polysph_cross(x=sc1, comp_x, ind_dj = comp_ind_dj(sum(dims)))
  nw_pc <- nw_sph_polysph(x = sc1, X = comp_x, Y = x, radius=radius, dims = dims, dist_x = dist_x, h = opt_h, K = dnorm)
  return(nw_pc)
  
}

# Inverts PNS onto the polysphere using either Nadaraya-Watson or inverse MDS, takes in polysphere with cartesian coordinates
inverse_mds <- function(x, result, dims = rep(1, ncol(x)), col = rainbow(nrow(x)), grid_size = 100, method = "NW", angular = TRUE, plot = TRUE) {
  
  d_sphere <- sum(dims)
  sc1 <- march_pns1_shapes(seq(0, result$pns_info$radii[d_sphere] * 2 * pi, l = grid_size),
    PNS = result$pns_info
  )
  comp_x <- result$mds_data
  radius <- result$radius
  nw_pc <- NULL
  poly_pc <- NULL
  if (method == "NW") {
    
    # Calculate optimum bandwidth for Nadaraya-Watson
    opt_h <- opt_nw_sph_polysph(x = sc1, X = comp_x, Y = x, radius = radius, dims = dims, K = dnorm)
    
    # Distance matrix on the sphere for the principal curve and MDS mapping
    dist_x <- dist_polysph_cross(x = sc1, comp_x, ind_dj = comp_ind_dj(sum(dims)))
    
    nw_pc <- nw_sph_polysph(x = sc1, X = comp_x, Y = x, radius=radius, dims = dims, dist_x = dist_x, h = opt_h, K = dnorm)
    out <- nw_pc
    
  } else if (method == "both") {
    
    opt_h <- opt_nw_sph_polysph(x = sc1, X = comp_x, Y = x, radius = radius, dims = dims, K = dnorm)
    dist_x <- dist_polysph_cross(x = sc1, comp_x, ind_dj = comp_ind_dj(sum(dims)))
    nw_pc <- nw_sph_polysph(x = sc1, X = comp_x, Y = x, radius=radius, dims = dims, dist_x = dist_x, h = opt_h, K = dnorm)
    out <- nw_pc
    poly_pc_angle <- inverse_mds_opt(x = x, angular = FALSE, dims = dims, grid_size = grid_size, mds = result, plot_pc = FALSE)
    poly_pc <- angles_to_polysph(poly_pc_angle, dims)
    
  } else {
    
    poly_pc_angle <- inverse_mds_opt(x = x, angular = FALSE, dims = dims, grid_size = grid_size, mds = result, plot_pc = FALSE)
    poly_pc <- angles_to_polysph(poly_pc_angle, dims)
    out <- poly_pc
  }
  
  if (plot) {
    
    plot_inverse(x = x, poly_pc = out, dims = dims, col = col, grid_size = grid_size, angular = angular, method = method)
  
  }
  
  return(out)
  
}

# Inverts PNS onto the polysphere via Nadaraya-Watson with a range of bandwidths chosen by scalar multiples of the optimum
inverse_mds_bw <- function(x, result, dims = rep(1, ncol(x)), col = rainbow(nrow(x)), bws, grid_size = 20, method = "NW", angular = TRUE, plot = TRUE) {
  
  d_sphere <- sum(dims)
  sc1 <- march_pns1_shapes(seq(0, result$pns_info$radii[d_sphere] * 2 * pi, l = grid_size),
                           PNS = result$pns_info
  )
  comp_x <- result$mds_data
  radius <- result$radius
  out <- matrix(0, nrow = grid_size*length(bws), ncol = sum(dims)+length(dims))
  
  # Calculate optimum bandwidth for Nadaraya-Watson
  opt_h <- opt_nw_sph_polysph(x = sc1, X = comp_x, Y = x, radius = radius, dims = dims, K = dnorm)
  
  # Distance matrix on the sphere for the principal curve and MDS mapping
  dist_x <- dist_polysph_cross(x = sc1, comp_x, ind_dj = comp_ind_dj(sum(dims)))
  
  for(i in 1:length(bws)){
    
    nw_pc <- nw_sph_polysph(x = sc1, X = comp_x, Y = x, radius=radius, dims = dims, dist_x = dist_x, h = opt_h*bws[i], K = dnorm)
    out[((i-1)*grid_size+1):(i*grid_size),] <- nw_pc
    
  } 
  return(out)
}

# Returns projection of data to sphere based on weighted distance to landmarks
weighted_sphere <- function(rows, landmarks, dist_initial, d_sphere) {
  
  results <- matrix(0, rows, d_sphere + 1)
  
  for (i in 1:rows) {
    
    # Get normalized weights from distance matrix
    weights <- 1 / dist_initial[i, ]
    weights <- weights / sum(weights)
    res <- weights %*% landmarks

    weighted_mean <- res
    results[i, ] <- weighted_mean
    
  }
  
  results <- results / sqrt(rowSums(results^2))
  return(results)
  
}

# Returns weighted sum of squared distances from landmarks on sphere
sphdist <- function(x1, group, weights, ind_dj) {
  
  group_n <- nrow(group)
  x1 <- x1 / sqrt(sum(x1^2))
  squared_dists <- weights %*% polykde::dist_polysph(matrix(group, nrow = group_n), t(matrix(x1)), norm_x = TRUE, norm_y = TRUE, ind_dj = ind_dj)**2
  return(sum(squared_dists))
  
}

# Returns gradient for interpolation of a single point from landmarks
grad_mds <- function(dist, rho, x, landmarks) {
  
  d <- length(x)
  len_ls <- nrow(landmarks)
  grad_total <- rep(0, d)
  # Gradient equation
  for (j in 1:len_ls) {
    
    y <- landmarks[j, ]
    # Normalized inner product between point and landmark
    normed_inner <- (x %*% y) / (norm_vec(x) * norm_vec(y))
    entry <- 2 * rho * (dist[j] - rho * acos(normed_inner)) / sqrt(1 - normed_inner^2)
    entry <- entry[1] * ((diag(d) - outer(x / norm_vec(x), x / norm_vec(x))) %*% y) / (norm_vec(x) * norm_vec(y))
    grad_total <- grad_total + entry
    
  }
  
  out <- c(grad_total / len_ls)
  return(out)
  
}

# Returns projection of data to polysphere based on weighted distance to landmarks
weighted_polysphere <- function(rows, landmarks, dist_initial, dims) {
  
  results <- matrix(0, rows, sum(dims) + length(dims))
  for (i in 1:rows) {
    
    # Get normalized weights from distance matrix
    weights <- 1 / dist_initial[i, ]
    weights <- weights / sum(weights)
    res <- weights %*% landmarks
    
    results[i, ] <- res
  }
  
  results <- polykde::proj_polysph(results, comp_ind_dj(dims))
  
  return(results)
  
}

# Solves MDS from the sphere to the polysphere for a single point directly given landmark centroids ls
poly_projection_stress <- function(yp, d_sphere, dims, n, yl, r_star, dist) {
  
  # landmark points defined outside the function
  y1 <- matrix(normalize_poly(yp, dims), ncol = length(dims) + sum(dims))

  ind_dj <- comp_ind_dj(dims)
  polysph_mat <- dist_polysph_cross(x = y1, y = yl, ind_dj = ind_dj, norm_x = TRUE, norm_y = TRUE, std = FALSE)

  # Stress function (only sum the distinct distances and disregard diagonal)
  out <- sum((r_star * dist - polysph_mat)^2) / length(dist)
  attr(out, "gradient") <- grad_mds(dist = dist, rho = r_star, x = c(y1), landmarks = yl)
  return(out)
  
}
