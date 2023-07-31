### Load Libraries & Setup
# install.packages(c("rotasym", "sphunif", "shapes", "smacof" , "sdetorus", "rgl", "BiocManager"))
# BiocManager::install("waddR")

### Setup
library(rotasym)
library(polykde)
library(waddR)
library(sphunif)
library(shapes)
library(usedist)
library(smacof)
library(sdetorus)
library(rgl)
library(MASS)
library(sdetorus)
library(rgl)
library(plot3D)
library(testthat)
library(usethis)
library(reticulate)
source("Script1_Setup_and_Functions.R")
source("clust_centers.R")


scipy_optimize <- reticulate::import("scipy.optimize")
numpy <- reticulate::import("numpy")


### Functions

## MDS & PNS for polysphere data (steps 1 and 2 of ST-PCA)
MDS_PNS_py <-
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
           python = FALSE,
           landmarks_centroid = FALSE) {
    # Dimensions of dataset
    n <- nrow(x)
    d <- ncol(x)
    d_sphere <- sum(dims)

    ls <- n * landmark_size

    if (is.null(col)) {
      col <- rep("black", n)
    }

    # Initial guess for radius is number of spheres in polysphere
    r_star <- sqrt(length(dims)) # Compare sum of squares

    # Compute polysphere indices [s^d_1, s^d_2...]
    ind_dj <- comp_ind_dj(dims)

    # convert from hyperspherical to cartesian coordinates if needed
    cart_x <- x
    if (angular) {
      cart_x <- polyAnglesToSphere(x, dims)
    }
    # go straight to upper/lower part
    dist <- dist_polysph_matrix(cart_x, ind_dj = ind_dj, std = FALSE)
    d_ij_polysphere <- dist[upper.tri(dist)]
    if (init == "cmds") {
      # Initialization by cmdscale
      y <-
        cmdscale(
          d = dist,
          k = d_sphere + 1,
          eig = FALSE,
          x.ret = FALSE
        )
      y <- y / sqrt(rowSums(y^2))
      y_initial <- y
      r_star_initial <- r_star
    }
    # either take landmarks as the centroids of clusters of the nearest point to the centroid
    if (landmarks_centroid) {
      landmarks <- get_centroids(x = x, dims = dims, angular = angular, landmark_size = landmark_size)
      landmarks <- landmarks / sqrt(rowSums(landmarks^2))
      # project the landmarks onto the sphere
      dist_ls <- dist_polysph_matrix(landmarks, ind_dj = ind_dj, std = FALSE)
      d_ij_polysphere_ls <- dist_ls[upper.tri(dist_ls)]
      y_ls <-
        cmdscale(
          d = dist_ls,
          k = d_sphere + 1,
          eig = FALSE,
          x.ret = FALSE
        )
      d_sphere <- ncol(y_ls) - 1
      y_ls <- y_ls / sqrt(rowSums(y_ls^2))
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

      landmark_stress <- function(y, d_sphere, r_star) {
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

      # Solves MDS for y given landmark centroids ls
      remaining_stress <- function(y, d_sphere, dims, n, yl, r_star, dist) {
        # landmark points defined outside the function
        y1 <-
          matrix(y[1:(n * (d_sphere + 1))], nrow = n, ncol = d_sphere + 1)
        y1 <- y1 / sqrt(rowSums(y1^2)) # Normalization of y
        y2 <- y1
        dim(y2) <- c(dim(y2), 1) # As an array
        # Spherical distance
        ind_tri <- sphunif::upper_tri_ind(n)
        d_ij_sphere <- sphunif::Psi_mat( # coerce to dist class
          data = y2,
          use_ind_tri = TRUE,
          ind_tri = ind_tri
        )
        # sph_mat <- matrix(0, n, n)
        # sph_mat[upper.tri(sph_mat, diag = FALSE)] <- d_ij_sphere
        # polysphere-kde
        # sph_mat <- drop(sph_mat)[1:(n - ls), (n - ls + 1):n]
        ind_dj <- comp_ind_dj(d_sphere)
        sph_mat <- dist_polysph_cross(x = y1, y = yl, ind_dj = ind_dj, norm_x = FALSE, norm_y = FALSE, std = FALSE)
        landmark_polysphere <- dist
        # Stress function (only sum the distinct distances and disregard diagonal)
        return(sum((landmark_polysphere - r_star * sph_mat)^2) / length(landmark_polysphere))
      }
      # Solve landmark MDS first
      res_landmark <- optim(
        fn = landmark_stress,
        par = c(y_ls, r_star), d_sphere = d_sphere,
        method = "BFGS",
        control = list(reltol = landmark_tol)
      )
      r_star <- res_landmark$par[ls * (d_sphere + 1) + 1] # r_star in optimization
      yl <-
        matrix(res_landmark$par[1:(ls * (d_sphere + 1))], nrow = ls, ncol = d_sphere + 1)
      yl <- yl / sqrt(rowSums(yl^2)) # Normalization of y

      ind_dj_poly <- comp_ind_dj(dims)
      # distances between original points and landmarks on the polysphere
      dist_initial <- dist_polysph_cross(x = cart_x, y = landmarks, norm_x = TRUE, norm_y = TRUE, ind_dj = ind_dj_poly, std = FALSE)
      y <-
        cmdscale(
          d = dist,
          k = d_sphere + 1,
          eig = FALSE,
          x.ret = FALSE
        )
      y <- y / sqrt(rowSums(y^2))
      y_initial <- y
      # Compare landmarks to non-landmarks
      res <- optim(
        fn = remaining_stress,
        par = c(y), r_star = r_star, yl = yl, d_sphere = d_sphere, n = n, dist = dist_initial,
        method = "BFGS",
        control = list(reltol = tol)
      )
      y <-
        matrix(res$par, nrow = n, ncol = d_sphere + 1)
    } else {
      landmarks <- get_center(x = x, dims = dims, angular = angular, landmark_size = landmark_size)
      not_landmarks <- c(1:nrow(x))
      not_landmarks <- which(!not_landmarks %in% landmarks)
      # Main optimization
      # Objective function
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

      landmark_stress <- function(y, ls, d_sphere, dist, landmarks, r_star) {
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
        dist_sub <- dist[landmarks, landmarks]
        landmark_polysphere <- dist_sub[upper.tri(dist_sub)]
        # Stress function (only sum the distinct distances and disregard diagonal)
        return(sum((landmark_polysphere - r_star * d_ij_sphere)^2) / length(landmark_polysphere))
      }
      # Solves MDS for the non-landmark coordinates yn given landmark coordinates yl
      remaining_stress <- function(yn, yl, d_sphere, n, ls, r_star, landmarks, not_landmarks, dist) {
        # landmark points defined outside the function
        yn <-
          matrix(yn[1:((n - ls) * (d_sphere + 1))], nrow = n - ls, ncol = d_sphere + 1)
        y1 <- rbind(yn, yl)
        y1 <-
          matrix(y1[1:(n * (d_sphere + 1))], nrow = n, ncol = d_sphere + 1)
        y1 <- y1 / sqrt(rowSums(y1^2)) # Normalization of y
        y2 <- y1
        dim(y2) <- c(dim(y2), 1) # As an array
        # Spherical distance
        ind_tri <- sphunif::upper_tri_ind(n)
        d_ij_sphere <- sphunif::Psi_mat( # coerce to dist class
          data = y2,
          use_ind_tri = TRUE,
          ind_tri = ind_tri
        )
        # sph_mat <- matrix(0, n, n)
        # sph_mat[upper.tri(sph_mat, diag = FALSE)] <- d_ij_sphere
        # polysphere-kde
        # sph_mat <- drop(sph_mat)[1:(n - ls), (n - ls + 1):n]
        ind_dj <- comp_ind_dj(d_sphere)
        sph_mat <- dist_polysph_cross(x = y1[1:(n - ls), ], y = y1[(n - ls + 1):n, ], ind_dj = ind_dj, std = FALSE)
        landmark_polysphere <- dist[not_landmarks, landmarks]
        # Stress function (only sum the distinct distances and disregard diagonal)
        return(sum((landmark_polysphere - r_star * sph_mat)^2) / length(landmark_polysphere))
      }
      # Solve landmark MDS first
      res_landmark <- optim(
        fn = landmark_stress,
        par = c(y[landmarks, ], r_star), ls = ls, d_sphere = d_sphere, dist = dist, landmarks = landmarks,
        method = "BFGS",
        control = list(reltol = landmark_tol)
      )
      r_star <- res_landmark$par[ls * (d_sphere + 1) + 1] # r_star in optimization
      yl <-
        matrix(res_landmark$par[1:(ls * (d_sphere + 1))], nrow = ls, ncol = d_sphere + 1)
      y[landmarks, ] <- yl
      # If there are non-landmarks left, run the analysis on them
      if (n - ls > 0) {
        res <- optim(
          fn = remaining_stress,
          par = c(y[not_landmarks, ]), r_star = r_star, yl = yl, d_sphere = d_sphere, n = n, ls = ls, landmarks = landmarks, not_landmarks = not_landmarks, dist = dist,
          method = "BFGS",
          control = list(reltol = tol)
        )
        yn <-
          matrix(res$par, nrow = n - ls, ncol = d_sphere + 1)
        y[not_landmarks, ] <- yn
      }
    }
    # Visualize output
    stress_initial <- stress(
      y = c(y_initial, r_star_initial),
      r = r_star_initial, n = n, d_sphere = d_sphere, d_ij_polysphere = d_ij_polysphere
    )
    fun <- stress(c(y, r_star), r_star, n = n, d_sphere = d_sphere, d_ij_polysphere = d_ij_polysphere)
    # r.hat <- res$x[n * (d_sphere + 1) + 1]
    MDS_info <- list(
      stress_init = stress_initial,
      stress_final = fun,
      rstar = r_star_initial,
      r.hat = r_star
    )
    out <-
      matrix(y, nrow = n, ncol = d_sphere + 1)
    out <- out / sqrt(rowSums(out^2))

    general <-
      pns(t(out), sphere.type = "seq.test") # PNS with sequential testing for type of fitting sphere
    # Plotting
    if (d_sphere == 2) {
      plot3d(
        out,
        col = col,
        xlim = c(-1, 1),
        ylim = c(-1, 1),
        zlim = c(-1, 1)
      )
      plot(
        t(general$resmat),
        col = col,
        main = paste("First two Components, r* =", r_star),
        asp = 1
      )
    } else if (ncol(x) < pc) {
      pairs(
        t(general$resmat[1:ncol, ]),
        col = col,
        xlim = c(-pi, pi), ylim = c(-pi, pi),
        main = paste("First", ncol, "components, r* =", r_star, collapse = " ")
      )
    } else {
      pairs(
        t(general$resmat[1:pc, ]),
        col = col,
        xlim = c(-pi, pi), ylim = c(-pi, pi),
        main = paste("First", pc, "components, r* =", r_star, collapse = " ")
      )
    }

    return(
      list(
        MDS_data = out,
        MDS_info = MDS_info,
        pns_data = t(general$resmat),
        pns_info = general$PNS,
        pns_var = general$percent,
        radius = r_star
      )
    )
  }
