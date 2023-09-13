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
    if(landmarks_centroid == TRUE){
    # Either take landmarks as the centroids of clusters of the nearest point to the centroid
    landmarks <- get_centroids(x = x, dims = dims, angular = angular, landmark_size = landmark_size)
    #identify landmarks that are contained in the dataset to avoid redundancy
    not_landmarks <- unique_indices(cart_x,landmarks)
    landmark_subset <- subset_indices(cart_x,landmarks)
    landmark_subset_direct <- subset_indices(landmarks,cart_x)
    
    # Project the landmarks onto the sphere
    dist_ls <- dist_polysph_matrix(landmarks, ind_dj = ind_dj, std = FALSE)
    d_ij_polysphere_ls <- dist_ls[upper.tri(dist_ls)]
    #
    dist_total <- dist_polysph_matrix(rbind(cart_x,landmarks), ind_dj = ind_dj, std = FALSE)
    d_ij_polysphere_ls2 <- dist_total[upper.tri(dist_total)]
    # project landmarks onto the hypersphere with the original points
    
    y_total <-
      cmdscale(
        d = dist_total,
        k = d_sphere + 1,
        eig = FALSE,
        x.ret = FALSE
      )
    y_total <- y_total / sqrt(rowSums(y_total^2))
    y_initial <-y_total[1:nrow(cart_x),]
    y_ls <- y_total[(nrow(cart_x)+1):(nrow(cart_x)+nrow(landmarks)),]
    
    landmark_stress <- function(y, d_sphere, r_star,d_ij_polysphere_ls) {
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
    remaining_stress <- function(yp, yn, d_sphere, dims, n, yl, r_star, dist) {
      # landmark points defined outside the function
      y1 <-
        matrix(yp, ncol = d_sphere + 1)
      y1 <- y1 / sqrt(rowSums(y1^2)) # Normalization of y
      y2 <- y1
      dim(y2) <- c(dim(y2), 1) # As an array
      
      ind_dj <- comp_ind_dj(d_sphere)
      sph_mat <- dist_polysph_cross(x = y1, y = yl, ind_dj = ind_dj, norm_x = TRUE, norm_y = TRUE, std = FALSE)
      # Stress function (only sum the distinct distances and disregard diagonal)
      return(sum((dist - r_star * sph_mat)^2) / length(dist))
    }
    # Solve landmark MDS first
    res_landmark <- optim(
      fn = landmark_stress,
      par = c(y_ls, r_star), d_sphere = d_sphere,d_ij_polysphere_ls=d_ij_polysphere_ls,
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
    yp <- y_initial[not_landmarks,]
    
    if(length(not_landmarks) > 0){
      dist_initial <- dist_polysph_cross(x = cart_x[not_landmarks,], y = landmarks, norm_x = TRUE, norm_y = TRUE, ind_dj = ind_dj_poly, std = FALSE)
      # # Compare landmarks to non-landmarks
      # res <- optim(
      #   fn = remaining_stress,
      #   par = c(yp),yn=yn, r_star = r_star, yl = yl, d_sphere = d_sphere, n = n, dist = dist_initial,
      #   method = "BFGS",
      #   control = list(reltol = tol)
      # )
      # 
      # y[not_landmarks,] <-
      #   matrix(res$par, nrow = nrow(yp), ncol = d_sphere + 1)
      # y[landmark_subset,] <- yl[landmark_subset_direct,]
      res <- matrix(0, nrow = nrow(yp), ncol = d_sphere + 1)
      for(i in 1:length(not_landmarks)){
        # Compare landmarks to non-landmarks
        resi <- optim(
          fn = remaining_stress,
          par = c(yp[i,]),yn=yn, r_star = r_star, yl = yl, d_sphere = d_sphere, n = n, dist = dist_initial[i,],
          method = "BFGS",
          control = list(reltol = tol)
        )
        res[i,] <- resi$par
      }
      y[not_landmarks,] <- res
      y[landmark_subset,] <- yl[landmark_subset_direct,]
    }else{
      y <- yl
    }
    
    } else {
      landmarks <- get_center(x = x, dims = dims, angular = angular, landmark_size = landmark_size)
      not_landmarks <- c(1:nrow(x))
      not_landmarks <- which(!not_landmarks %in% landmarks)
      # Main optimization
      landmark_stress <- function(y, ls, d_sphere, dist, r_star) {
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
        landmark_polysphere <- dist[upper.tri(dist)]
        # Stress function (only sum the distinct distances and disregard diagonal)
        return(sum((landmark_polysphere - r_star * d_ij_sphere)^2) / length(landmark_polysphere))
      }
      # Solves MDS for the non-landmark coordinates yn given landmark coordinates yl
      remaining_stress <- function(yn, yl, d_sphere, n, r_star, dist) {
        # landmark points defined outside the function
        yn <-
          matrix(yn, ncol = d_sphere + 1)
        y1 <- yn / sqrt(rowSums(yn^2)) # Normalization of y
        y2 <- y1
        dim(y2) <- c(dim(y2), 1) # As an array

        ind_dj <- comp_ind_dj(d_sphere)
        sph_mat <- dist_polysph_cross(x = y1, y = yl, ind_dj = ind_dj, std = FALSE)
        landmark_polysphere <- dist
        # Stress function (only sum the distinct distances and disregard diagonal)
        return(sum((landmark_polysphere - r_star * sph_mat)^2) / length(landmark_polysphere))
      }
      # Solve landmark MDS first
      res_landmark <- optim(
        fn = landmark_stress,
        par = c(y[landmarks, ], r_star), ls = ls, d_sphere = d_sphere, dist = dist[landmarks,landmarks],
        method = "BFGS",
        control = list(reltol = landmark_tol)
      )
      r_star <- res_landmark$par[ls * (d_sphere + 1) + 1] # r_star in optimization
      yl <-
        matrix(res_landmark$par[1:(ls * (d_sphere + 1))], nrow = ls, ncol = d_sphere + 1)
      yl <- yl / sqrt(rowSums(yl^2))
      y[landmarks, ] <- yl
      # If there are non-landmarks left, run the analysis on them
      if (n - ls > 0) {
        # res <- optim(
        #   fn = remaining_stress,
        #   par = c(y[not_landmarks, ]), r_star = r_star, yl = yl, d_sphere = d_sphere, n = n, ls = ls, landmarks = landmarks, not_landmarks = not_landmarks, dist = dist,
        #   method = "BFGS",
        #   control = list(reltol = tol)
        # )
        # yn <-
        #   matrix(res$par, nrow = n - ls, ncol = d_sphere + 1)
        # y[not_landmarks, ] <- yn
        yp <-y_initial[not_landmarks, ]
        res <- matrix(0, nrow = nrow(yp), ncol = d_sphere + 1)
        for(i in 1:(n - ls)){
          # Compare landmarks to non-landmarks
          resi <- optim(
            fn = remaining_stress,
            par = c(yp[i,]), r_star = r_star, yl = yl, d_sphere = d_sphere, n = n, dist = dist[not_landmarks[i],landmarks],
            method = "BFGS",
            control = list(reltol = tol)
          )
          res[i,] <- resi$par
        }
        y[not_landmarks,] <- res
      }
    }
    y <- y / sqrt(rowSums(y^2))
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
    #Plotting
    if (ncol(x) < pc && d_sphere > 2) {
      pairs(
        t(general$resmat[1:ncol, ]),
        col = col,
        xlim = c(-pi, pi), ylim = c(-pi, pi),
        main = paste("First", ncol, "components, r* =", r_star, collapse = " ")
      )
    } else if (d_sphere >2) {
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

# returns indices of rows in first array that are not contained in second array
unique_indices <- function(arr1, arr2){
  inds <- 1:nrow(arr1)
  diff <- intersect(arr1,arr2)
  for (i in inds){
    if(length(intersect(arr1[i,],diff)) == length(arr1[i,])){
      inds[i] = 0
    }
  }
  return(inds[inds != 0])
}

# returns indices of rows in first array that are contained in second array
subset_indices <- function(arr1, arr2){
  inds <- 1:nrow(arr1)
  diff <- intersect(arr1,arr2)
  for (i in inds){
    if(length(intersect(arr1[i,],diff))== 0){
      inds[i] = 0
    }
  }
  return(inds[inds != 0])
}
