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

scipy_optimize <- reticulate::import("scipy.optimize")
numpy <- reticulate::import("numpy")


### Functions

## Angles to Spheres transformation
anglesToSphere <- function(theta) {
  if (!is.matrix(theta)) {
    theta <- matrix(theta, nrow = 1)
  }
  d <- ncol(theta)
  x <- matrix(0, nrow = nrow(theta), ncol = d + 1)
  cosTheta <- cos(theta)
  sinTheta <- sin(theta)
  if (d == 1) {
    x <- cbind(cosTheta, sinTheta)
  } else {
    prodSinTheta <- t(apply(sinTheta, 1, cumprod))
    x <- cbind(
      cosTheta[, 1, drop = FALSE], prodSinTheta[, -d, drop = FALSE] *
        cosTheta[, -1, drop = FALSE],
      prodSinTheta[, d, drop = FALSE]
    )[, , drop = TRUE]
  }
  return(x)
}
## Rectangular to polysphere coordinate transformation
polySphereToAngles <- function(x, dims) {
  if (!is.matrix(theta)) {
    x <- matrix(x, nrow = 1)
  }
  theta <-
    matrix(0, nrow = nrow(x), ncol = sum(dims))
  ind_p <- 0
  ind_r <- 0
  for (d in dims) {
    input <- as.matrix(x[, (ind_p + 1):(ind_p + d+1)])
    if (nrow(x) == 1) {
      input <- t(input)
    }
    theta[, (ind_r + 1):(ind_r + d)] <- sphere_to_angles(input)
    ind_p <- ind_p + d +1
    ind_r <- ind_r + d
  }
  return(theta)
}

# turn degrees into radians
deg2rad <- function(deg) {(deg * pi) / (180)}


## Polysphere to rectangular coordinate transformation
polyAnglesToSphere <- function(theta, dims) {
  if (!is.matrix(theta)) {
    theta <- matrix(theta, nrow = 1)
  }
  x <-
    matrix(0, nrow = nrow(theta), ncol = length(dims) + sum(dims))
  ind_p <- 0
  ind_r <- 0
  for (d in dims) {
    input <- as.matrix(theta[, (ind_p + 1):(ind_p + d)])
    if (nrow(theta) == 1) {
      input <- t(input)
    }
    x[, (ind_r + 1):(ind_r + d + 1)] <- angles_to_sphere(input)
    ind_p <- ind_p + d
    ind_r <- ind_r + d + 1
  }
  return(x)
}

# find subset of polysphere points to run interpolation on.
get_center <-
  function(x,
           angular = TRUE,
           dims = rep(1, ncol(x)),
           landmark_size = .1) {
    n <- nrow(x)
    d <- ncol(x)
    d_sphere <- sum(dims)
    
    if (is.null(col)) {
      col <- rep("black", n)
    }
    
    # initial guess for radius is number of spheres in polysphere
    r_star <- sqrt(length(dims)) # compare sum of squares
    
    # Compute polysphere distances [s^d_1, s^d_2...]
    ind_dj <- comp_ind_dj(dims)
    
    # convert from hyperspherical to cartesian coordinates if needed.
    cart_x <- x
    if (angular) {
      cart_x <- polyAnglesToSphere(x, dims)
    }
    # go straight to upper/lower part
    dist <- dist_polysph_matrix(cart_x, ind_dj = ind_dj, std = FALSE)
    # make sure count of landmarks is integral
    ls <- floor(n*landmark_size)
    # divide data to clusters
    hc <- hclust(as.dist(dist))
    ct <- cutree(hc,k=ls)
    # pick the point closest to the center of each cluster as a landmark
    landmarks <- numeric(ls)
    for(i in 1:ls){
      grpi <- x[ct == i,]
      #convert to data frame for the dist_to_centroids function
      if(is.vector(grpi)){
        landmarks[i] <- which(ct ==i)
      } else{
        dgrpi <- data.frame(grpi, Group = rep(1,nrow(grpi)))
        grp_dists <- dist_to_centroids(dist_subset(dist, which(ct==i)),dgrpi$Group)
        landmarks[i] <- strtoi(grp_dists$Item[which.min(grp_dists$CentroidDistance)])
      }
    }
    return(landmarks)
  }

# find subset of polysphere points to run interpolation on.
get_centers <-
  function(x,dist,landmark_size){
    n <- nrow(x)
    d <- ncol(x)
    d_sphere <- sum(dims)
        # initial guess for radius is number of spheres in polysphere
        r_star <- sqrt(length(dims)) # compare sum of squares
        
        # Compute polysphere distances [s^d_1, s^d_2...]
        ind_dj <- comp_ind_dj(dims)
        
        # convert from hyperspherical to cartesian coordinates if needed.
        cart_x <- x
        if (angular) {
          cart_x <- polyAnglesToSphere(x, dims)
        }
        # go straight to upper/lower part
        dist <- dist_polysph_matrix(cart_x, ind_dj = ind_dj, std = FALSE)
        # make sure count of landmarks is integral
        ls <- floor(n*landmark_size)
        # divide data to clusters
        hc <- hclust(as.dist(dist))
        ct <- cutree(hc,k=ls)
        # pick the point closest to the center of each cluster as a landmark
        landmarks <- numeric(ls)
        grpi <- x[ct == i,]
        
        #convert to data frame for the dist_to_centroids function
        if(is.vector(grpi)){
          landmarks[i] <- which(ct == i)
        } else{
          dgrpi <- data.frame(grpi, Group = rep(1,nrow(grpi)))
          grp_dists <- dist_to_centroids(dist[which(ct == i),which(ct == i)],dgrpi$Group)
          
          landmarks[i] <- strtoi(grp_dists$Item[which.min(grp_dists$CentroidDistance)])
        }

    return(landmarks)
  }

## MDS & PNS for polysphere data (steps 1 and 2 of ST-PCA)
MDS_PNS_py <-
  function(x,
           init = "cmds",
           angular = TRUE,
           dims = rep(1, ncol(x)),
           col = NULL,
           method = "Main",
           tol = 5e-3,
           landmark_tol = 5e-5,
           pc = min(3,ncol(x)),
           landmark_size = .1,
           python = FALSE) {
    # dimensions of dataset
    n <- nrow(x)
    d <- ncol(x)
    d_sphere <- sum(dims)
    
    ls <- n*landmark_size
    
    if (is.null(col)) {
      col <- rep("black", n)
    }
    
    # initial guess for radius is number of spheres in polysphere
    r_star <- sqrt(length(dims)) # compare sum of squares
    
    # Compute polysphere distances [s^d_1, s^d_2...]
    ind_dj <- comp_ind_dj(dims)
    
    # convert from hyperspherical to cartesian coordinates if needed.
    cart_x <- x
    if (angular) {
      cart_x <- polyAnglesToSphere(x, dims)
    }
    # go straight to upper/lower part
    dist <- dist_polysph_matrix(cart_x, ind_dj = ind_dj, std = FALSE)
    d_ij_polysphere <- dist[upper.tri(dist)]
    landmarks <- get_center(x=x,dims=dims, angular = angular, landmark_size=landmark_size)
    not_landmarks <- c(1:nrow(x))
    not_landmarks <- which(!not_landmarks %in% landmarks)
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
    # Main optimization
    if (method == "Main") {
      # Objective function
      stress <- function(y,r_star) {
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
      
      landmark_stress <- function(y,r_star) {
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
        dist_sub <-  dist[landmarks,landmarks]
        landmark_polysphere <- dist_sub[upper.tri(dist_sub)]
        # Stress function (only sum the distinct distances and disregard diagonal)
        return(sum((landmark_polysphere - r_star * d_ij_sphere)^2) / length(landmark_polysphere))
        
      }
      # Solves MDS for the non-landmark coordinates yn given landmark coordinates yl
      remaining_stress <- function(yn) {
        #landmark points defined outside the function
        yn <-
          matrix(yn[1:((n-ls) * (d_sphere + 1))], nrow = n-ls, ncol = d_sphere + 1)
        y <- rbind(yn,yl)
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
        sph_mat <- matrix(0,n,n)
        #sph_mat[lower.tri(sph_mat, diag=FALSE)] <- d_ij_sphere
        sph_mat[upper.tri(sph_mat, diag=FALSE)] <- d_ij_sphere
        sph_mat <- drop(sph_mat)[1:(n-ls),(n-ls+1):n]
        landmark_polysphere <- dist[not_landmarks,landmarks]
        # Stress function (only sum the distinct distances and disregard diagonal)
        return(sum((landmark_polysphere - r_star * sph_mat)^2) / length(landmark_polysphere))
        
      }
      if(python){
        bounds <- list(c(-pi,pi))
        bounds <- append(bounds,rep(list(c(-pi/2,pi/2)),(d_sphere+1)*ls-1))
        bounds <- append(bounds, list(c(0,3*r_star)))
        
        # solve landmark MDS first
       # res_landmark <- scipy_optimize$minimize(landmark_stress,x0 = c(y[landmarks,],r_star), bounds=bounds,tol = tol)
        res_landmark <- scipy_optimize$minimize(landmark_stress,x0 = c(y[landmarks,],r_star), method = "BFGS",tol = tol)
        
        r_star <- res_landmark$x[ls * (d_sphere + 1) + 1] # r_star in optimization
        yl <-
          matrix(res_landmark$x[1:(ls * (d_sphere + 1))], nrow = ls, ncol = d_sphere + 1)
        y[landmarks,] <- yl
        # if there are non-landmarks left, run the analysis on them
        if (n -ls > 0){
          bounds <- rep(list(c(-pi/2,pi/2)),(d_sphere+1)*(n-ls))
        #  res <- scipy_optimize$minimize(remaining_stress,x0 = c(y[not_landmarks,]), bounds=bounds,tol = tol)
          res <- scipy_optimize$minimize(remaining_stress,x0 = c(y[not_landmarks,]), method="BFGS",tol = tol)
          yn <-
            matrix(res$x, nrow = n-ls, ncol = d_sphere + 1)
          y[not_landmarks,] <- yn
        }
      } else{
        # solve landmark MDS first
        res_landmark <- optim(
          fn = landmark_stress,
          par = c(y[landmarks,],r_star),
          method = "BFGS",
          control = list(reltol = landmark_tol)
        )
        r_star <- res_landmark$par[ls * (d_sphere + 1) + 1] # r_star in optimization
        yl <-
          matrix(res_landmark$par[1:(ls * (d_sphere + 1))], nrow = ls, ncol = d_sphere + 1)
        y[landmarks,] <- yl
        # if there are non-landmarks left, run the analysis on them
        if (n -ls > 0){
          res <- optim(
            fn = remaining_stress,
            par = c(y[not_landmarks,]),
            method = "BFGS",
            control = list(reltol = tol)
          )
          yn <-
            matrix(res$par, nrow = n-ls, ncol = d_sphere + 1)
          y[not_landmarks,] <- yn
        }
      }
      
      # Visualize output
      stress_initial <- stress(
        y = c(y_initial,r_star_initial),
        r = r_star_initial
      )
      fun <- stress(c(y,r_star),r_star)
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
    }
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
        xlim=c(-pi,pi), ylim=c(-pi,pi),
        main = paste("First", ncol, "components, r* =", r_star,collapse = " ")
      )
    }
    else {
      pairs(
        t(general$resmat[1:pc, ]),
        col = col,
        xlim=c(-pi,pi), ylim=c(-pi,pi),
        main = paste("First", pc, "components, r* =", r_star,collapse = " ")
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

## Transformation from Space of Principal scores to original sphere
march_PNS1_shapes_py <- function(theta, PNS) {
  N <- length(theta)
  d <- length(PNS$mean) - 1
  t(shapes::PNSe2s(matrix(
    c(theta, rep(0, N * (d - 1))),
    ncol = N,
    nrow = d,
    byrow = TRUE
  ), PNS))
}

## Toroidal distance
torus_distance <- function(a, b) {
  dist.coord <- rep(NA, length(a))
  for (i in seq_along(a)) {
    dist1 <- abs(a[i] - b[i])
    dist2 <- 2 * pi - abs(a[i] - b[i])
    dist.coord[i] <- min(dist1, dist2)^2
  }
  dist.torus <- sqrt(sum(dist.coord))
  return(dist.torus)
}

## compute distance between two points on polysphere
polysphere_distance <-
  function(a,
           b,
           angular = TRUE,
           dims = rep(1, length(a))) {
    if (angular) {
      a <- polyAnglesToSphere(a, dims)
      b <- polyAnglesToSphere(b, dims)
    }
    dist.coord <- rep(NA, length(dims))
    ind <- 0
    for (i in seq_along(dims)) {
      d <- dims[i]
      a1 <- a[(ind + 1):(ind + d + 1)]
      b1 <- b[(ind + 1):(ind + d + 1)]
      ind <- ind + d + 1
      dist.coord[i] <- acos(a1 %*% b1)
    }
    
    dist.sphere <- sqrt(sum(dist.coord**2))
    return(dist.sphere)
  }

## Ancillary functions

sphere1.f <- function(x0 = 0,
                      y0 = 0,
                      z0 = 0,
                      r = 1,
                      n = 101,
                      ...) {
  f <- function(s, t) {
    cbind(
      r * cos(t) * cos(s) + x0,
      r * sin(s) + y0,
      r * sin(t) * cos(s) + z0
    )
  }
  persp3d(
    f,
    slim = c(-pi / 2, pi / 2),
    tlim = c(0, 2 * pi),
    n = n,
    add = TRUE,
    ...
  )
}

lin_interp <- function(X, N = 10) {
  n_X <- nrow(X)
  X_new <- matrix(nrow = N * (n_X - 1) + 1, ncol = ncol(X))
  t0 <- seq(0, 1, l = N + 1)[-(N + 1)]
  t1 <- 1 - t0
  for (i in 1:(n_X - 1)) {
    X_new[N * (i - 1) + 1:N, ] <- t1 %x% X[i, , drop = FALSE] +
      t0 %x% X[i + 1, , drop = FALSE]
  }
  X_new[N * (n_X - 1) + 1, ] <- X[n_X, ]
  return(X_new)
}

# run MDS on polysphere, plot the first three principal components, and plot the residuals in the first 3 dimensions
plot_MDS_py <- function(x, dims=rep(1,ncol(x)), landmark_size=.1, angular=TRUE, pc =min(3,ncol(x)),tol=5e-3, landmark_tol = 5e-5,python=FALSE){
  n = nrow(x)
  d_sphere = sum(dims)
  simulation_2d <- MDS_PNS_py(x = x, angular=angular,dims=dims,col = rainbow(n),landmark_size = landmark_size,tol = tol, landmark_tol = landmark_tol, method = "Main",pc=pc,python=python)
  print(simulation_2d$MDS_info)
  # percent_var_explained(simulation_2d, x = x)
  sc <- march_PNS1_shapes_py(seq(0, simulation_2d$pns_info$radii[d_sphere] * 2 * pi, l = 100),
                          PNS = simulation_2d$pns_info
  )
  
  if (sum(dims) < 3){
    plot(rbind(simulation_2d$pns_data, cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),main = paste("First two Components"),
         col = c(rainbow(n), rep("black", 100)), xlim = c(-pi, pi), ylim = c(-pi/2, pi/2), axes = FALSE,
         xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
    )
    points(lin_interp(cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)), N = 10))
    sdetorus::torusAxis()
  } else{
    # compare first two principal components
    plot(rbind(simulation_2d$pns_data[, c(1,2)], cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),
         main = paste("First two Components"), col = c(rainbow(n), rep("black", 100)), xlim = c(-pi, pi), ylim = c(-pi / 2, pi / 2), axes = FALSE,
         xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
    )
    points(lin_interp(cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)), N = 10))
    sdetorus::torusAxis(sides = 1)
    axis(side = 2, at = c(-pi, -pi / 2, 0, pi / 2, pi), labels = expression(-pi / 2, -pi / 4, 0, pi / 4, pi / 2))
    
    sc <- march_PNS1_shapes(seq(0, simulation_2d$pns_info$radii[d_sphere] * 2 * pi, l = 100),
                            PNS = simulation_2d$pns_info
    )
    plot3d(rbind(simulation_2d$MDS_data, sc)[,1:3],
           col = c(rainbow(n), rep("black", 100)),
           xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), box = FALSE, axes = FALSE,
           xlab = "", ylab = "", zlab = "", pch = 19
    )
    sphere1.f(x0 = 0, y0 = 0, z0 = 0, r = 1, col = "azure2", lit = FALSE, alpha = 0.5)
    
  }
  return(simulation_2d$pns_data)
}

  # run MDS on polysphere, return results
  run_MDS_py <- function(x, dims=rep(1,ncol(x)), landmark_size=.1, angular=TRUE, pc =min(3,ncol(x)),tol=5e-3, landmark_tol = 5e-5,python=FALSE){
    n = nrow(x)
    d_sphere = sum(dims)
    simulation_2d <- MDS_PNS_py(x = x, angular=angular,dims=dims,col = rainbow(n),landmark_size = landmark_size,tol = tol, landmark_tol = landmark_tol, method = "Main",pc=pc,python=python)
    print(simulation_2d$MDS_info$stress_final)
    return(simulation_2d$MDS_info)
  }
  
  
