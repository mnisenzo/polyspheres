## Modification of the prediction algorithm to generate contour plots
## Attention: this function will SAVE all contour frames in the user's directory
## Figure 4c will be saved as a pdf titled simul_2d_contour_56

prediction <- function(x_old, y_old, y_new, r_star, col, contour) {
  all_sphere <- rbind(y_old, y_new)
  n_new <- nrow(y_new) # Assume for now that y_new is a single added point, then n_new = 1
  n_old <- nrow(y_old)
  n <- n_new + n_old
  d <- ncol(y_new)
  x_new_final <- matrix(0, nr = n_new, nc = d - 1)
  
  for (i in c(1:n_new)) {
    y <- rbind(y_old, y_new[i, ])
    
    # Calculate Spherical Distances
    
    # Normalization of y
    y <- y / sqrt(rowSums(y^2))
    
    # As an array
    dim(y) <- c(dim(y), 1)
    
    # Spherical distance
    ind_tri <- sphunif::upper_tri_ind(n_old + 1)
    d_ij_sphere <- sphunif::Psi_mat(
      data = y, use_ind_tri = TRUE,
      ind_tri = ind_tri
    )
    d_ij_sphere <- drop(d_ij_sphere)
    
    dist <- matrix(0, nrow = n_old + 1, ncol = n_old + 1)
    dist[upper.tri(dist)] <- d_ij_sphere # Sets the Upper triangular part of the matrix
    dist <- dist + t(dist) # adds the lower to make it symmetric
    dist3 <- dist[, n_old + 1][1:n_old]
    dist_torus <- matrix(0, nr = n_old + 1, nc = n_old + 1)
    
    # Objective function
    obj <- function(x_new_torus, x_old_torus, n, r, dist) {
      x_new_torus <- (x_new_torus + pi) %% (2 * pi) - pi
      # Calculate torus distances
      x <- rbind(x_old_torus, x_new_torus)
      # Toroidal Distances for full dataset
      ind_tri <- sphunif::upper_tri_ind(n)
      x_array <- array(dim = c(n, 1, d - 1))
      x_array[, 1, ] <- x
      d_ij_circ_each_dim <- sphunif::Psi_mat(
        data = x_array, use_ind_tri = TRUE,
        ind_tri = ind_tri
      )
      d_ij_torus <- sqrt(rowSums(d_ij_circ_each_dim^2))
      
      dist_torus[upper.tri(dist_torus)] <- d_ij_torus # Sets the Upper triangular part of the matrix
      dist_torus <- dist_torus + t(dist_torus)
      dist2 <- dist_torus[, n_old + 1][1:n_old]
      
      # Stress function (only sum the distinct distances and disregard diagonal)
      return(2 * sum((r * dist3 - dist2)^2) / n)
    }
    
    
    # Initialization
    if (i == 1) {
      "      x_new_torus = x_old[which.min(dist3[1:(length(dist2)-1)]),] "
      x_new_torus <- x_old[which.min(dist3), ]
      " x_new_torus = x_old[ceil(rand()*n_old),]"
    } else {
      # x_new_torus = x_new_final[i-1,]
      x_new_torus <- x_old[which.min(dist3), ]
    }
    
    pred <- optim(
      fn = obj, par = c(x_new_torus), x_old_torus = x_old,
      n = n_old + 1, r = r_star, dist = dist3, method = "BFGS",
      control = list(trace = 3)
    )$par
    
    pred <- (pred + pi) %% (2 * pi) - pi
    
    x_new_final[i, ] <- pred
    
    if (contour == 1) {
      k <- 40
      theta <- seq(from = -pi, to = pi, by = (2 * pi / (k - 1)))
      
      B <- matrix(NA, nrow = k, ncol = k)
      for (v in 1:k) {
        for (w in 1:k) {
          B[v, w] <- obj(x_new_torus = c(theta[v], theta[w]), x_old_torus = x_old, n = n_old + 1, r = r_star, dist = dist3)
        }
      }
      
      minimizer <- which(B == min(B), arr.ind = TRUE) # find the minimizer on the grid
      # Save contour plots
      title <- paste("simul_2d_contour", i, sep = "_")
      pdf(title, width = 7, height = 7)
      plot(x_old,
           pch = 19, col = col[1:(n_old + i)], xlab = expression(paste(theta)), ylab = expression(paste(phi)),
           xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE, xaxs = "i", yaxs = "i"
      )
      contour(x = theta, y = theta, z = B, add = TRUE, levels = seq(from = 0, to = 7, by = (10 - 0) / 25), col = rainbow(26))
      points(x_new_final[i, 1], x_new_final[i, 2], pch = 4, col = "darkred") # Plot with red the optimization result
      points(theta[minimizer[1]], theta[minimizer[2]], pch = 8, col = "forestgreen") # Contour minimum
      points(x_new_final[1:i, ], pch = 0, col = rep("black", i))
      sdetorus::torusAxis()
      dev.off()
    }
  }
  return(x_new_final)
}