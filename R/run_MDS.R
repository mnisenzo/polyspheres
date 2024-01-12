library(rgl)
library(shapes)
library(rotasym)
library(polykde)
library(waddR)
library(plot3D)
library(usethis)

# Run MDS on polysphere, plot the first three principal components, and plot the residuals in the first 3 dimensions
plot_mds <- function(x, dims = rep(1, ncol(x)), landmark_size = .1, angular = TRUE, pc = min(3, ncol(x)), col = rainbow(nrow(x)), tol = 5e-3, landmark_tol = 5e-4, gradient = FALSE) {
  
  d_sphere <- sum(dims)
  simulation_2d <- mds_pns(x = x, angular = angular, dims = dims, col = col, landmark_size = landmark_size, tol = tol, landmark_tol = landmark_tol, method = "Main", pc = pc, gradient = gradient)
  print(simulation_2d$mds_info)
  sc <- march_pns1_shapes(seq(0, simulation_2d$pns_info$radii[d_sphere] * 2 * pi, l = 100),
    PNS = simulation_2d$pns_info
  )
  if (sum(dims) < 3) {
    
    plot(rbind(simulation_2d$pns_data, cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),
      main = paste("First two Components"),
      col = c(col, rep("black", 100)), xlim = c(-pi, pi), ylim = c(-pi / 2, pi / 2), axes = FALSE,
      xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
    )
    points(lin_interp(cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)), N = 10))
    sdetorus::torusAxis()
    
  } else {
    
    # Compare first two principal components
    plot(rbind(simulation_2d$pns_data[, c(1, 2)], cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),
      main = paste("First two Components"), col = c(col, rep("black", 100)), xlim = c(-pi, pi), ylim = c(-pi / 2, pi / 2), axes = FALSE,
      xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
    )
    points(lin_interp(cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)), N = 10))
    sdetorus::torusAxis(sides = 1)
    axis(side = 2, at = c(-pi, -pi / 2, 0, pi / 2, pi), labels = expression(-pi / 2, -pi / 4, 0, pi / 4, pi / 2))

    sc <- march_pns1_shapes(seq(0, simulation_2d$pns_info$radii[d_sphere] * 2 * pi, l = 100),
      PNS = simulation_2d$pns_info
    )
    plot3d(rbind(simulation_2d$mds_data, sc)[, 1:3],
      col = c(col, rep("black", 100)),
      xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), box = FALSE, axes = FALSE,
      xlab = "", ylab = "", zlab = "", pch = 19
    )
    sphere1.f(x0 = 0, y0 = 0, z0 = 0, r = 1, col = "azure2", lit = FALSE, alpha = 0.5)
  
  }
  
  return(simulation_2d)
  
}

# Plot the inverted principle curve with the original dataset.
plot_inverse <- function(x, poly_pc, col = rainbow(nrow(x)), dims = rep(1, ncol(x)), grid_size = 100, angular = TRUE, method = "NW") {
  
  if (angular) {
    
    clusters <- polysph_to_angles(x, dims)
    if (!is.null(poly_pc)) {
      
      poly_pc <- polysph_to_angles(poly_pc, dims)
      
    }
    
    xlim <- c(-pi, pi)
    ylim <- c(-pi, pi)
    zlim <- c(-pi, pi)
    
  } else {
    
    clusters <- x
    xlim <- c(-1, 1)
    ylim <- c(-1, 1)
    zlim <- c(-1, 1)
    
  }
  
  if (method == "NW") {
    
    out <- rbind(clusters, poly_pc)
    if (angular) {
      
      out <- out - pi
      
    }
    
    if (ncol(out) < 3) {
      
      plot(out,
        col = c(col, rep("black", grid_size)), xlim = xlim, ylim = ylim,
        xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
      )
      
    } else {
      
      plot3d(out, xlim = xlim, ylim = ylim, zlim = zlim, col = c(col, rep("black", grid_size)))
    
    }
    
  } else {
    
    out <- rbind(clusters, poly_pc)
    if (angular) {
      
      out <- out - pi
      
    }
    if (ncol(out) < 3) {
      
      plot(out,
        col = c(col, rep("black", grid_size)), xlim = xlim, ylim = ylim,
        xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
      )
    
    } else {
      
      plot3d(out, xlim = xlim, ylim = ylim, zlim = zlim, col = c(col, rep("black", grid_size)))
    
    }
    
  }
  
  return(out)
  
}

# Run MDS on polysphere, return results
run_mds <- function(x, dims = rep(1, ncol(x)), landmark_size = .1, angular = TRUE, pc = min(3, ncol(x)), tol = 5e-3, landmark_tol = 5e-4, plot = FALSE, gradient = FALSE) {
  
  n <- nrow(x)
  d_sphere <- sum(dims)
  simulation_2d <- mds_pns(x = x, angular = angular, dims = dims, col = rainbow(n), landmark_size = landmark_size, tol = tol, landmark_tol = landmark_tol, method = "Main", pc = pc, gradient = gradient)
  if (plot) {
    
    print(simulation_2d$mds_info)
    sc <- march_pns1_shapes(seq(0, simulation_2d$pns_info$radii[d_sphere] * 2 * pi, l = 100),
      PNS = simulation_2d$pns_info
    )
    if (sum(dims) < 3) {
      
      plot(rbind(simulation_2d$pns_data, cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),
        main = paste("First two Components"),
        col = c(rainbow(n), rep("black", 100)), xlim = c(-pi, pi), ylim = c(-pi / 2, pi / 2), axes = FALSE,
        xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
      )
      points(lin_interp(cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)), N = 10))
      sdetorus::torusAxis()
      
    } else {
      
      # Compare first two principal components
      plot(rbind(simulation_2d$pns_data[, c(1, 2)], cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),
        main = paste("First two Components"), col = c(rainbow(n), rep("black", 100)), xlim = c(-pi, pi), ylim = c(-pi / 2, pi / 2), axes = FALSE,
        xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
      )
      points(lin_interp(cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)), N = 10))
      sdetorus::torusAxis(sides = 1)
      axis(side = 2, at = c(-pi, -pi / 2, 0, pi / 2, pi), labels = expression(-pi / 2, -pi / 4, 0, pi / 4, pi / 2))

      sc <- march_pns1_shapes(seq(0, simulation_2d$pns_info$radii[d_sphere] * 2 * pi, l = 100),
        PNS = simulation_2d$pns_info
      )
      plot3d(rbind(simulation_2d$mds_data, sc)[, 1:3],
        col = c(rainbow(n), rep("black", 100)),
        xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), box = FALSE, axes = FALSE,
        xlab = "", ylab = "", zlab = "", pch = 19
      )
      sphere1.f(x0 = 0, y0 = 0, z0 = 0, r = 1, col = "azure2", lit = FALSE, alpha = 0.5)
  
    }
    
  }
  
  return(list(pns_data = simulation_2d$pns_data, mds_info = simulation_2d$mds_info))

}

# Compare centroid and interpolation MDS accuracy and runtime with default
compare_mds <- function(x, dims = rep(1, ncol(x)), landmark_size = .1, angular = TRUE, pc = min(3, ncol(x)), tol = 5e-3, landmark_tol = 5e-4, plot = FALSE, default = FALSE) {

  time <- Sys.time()
  res_gradient <- run_mds(x, landmark_size = landmark_size, dims = dims, angular = angular, plot = plot, gradient = TRUE)
  gradient_time <- Sys.time() - time

  time <- Sys.time()
  res_centroid <- run_mds(x, landmark_size = landmark_size, dims = dims, angular = angular, plot = plot, gradient = FALSE)
  centroid_time <- Sys.time() - time

  # Return function values for residual and centroid MDS
  value_gradient <- res_gradient$mds_info$stress_final
  value_centroids <- res_centroid$mds_info$stress_final
  function_values <- c(value_gradient, value_centroids)
  
  # Account for minimum R^2 scores up to a flip
  times <- c(gradient_time, centroid_time)
  if (default) {
    
    time <- Sys.time()
    res <- run_mds(x, dims = dims, angular = angular, tol = tol, plot = plot)
    normal_time <- Sys.time() - time
    value_mds <- res$mds_info$stress_final
    function_values <- c(value_gradient, value_centroids, value_mds)
    times <- c(gradient_time, centroid_time, normal_time)
    
  }
  
  return(c(function_values, times))
  
}

# Compare runtimes, function values, and R2 scores of interpolation MDS over multiple trials
compare_trials <- function(n, nspheres, landmark_size, trials, angular = FALSE, dims = rep(1, nspheres), default = FALSE) {
  
  n <- n
  nspheres <- nspheres
  dims <- dims
  angular <- angular
  landmark_size <- landmark_size
  if (default) {
    
    res <- rep(0, 6)
    
  } else {
    
    res <- rep(0, 4)
    
  }
  
  for (i in 1:trials) {
    
    samp_1 <- r_path_s1r(n = n, r = nspheres, k = rep(1, nspheres), angles = TRUE)
    clusters <- samp_1
    trial_res <- compare_mds(clusters, landmark_size = landmark_size, dims = dims, angular = angular, landmark_tol = 5e-4, tol = 5e-4, plot = FALSE, default = default)

    res <- res + trial_res
    
  }
  
  return(res / trials)
  
}

## Transformation from Space of Principal scores to original sphere
march_pns1_shapes <- function(theta, PNS) {
  
  N <- length(theta)
  d <- length(PNS$mean) - 1
  t(shapes::PNSe2s(matrix(
    c(theta, rep(0, N * (d - 1))),
    ncol = N,
    nrow = d,
    byrow = TRUE
  ), PNS))
  
}
