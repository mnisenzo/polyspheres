library(rgl)
library(shapes)
library(rotasym)
library(polykde)
library(waddR)
library(plot3D)
library(usethis)


# run MDS on polysphere, plot the first three principal components, and plot the residuals in the first 3 dimensions
plot_MDS_py <- function(x, dims = rep(1, ncol(x)), landmark_size = .1, angular = TRUE, pc = min(3, ncol(x)), col = rainbow(nrow(x)), tol = 5e-3, landmark_tol = 5e-4, gradient = FALSE) {
  d_sphere <- sum(dims)
  simulation_2d <- MDS_PNS_py(x = x, angular = angular, dims = dims, col = col, landmark_size = landmark_size, tol = tol, landmark_tol = landmark_tol, method = "Main", pc = pc, gradient = gradient)
  print(simulation_2d$MDS_info)
  sc <- march_PNS1_shapes(seq(0, simulation_2d$pns_info$radii[d_sphere] * 2 * pi, l = 100),
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
    # compare first two principal components
    plot(rbind(simulation_2d$pns_data[, c(1, 2)], cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),
      main = paste("First two Components"), col = c(col, rep("black", 100)), xlim = c(-pi, pi), ylim = c(-pi / 2, pi / 2), axes = FALSE,
      xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
    )
    points(lin_interp(cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)), N = 10))
    sdetorus::torusAxis(sides = 1)
    axis(side = 2, at = c(-pi, -pi / 2, 0, pi / 2, pi), labels = expression(-pi / 2, -pi / 4, 0, pi / 4, pi / 2))

    sc <- march_PNS1_shapes(seq(0, simulation_2d$pns_info$radii[d_sphere] * 2 * pi, l = 100),
      PNS = simulation_2d$pns_info
    )
    plot3d(rbind(simulation_2d$MDS_data, sc)[, 1:3],
      col = c(col, rep("black", 100)),
      xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), box = FALSE, axes = FALSE,
      xlab = "", ylab = "", zlab = "", pch = 19
    )
    sphere1.f(x0 = 0, y0 = 0, z0 = 0, r = 1, col = "azure2", lit = FALSE, alpha = 0.5)
  }
  return(simulation_2d)
}

plot_inverse <- function(x, nw_pc = NULL, poly_pc = NULL, col = rainbow(nrow(x)), dims = rep(1, ncol(x)), grid_size = 100, angular = TRUE, method = "NW") {
  if (angular) {
    clusters <- polySphereToAngles(x, dims)
    if (!is.null(nw_pc)) {
      nw_pc <- polySphereToAngles(nw_pc, dims)
    }
    if (!is.null(poly_pc)) {
      poly_pc <- polySphereToAngles(poly_pc, dims)
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
    out <- rbind(clusters, nw_pc)
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
  } else if (method == "both") {
    out <- rbind(clusters, nw_pc, poly_pc)
    if (angular) {
      out <- out - pi
    }
    if (ncol(out) < 3) {
      plot(out,
        col = c(col, rep("black", grid_size), rep("grey", grid_size)),
        xlab = "", ylab = "", xlim = xlim, ylim = ylim, pch = 19, xaxs = "i", yaxs = "i"
      )
    } else {
      plot3d(out, xlim = xlim, ylim = ylim, zlim = zlim, col = c(col, rep("black", grid_size), rep("grey", grid_size)))
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
run_MDS_py <- function(x, dims = rep(1, ncol(x)), landmark_size = .1, angular = TRUE, pc = min(3, ncol(x)), tol = 5e-3, landmark_tol = 5e-4, plot = FALSE, gradient = FALSE) {
  n <- nrow(x)
  d_sphere <- sum(dims)
  simulation_2d <- MDS_PNS_py(x = x, angular = angular, dims = dims, col = rainbow(n), landmark_size = landmark_size, tol = tol, landmark_tol = landmark_tol, method = "Main", pc = pc, gradient = gradient)
  if (plot) {
    print(simulation_2d$MDS_info)
    sc <- march_PNS1_shapes(seq(0, simulation_2d$pns_info$radii[d_sphere] * 2 * pi, l = 100),
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
      # compare first two principal components
      plot(rbind(simulation_2d$pns_data[, c(1, 2)], cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),
        main = paste("First two Components"), col = c(rainbow(n), rep("black", 100)), xlim = c(-pi, pi), ylim = c(-pi / 2, pi / 2), axes = FALSE,
        xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
      )
      points(lin_interp(cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)), N = 10))
      sdetorus::torusAxis(sides = 1)
      axis(side = 2, at = c(-pi, -pi / 2, 0, pi / 2, pi), labels = expression(-pi / 2, -pi / 4, 0, pi / 4, pi / 2))

      sc <- march_PNS1_shapes(seq(0, simulation_2d$pns_info$radii[d_sphere] * 2 * pi, l = 100),
        PNS = simulation_2d$pns_info
      )
      plot3d(rbind(simulation_2d$MDS_data, sc)[, 1:3],
        col = c(rainbow(n), rep("black", 100)),
        xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), box = FALSE, axes = FALSE,
        xlab = "", ylab = "", zlab = "", pch = 19
      )
      sphere1.f(x0 = 0, y0 = 0, z0 = 0, r = 1, col = "azure2", lit = FALSE, alpha = 0.5)
    }
  }
  return(list(pns_data = simulation_2d$pns_data, MDS_info = simulation_2d$MDS_info))
}


# Compare centroid and interpolation MDS accuracy and runtime with default
compare_MDS <- function(x, dims = rep(1, ncol(x)), landmark_size = .1, angular = TRUE, pc = min(3, ncol(x)), tol = 5e-3, landmark_tol = 5e-4, plot = FALSE, default = FALSE) {

  time <- Sys.time()
  res_gradient <- run_MDS_py(x, landmark_size = landmark_size, dims = dims, angular = angular, plot = plot, gradient = TRUE)
  gradient_time <- Sys.time() - time

  time <- Sys.time()
  res_centroid <- run_MDS_py(x, landmark_size = landmark_size, dims = dims, angular = angular, plot = plot, gradient = FALSE)
  centroid_time <- Sys.time() - time

  # Return function values for residual and centroid MDS
  value_gradient <- res_gradient$MDS_info$stress_final
  value_centroids <- res_centroid$MDS_info$stress_final
  function_values <- c(value_gradient, value_centroids)
  # Account for minimum R^2 scores up to a flip
  times <- c(gradient_time, centroid_time)
  if (default) {
    time <- Sys.time()
    res <- run_MDS(x, dims = dims, angular = angular, tol = tol, plot = plot)
    normal_time <- Sys.time() - time
    value_MDS <- res$MDS_info$stress_final
    function_values <- c(value_gradient, value_centroids, value_MDS)
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
    trial_res <- compare_MDS(clusters, landmark_size = landmark_size, dims = dims, angular = angular, landmark_tol = 5e-4, tol = 5e-4, plot = FALSE, default = default)

    res <- res + trial_res
  }
  return(res / trials)
}
# Compare predicted 1D data with a given function using the formula RSS/TSS up to a flip
modified_R2 <- function(y, f) {
  n <- length(y)
  flip_y <- y[c(n:1)]
  result <- min(sum((y - f)^2), sum((flip_y - f)^2)) / (n * var(f))
  return(result)
}
## Transformation from Space of Principal scores to original sphere
march_PNS1_shapes <- function(theta, PNS) {
  N <- length(theta)
  d <- length(PNS$mean) - 1
  t(shapes::PNSe2s(matrix(
    c(theta, rep(0, N * (d - 1))),
    ncol = N,
    nrow = d,
    byrow = TRUE
  ), PNS))
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
