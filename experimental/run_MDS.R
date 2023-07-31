source("Script1_Setup_and_Functions.R")
source("efficient_mds.R")
source("STPCA_prediction.R")

# run MDS on polysphere, plot the first three principal components, and plot the residuals in the first 3 dimensions
plot_MDS <- function(x, dims = rep(1, ncol(x)), angular = TRUE, pc = min(3, ncol(x)), tol = 5e-3) {
  n <- nrow(x)
  d_sphere <- sum(dims)
  simulation_2d <- MDS_PNS(x = x, angular = angular, dims = dims, col = rainbow(n), tol = tol, method = "Main", pc = pc)
  print(simulation_2d$MDS_info)
  # percent_var_explained(simulation_2d, x = x)
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
  return(simulation_2d$pns_data)
}
# run MDS on polysphere, return results
run_MDS <- function(x, dims = rep(1, ncol(x)), angular = TRUE, pc = min(3, ncol(x))) {
  n <- nrow(x)
  d_sphere <- sum(dims)
  simulation_2d <- MDS_PNS(x = x, angular = angular, dims = dims, col = rainbow(n), method = "Main", pc = pc)
  return(simulation_2d$pns_data)
}


# run MDS on polysphere, plot the first three principal components, and plot the residuals in the first 3 dimensions
plot_MDS_py <- function(x, dims = rep(1, ncol(x)), landmark_size = .1, landmarks_centroid = FALSE, angular = TRUE, pc = min(3, ncol(x)), tol = 5e-3, landmark_tol = 5e-5, python = FALSE) {
  n <- nrow(x)
  d_sphere <- sum(dims)
  simulation_2d <- MDS_PNS_py(x = x, angular = angular, dims = dims, col = rainbow(n), landmarks_centroid = landmarks_centroid, landmark_size = landmark_size, tol = tol, landmark_tol = landmark_tol, method = "Main", pc = pc, python = python)
  print(simulation_2d$MDS_info)
  # percent_var_explained(simulation_2d, x = x)
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
  return(simulation_2d$pns_data)
}

# run MDS on polysphere, return results
run_MDS_py <- function(x, dims = rep(1, ncol(x)), landmark_size = .1, angular = TRUE, pc = min(3, ncol(x)), tol = 5e-3, landmark_tol = 5e-5, python = FALSE) {
  n <- nrow(x)
  d_sphere <- sum(dims)
  simulation_2d <- MDS_PNS_py(x = x, angular = angular, dims = dims, col = rainbow(n), landmark_size = landmark_size, tol = tol, landmark_tol = landmark_tol, method = "Main", pc = pc, python = python)
  print(simulation_2d$MDS_info$stress_final)
  return(simulation_2d$MDS_info)
}
