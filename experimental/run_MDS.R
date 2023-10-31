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
run_MDS <- function(x, dims = rep(1, ncol(x)), angular = TRUE, pc = min(3, ncol(x)),tol = 5e-3,plot=FALSE) {
  n <- nrow(x)
  d_sphere <- sum(dims)
  simulation_2d <- MDS_PNS(x = x, angular = angular, dims = dims, col = rainbow(n), method = "Main", pc = pc,tol = tol)
  if(plot){
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
  return(list(pns_data = simulation_2d$pns_data,MDS_info = simulation_2d$MDS_info))
}


# run MDS on polysphere, plot the first three principal components, and plot the residuals in the first 3 dimensions
plot_MDS_py <- function(x, dims = rep(1, ncol(x)), landmark_size = .1, angular = TRUE, pc = min(3, ncol(x)), tol = 5e-3, landmark_tol = 5e-4, gradient = FALSE) {
  n <- nrow(x)
  d_sphere <- sum(dims)
  simulation_2d <- MDS_PNS_py(x = x, angular = angular, dims = dims, col = rainbow(n), landmark_size = landmark_size, tol = tol, landmark_tol = landmark_tol, method = "Main", pc = pc, gradient= gradient)
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
  return(simulation_2d)
}

# Inverts PNS onto the polysphere using either Nadaraya-Watson or inverse MDS, takes in polysphere with cartesian coordinates

invert_MDS <- function(x, result, dims = rep(1, ncol(x)), grid_size =100, method = "NW", angular= TRUE, plot = TRUE){
  n <- nrow(x)
  d_sphere <- sum(dims)
  sc1 <- march_PNS1_shapes(seq(0, result$pns_info$radii[d_sphere] * 2 * pi, l = grid_size),
                           PNS = result$pns_info
  )
  comp_x <- result$MDS_data
  if (method == "NW"){
    nw_pc <- opt_mNW(x=sc1, X=comp_x, Y=x, dims = dims, K = dnorm)
    nw_angle <- polySphereToAngles(nw_pc,dims) 
    out <- nw_angle
  } else if(method == "both"){
    nw_pc <- opt_mNW(x=sc1, X=comp_x, Y=x, dims = dims, K = dnorm)
    nw_angle <- polySphereToAngles(nw_pc,dims)
    out <- nw_angle
    poly_pc_angle <- inverse_MDS(x, angular = FALSE, dims = dims, grid_size = grid_size,mds = result, plot_pc = FALSE)
  } else {
    poly_pc_angle <- inverse_MDS(x, angular = FALSE, dims = dims, grid_size = grid_size, mds=result, plot_pc = FALSE)
    out <- poly_pc_angle
  }
  if(plot){
    plot_inverse(x, nw_angle, poly_pc_angle, dims = dims, grid_size =grid_size, angular = angular, method = method)
  }
  return(out)
}

plot_inverse <- function(x, nw_pc, poly_pc, dims = rep(1, ncol(x)), grid_size =100, angular = TRUE, method = "NW"){
  if(angular){
    clusters <- polySphereToAngles(x,dims)
  } else{
    nw_pc <- polyAnglesToSphere(nw_pc,dims)
    poly_pc <- polyAnglesToSphere(poly_pc,dims)
  }
  if(method == "NW"){
    out <- nw_pc
    if(d_sphere < 3){
      plot(rbind(clusters,out),col=c(rainbow(n), rep("black", grid_size)),
           xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i")
    } else{
      plot3d(rbind(clusters,out),col=c(rainbow(n), rep("black", grid_size)))
    }
  } else if(method == "both"){
    if(d_sphere < 3){
      plot(rbind(clusters,nw,poly_pc),col=c(rainbow(n), rep("black", grid_size), rep("grey", grid_size)),
           xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i")
    } else{
      plot3d(rbind(clusters,nw_pc,poly_pc),col=c(rainbow(n), rep("black", grid_size), rep("grey", grid_size)))
    }
  }
  else {
    out <- poly_pc
    if(d_sphere < 3){
      plot(rbind(clusters,out),col=c(rainbow(n), rep("black", grid_size)),
           xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i")
    } else{
      plot3d(rbind(clusters,out),col=c(rainbow(n), rep("black", grid_size)))
    }
  }
}

# Run MDS on polysphere, return results
run_MDS_py <- function(x, dims = rep(1, ncol(x)), landmark_size = .1, angular = TRUE,  pc = min(3, ncol(x)), tol = 5e-3, landmark_tol = 5e-4, plot = FALSE, gradient = FALSE) {
  n <- nrow(x)
  d_sphere <- sum(dims)
  simulation_2d <- MDS_PNS_py(x = x, angular = angular, dims = dims, col = rainbow(n), landmark_size = landmark_size, tol = tol, landmark_tol = landmark_tol, method = "Main", pc = pc, gradient = gradient)
  #print(simulation_2d$MDS_info$stress_final)
  if(plot){
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
  }
  return(list(pns_data = simulation_2d$pns_data,MDS_info = simulation_2d$MDS_info))

}


# compare centroid and interpolation MDS accuracy and runtime with default
compare_MDS <- function(x,dims = rep(1, ncol(x)), landmark_size = .1, angular = TRUE, pc = min(3, ncol(x)), tol = 5e-3, landmark_tol = 5e-4,plot=FALSE,default=FALSE){
  n <- nrow(x)
  d_sphere <- sum(dims)

  time <- Sys.time()
  res_gradient <- run_MDS_py(x, landmark_size = landmark_size, dims=dims, angular=angular,plot=plot,gradient = TRUE)
  gradient_time <- Sys.time() - time
  
  time <- Sys.time()
  res_centroid <- run_MDS_py(x, landmark_size = landmark_size, dims=dims, angular=angular,plot=plot,gradient = FALSE)
  centroid_time <- Sys.time() - time
  
  # Return function values for residual and centroid MDS
  value_gradient <- res_gradient$MDS_info$stress_final
  value_centroids <- res_centroid$MDS_info$stress_final
  function_values <- c(value_gradient,value_centroids)
  # Account for minimum R^2 scores up to a flip
  # residual_centroids <- modified_R2(res_centroid$pns_data[,2],res$pns_data[,2])
  # score_residuals <- c(residual_centroids,residual_subset)
  # Keep track of runtimes
  times <- c(gradient_time,centroid_time)
  if(default){
    time <- Sys.time()
    res <- run_MDS(x,dims=dims, angular=angular,tol=tol,plot=plot)
    normal_time<- Sys.time() - time
    value_MDS <- res$MDS_info$stress_final
    function_values <- c(value_gradient,value_centroids,value_MDS)
    times <- c(gradient_time,centroid_time,normal_time)
    
  }
  return(c(function_values,times))
}
# Compare runtimes, function values, and R2 scores of interpolation MDS over multiple trials
compare_trials <- function(n,nspheres, landmark_size,trials, angular= FALSE, dims = rep(1,nspheres),default=FALSE){
  n <- n
  nspheres <- nspheres
  dims = dims
  angular = angular
  landmark_size = landmark_size
  if(default){
    res <- rep(0,6) 
  }else{
    res <- rep(0,4)
  }
  for(i in 1:trials){
    samp_1 <- r_path_s1r(n = n, r = nspheres, k = rep(1, nspheres), angles = TRUE)
    clusters <- samp_1
    trial_res <- compare_MDS(clusters, landmark_size = landmark_size, dims=dims, angular=angular, landmark_tol = 5e-4, tol = 5e-4, plot =FALSE,default=default)
    
    res <- res + trial_res
  }
  return(res/trials)
}
# Compare predicted 1D data with a given function using the formula RSS/TSS up to a flip 
modified_R2 <- function(y,f){
  n <- length(y)
  flip_y <- y[c(n:1)]
  result<-min(sum((y-f)^2),sum((flip_y-f)^2))/(n*var(f))
  return(result)
}