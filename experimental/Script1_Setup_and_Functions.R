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
library(MASS)
library(sdetorus)
library(rgl)
library(plot3D)
library(testthat)
library(usethis)


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
    # print(as.matrix(theta[,(ind_p+1):(ind_p+d)]))
    # print(anglesToSphere(as.matrix(theta[,(ind_p+1):(ind_p+d)])))
    if (nrow(theta) == 1) {
      input <- t(input)
    }
    x[, (ind_r + 1):(ind_r + d + 1)] <- anglesToSphere(input)
    ind_p <- ind_p + d
    ind_r <- ind_r + d + 1
  }
  return(x)
}


## MDS & PNS for toroidal data (steps 1 and 2 of ST-PCA)
MDS_PNS <-
  function(x,
           init = "cmds",
           angular = TRUE,
           dims = rep(1, ncol(x)),
           col = NULL,
           method = "Main",
           tol = 5e-2) {
    # dimensions of dataset
    n <- nrow(x)
    d <- ncol(x)
    d_sphere <- d

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
    }
    # Main optimization
    if (method == "Main") {
      # Objective function
      stress <- function(y, d_ij_poly, n, d_sphere, r_star) {
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
        return(sum((d_ij_poly - r_star * d_ij_sphere)^2) / length(d_ij_poly))
      }
      # Optimization
      res <- optim(
        fn = stress,
        par = c(y, r_star),
        d_ij_poly = d_ij_polysphere,
        n = n,
        d_sphere = d_sphere,
        method = "BFGS",
        control = list(reltol = tol)
      )


      # Visualize output
      stress_initial <- stress(
        y = y,
        d_ij_poly = d_ij_polysphere,
        n = n,
        d_sphere = d_sphere,
        r = r_star
      )
      r.hat <- res$par[n * (d_sphere + 1) + 1]
      MDS_info <- list(
        stress_init = stress_initial,
        stress_final = res$value,
        rstar = r_star,
        r.hat = r.hat
      )
      out <-
        matrix(res$par[1:(n * (d_sphere + 1))], nrow = n, ncol = d_sphere + 1)
      out <- out / sqrt(rowSums(out^2))
    }
    general <-
      pns(t(out), sphere.type = "seq.test") # PNS with sequential testing for type of fitting sphere

    ## Plotting
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
    } else {
      plot(
        t(general$resmat),
        col = col,
        main = paste("First two Components, r* =", r_star),
        asp = 1
      )
      pairs(
        t(general$resmat[1:3, ]),
        col = col,
        main = paste("First three Components, r* =", r_star)
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

## Prediction (Step 3 of ST-PCA)
prediction <- function(x_old, y_old, y_new, r_star, d_sphere) {
  if (is.matrix(y_new)) {
    n_new <- nrow(y_new)
    d <- ncol(y_new)
  } else {
    n_new <- 1
    d <- length(y_new)
  }
  n_old <- nrow(y_old)
  x_new_final <- matrix(0, nrow = n_new, ncol = d - 1)

  # Objective function
  obj <- function(x_new_torus, x_old_torus, n, r, dist) {
    x_new_torus <- (x_new_torus + pi) %% (2 * pi) - pi
    # Calculate torus distances
    x <- rbind(x_old_torus, x_new_torus)
    # Spherical Distances for full dataset
    ind_tri <- sphunif::upper_tri_ind(n)
    x_array <- array(dim = c(n, 1, d - 1))
    x_array[, 1, ] <- x
    d_ij_circ_each_dim <- sphunif::Psi_mat(
      data = x_array,
      use_ind_tri = TRUE,
      ind_tri = ind_tri
    )
    d_ij_torus <- sqrt(rowSums(d_ij_circ_each_dim^2))

    dist_torus[upper.tri(dist_torus)] <-
      d_ij_torus # Sets the Upper triangular part of the matrix
    dist_torus <- dist_torus + t(dist_torus)
    dist2 <- dist_torus[, n_old + 1][1:n_old]

    # Stress function (only sum the distinct distances and disregard diagonal)
    return(2 * sum((r * dist3 - dist2)^2) / n)
  }

  for (i in c(1:n_new)) {
    if (n_new == 1) {
      y <- rbind(y_old, y_new)
    } else {
      y <- rbind(y_old, y_new[i, ])
    }
    # Calculate Spherical Distances
    y <- y / sqrt(rowSums(y^2))
    dim(y) <- c(dim(y), 1)
    ind_tri <- sphunif::upper_tri_ind(n_old + 1)
    d_ij_sphere <- sphunif::Psi_mat(
      data = y,
      use_ind_tri = TRUE,
      ind_tri = ind_tri
    )
    d_ij_sphere <- drop(d_ij_sphere)

    dist <- matrix(0, nrow = n_old + 1, ncol = n_old + 1)
    dist[upper.tri(dist)] <-
      d_ij_sphere # Sets the Upper triangular part of the matrix
    dist <- dist + t(dist) # adds the lower to make it symmetric
    dist3 <- dist[, n_old + 1][1:n_old]
    dist_torus <- matrix(0, nrow = n_old + 1, ncol = n_old + 1)

    # Initialization
    if (i == 1) {
      x_new_torus <- x_old[which.min(dist3), ]
    } else {
      x_new_torus <- x_new_final[i - 1, ]
    }

    # Optimization
    pred <- optim(
      fn = obj,
      par = c(x_new_torus),
      x_old_torus = x_old,
      n = n_old + 1,
      r = r_star,
      dist = dist3,
      method = "BFGS",
      control = list(trace = 3)
    )$par
    pred <- (pred + pi) %% (2 * pi) - pi
    x_new_final[i, ] <- pred
  }
  return(x_new_final)
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
## Percent of variance explained by each component
percent_var_explained <-
  function(MDSPNS, x, dims = rep(1, ncol(x))) {
    d_sphere <- dim(x)[2]
    testsphere <- array(NA, dim = c(
      length(MDSPNS$MDS_data[, 1]), d_sphere + 1,
      d_sphere
    ))
    testtorus <- array(NA, dim = c(
      length(MDSPNS$MDS_data[, 1]), d_sphere,
      d_sphere
    ))
    distances <-
      matrix(NA, nrow = length(MDSPNS$MDS_data[, 1]), ncol = d_sphere)
    for (d in 1:d_sphere) {
      if (d < d_sphere) {
        d_s <- d_sphere - d
        testsphere[, , d] <- t(shapes::PNSe2s(rbind(
          t(MDSPNS$pns_data)[1:d_s, ],
          matrix(0,
            ncol = ncol(t(
              MDSPNS$pns_data
            )),
            nrow = d
          )
        ), MDSPNS$pns_info))
        if (d == 1) {
          testtorus[, , d] <- prediction(
            y_new = testsphere[, , d],
            x_old = x,
            y_old = MDSPNS$MDS_data,
            r_star = MDSPNS$MDS_info$r.hat,
            d_sphere = d_sphere
          )
        } else {
          testtorus[, , d] <- prediction(
            y_new = testsphere[, , d],
            x_old = testtorus[, , d - 1],
            y_old = testsphere[, , d - 1],
            r_star = MDSPNS$MDS_info$r.hat,
            d_sphere = d_sphere
          )
        }
      } else {
        testsphere[1, , d] <- MDSPNS$pns_info$mean
        testtorus[1, , d] <- prediction(
          y_new = testsphere[1, , d],
          x_old = testtorus[, , d - 1],
          y_old = MDSPNS$MDS_data,
          r_star = MDSPNS$MDS_info$r.hat,
          d_sphere = d_sphere
        )
        # had x_old = x.
      }
      # next can be done with sapply
      for (j in seq_along(MDSPNS$MDS_data[, 1])) {
        if (d == 1) {
          # distances[j, d] <- torus_distance(x[j, ], testtorus[j, , 1])
          print(c(
            x[j, ],
            testtorus[j, , 1],
            torus_distance(x[j, ], testtorus[j, , 1]),
            polysphere_distance(x[j, ], testtorus[j, , 1], dims = dims)
          ))
          distances[j, d] <-
            polysphere_distance(x[j, ], testtorus[j, , 1], dims = dims)
        } else if (d < d_sphere) {
          # distances[j, d] <- torus_distance(
          #   testtorus[j, , d - 1],
          #   testtorus[j, , d]
          # )

          distances[j, d] <- polysphere_distance(testtorus[j, , d - 1],
            testtorus[j, , d],
            dims = dims
          )
        } else {
          # distances[j, d] <- torus_distance(
          #   testtorus[j, , d - 1],
          #   testtorus[1, , d]
          # )
          distances[j, d] <- polysphere_distance(testtorus[j, , d - 1],
            testtorus[1, , d],
            dims = dims
          )
        }
      }
    }
    var.total <- sum(distances^2)
    var.comp <- rep(NA, length(distances[1, ]))
    for (j in seq_along(distances[1, ])) {
      var.comp[j] <- sum(distances[, length(distances[1, ]) - j + 1]^2) /
        var.total
    }

    return(var.comp)
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
plot_MDS <- function(x, dims=rep(1,ncol(x)), angular=TRUE){
  n = nrow(x)
  d_sphere = sum(dims)
  simulation_2d <- MDS_PNS(x = x, angular=angular,dims=dims,col = rainbow(n), method = "Main")
  percent_var_explained(simulation_2d, x = x)
  
  sc <- march_PNS1_shapes(seq(0, simulation_2d$pns_info$radii[d_sphere] * 2 * pi, l = 100),
                          PNS = simulation_2d$pns_info
  )
  
  if (sum(dims) < 3){
    plot(rbind(simulation_2d$pns_data, cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),
         col = c(rainbow(n), rep("black", 100)), xlim = c(-pi, pi), ylim = c(-pi/2, pi/2), axes = FALSE,
         xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
    )
    points(lin_interp(cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)), N = 10))
    sdetorus::torusAxis()
  } else{
    # compare first two principal components
    plot(rbind(simulation_2d$pns_data[, c(1,2)], cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),
         col = c(rainbow(n), rep("black", 100)), xlim = c(-pi, pi), ylim = c(-pi / 2, pi / 2), axes = FALSE,
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
  
  
}
