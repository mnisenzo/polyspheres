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
source("angles_to_spheres.R")



### Functions

## MDS & PNS for polysphere data (steps 1 and 2 of ST-PCA)
MDS_PNS <-
  function(x,
           init = "cmds",
           angular = TRUE,
           dims = rep(1, ncol(x)),
           col = NULL,
           method = "Main",
           tol = 5e-3,
           pc = min(3, ncol(x))) {
    # Dimensions of dataset
    n <- nrow(x)
    d <- ncol(x)
    d_sphere <- sum(dims)

    if (is.null(col)) {
      col <- rep("black", n)
    }

    # Initial guess for radius is number of spheres in polysphere
    r_star <- sqrt(length(dims)) # compare sum of squares

    # Compute polysphere distances [s^d_1, s^d_2...]
    ind_dj <- comp_ind_dj(dims)

    # Convert from hyperspherical to cartesian coordinates if needed
    cart_x <- x
    if (angular) {
      cart_x <- polyAnglesToSphere(x, dims)
    }
    # Go straight to upper/lower part
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
        y = c(y, r_star),
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
    sink("NUL")
    general <-
      pns(t(out), sphere.type = "seq.test") # PNS with sequential testing for type of fitting sphere
    sink()
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
      # plot(
      #   t(general$resmat),
      #   col = col,
      #   main = paste("First two Components, r* =", r_star),
      #   asp = 1
      # )
      pairs(
        t(general$resmat[1:ncol, ]),
        col = col,
        xlim = c(-pi, pi), ylim = c(-pi, pi),
        main = paste("First", ncol, "components, r* =", r_star, collapse = " ")
      )
    } else {
      # plot(
      #   t(general$resmat),
      #   col = col,
      #   main = paste("First two Components, r* =", r_star),
      #   asp = 1
      # )
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
