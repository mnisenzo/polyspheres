### Supplementary Material - Further simulations

source("C:/Users/minise/Downloads/STPCA/Script1_Setup_and_Functions.R")
source("C:/Users/minise/Downloads/STPCA/Script5_Modification_of_function_for_SM.R")


## Sum Subsection B.1: robustness to outliers
## when tolerance high it's robust, balance the speed. recommendation when there is concerns
## about outliers, more precision is recommended.
n <- 300
unif_square1 <- runif(n, min = 3 * pi / 4, max = 5 * pi / 4)
unif_square2 <- runif(n, min = 3 * pi / 4, max = 5 * pi / 4)
unif_square <- cbind(unif_square1, unif_square2)
outliers <- rbind(c(pi / 4, pi / 3), c(pi / 4, pi), c(pi / 4, 5 * pi / 3))
outlier_data <- rbind(unif_square, outliers)

outlier_data <- outlier_data - pi
plot(outlier_data, xlim = c(-pi, pi), ylim = c(-pi, pi), col = c(rep("black", n), rep("red", 3)))

Analysis <- MDS_PNS(outlier_data, col = c(rep("black", n), rep("red", 3)), tol = 1e-4)

# Generate plots

# SM Figure 1a
plot(outlier_data,
  col = c(rep("black", 300), rep("red", 3)), xlab = expression(paste(theta)), ylab = expression(paste(phi)),
  xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE, xaxs = "i", pch = 20, yaxs = "i"
)
sdetorus::torusAxis()

Analysis$pns_data[, 1] <- -Analysis$pns_data[, 1]
# SM Figure 1c
plot(Analysis$pns_data,
  col = c(rep("black", 300), rep("red", 3)), xlim = c(-pi, pi), ylim = c(-pi / 2, pi / 2), axes = FALSE,
  xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
)
sdetorus::torusAxis(sides = 1)
axis(side = 2, at = c(-pi / 2, -pi / 4, 0, pi / 4, pi / 2), labels = expression(-pi / 2, -pi / 4, 0, pi / 4, pi / 2))

## SM Subsection B.2: higher dimensional data

set.seed(100)
n <- 300
d <- 50
Sigma <- 0.1 * diag(rep(1, 50))
cluster1 <- mvrnorm(n = n / 2, mu = c(-1, rep(0, 49)), Sigma = Sigma)
cluster2 <- mvrnorm(n = n / 2, mu = c(1, rep(0, 49)), Sigma = Sigma)
clusters <- rbind(cluster1, cluster2)
clusters <- (clusters + pi) %% (2 * pi) - pi
col1 <- rep("#009999", (n / 2))
col2 <- rep("#FF0033", (n / 2))
col <- c(col1, col2)
simul_50 <- MDS_PNS(x = clusters, col = col)

## Generate the plots

# SM Figure 2a
plot(
  x = clusters[, 1], y = clusters[, 2], col = col, xlab = expression(paste(theta)), ylab = expression(paste(phi)),
  xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE, xaxs = "i", pch = 20, yaxs = "i"
)
sdetorus::torusAxis()

# SM Figure 2c
simul_50$pns_data[, 1] <- -simul_50$pns_data[, 1]
plot(
  x = simul_50$pns_data[, 1], y = simul_50$pns_data[, 2],
  col = col, xlim = c(-pi, pi), ylim = c(-pi / 2, pi / 2), axes = FALSE,
  xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
)
sdetorus::torusAxis(sides = 1)
axis(side = 2, at = c(-pi / 2, -pi / 4, 0, pi / 4, pi / 2), labels = expression(-pi / 2, -pi / 4, 0, pi / 4, pi / 2))


d_sphere <- 50
sc <- march_PNS1_shapes(seq(0, simul_50$pns_info$radii[d_sphere] * 2 * pi, l = 100),
  PNS = simul_50$pns_info
)
simul50_allsphere <- anglesToSphere(rbind(simul_50$pns_data[, 1:2], cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))))

plot3d(simul50_allsphere,
  col = c(col, rep("black", 100)), xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), box = FALSE, axes = FALSE,
  xlab = "", ylab = "", zlab = "", pch = 19
)
sphere1.f(x0 = 0, y0 = 0, z0 = 0, r = 1, col = "azure2", lit = FALSE, alpha = 0.5)


### SM Subsection B.5: contour analysis

## Generate the data
set.seed(1)
n <- 300
d <- 2
Sigma <- 0.1 * diag(c(1, 1))
bivar1 <- mvrnorm(n = n / 3, mu = c(-1, -2), Sigma = Sigma / 2)
bivar2 <- mvrnorm(n = n / 3, mu = c(3, 0.5), Sigma = Sigma / 2)
bivar3 <- mvrnorm(n = n / 3, mu = c(-0.8, 2.5), Sigma = Sigma / 2)
x3 <- rbind(bivar1, bivar2, bivar3)
x3 <- (x3 + pi) %% (2 * pi) - pi
plot(x3)
col3 <- c(rep("#009999", (n / 3)), rep("#FF0033", (n / 3)), rep("#66FF00", (n / 3)), rainbow(100))

ST_PCA_contour <- MDS_PNS(x3, col = c(rep("black", n), rep("red", 3)), tol = 1e-4)
d_sphere <- 2
sc <- march_PNS1_shapes(seq(0, ST_PCA_contour$pns_info$radii[d_sphere] * 2 * pi, l = 100),
  PNS = ST_PCA_contour$pns_info
)
## Figure 4c will be saved as a pdf titled simul_2d_contour_56
pred <- prediction(
  y_new = sc, x_old = x3, y_old = ST_PCA_contour$MDS_data,
  r_star = ST_PCA_contour$MDS_info$rstar, contour = 1, col = col3
)


all_torus <- rbind(x3, pred)
all_sphere <- rbind(ST_PCA_contour$MDS_data, sc)
group <- c(rep(1, n), rep(2, 100))

x_unwrapped <- unwrapCircSeries(x = rbind(pred, pred[1, ]))
plot(x_unwrapped)
x_interp <- lin_interp(rbind(x_unwrapped), N = 50)
lines(x_interp, col = 2, type = "o")
points(x_unwrapped, pch = 16)
# Wrap the interpolated data
plot(pred, axes = FALSE, xlab = "", ylab = "") # Irregular data
sdetorus::torusAxis()
x_interp_wra <- toPiInt(x_interp)
# Nice plot
linesTorus(x_interp_wra[, 1], x_interp_wra[, 2], ltyCross = 3, col = 2)
points(x_interp_wra, col = 1)
points(pred, pch = 16)


# SM Figure 4a
plot(all_torus,
  col = col3, xlab = expression(paste(theta)), ylab = expression(paste(phi)),
  xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE, xaxs = "i", pch = 20, yaxs = "i"
)
linesTorus(x_interp_wra[, 1], x_interp_wra[, 2], ltyCross = 3, col = 2)
points(x_interp_wra, col = 1)
sdetorus::torusAxis()


# SM Figure 4b
plot(rbind(ST_PCA_contour$pns_data, cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),
  col = col3, xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE,
  xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
)
points(lin_interp(cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)), N = 10))
sdetorus::torusAxis()
