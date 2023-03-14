### Numerical Experiments of the main paper (Section 4)

source("C:/Users/minise/Downloads/STPCA/Script1_Setup_and_Functions.R")

### Section 4.1 T^2 clusters toy data

set.seed(100)
n <- 90
d <- 3
Sigma <- 0.1 * diag(c(1, 1))
cluster1 <- mvrnorm(n = n / 3, mu = c(-1, -2), Sigma = Sigma)
cluster2 <- mvrnorm(n = n / 3, mu = c(3, 0.4), Sigma = Sigma)
cluster3 <- mvrnorm(n = n / 3, mu = c(-0.8, 2.5), Sigma = Sigma)
clusters <- rbind(cluster1, cluster2, cluster3)
clusters <- (clusters + pi) %% (2 * pi) - pi


a = anglesToSphere(clusters)- polyAnglesToSphere(clusters,c(2))

col1 <- rep("#009999", (n / 3))
col2 <- rep("#FF0033", (n / 3))
col3 <- rep("#66FF00", (n / 3))
col <- c(col1, col2, col3)



plot(clusters, col = col)

simulation_2d <- MDS_PNS(x = clusters, col = col, method = "Main")
percent_var_explained(simulation_2d, x = clusters)

d_sphere <- 2
sc <- march_PNS1_shapes(seq(0, simulation_2d$pns_info$radii[d_sphere] * 2 * pi, l = 100),
  PNS = simulation_2d$pns_info
)

pred <- prediction(
  y_new = sc, x_old = clusters, y_old = simulation_2d$MDS_data,
  r_star = simulation_2d$MDS_info$rstar, d_sphere = d_sphere
)

## Plotting

all_torus <- rbind(clusters, pred)
all_sphere <- rbind(simulation_2d$MDS_data, sc)
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

# Figure 2.a
plot(all_torus,
  col = c(col, rep("black", 100)), xlab = expression(paste(theta)), ylab = expression(paste(phi)),
  xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE, xaxs = "i", pch = 20, yaxs = "i"
)
linesTorus(x_interp_wra[, 1], x_interp_wra[, 2], ltyCross = 3, col = 2)
points(x_interp_wra, col = 1)
sdetorus::torusAxis()

# Figure 2.c
plot(rbind(simulation_2d$pns_data, cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),
  col = c(col, rep("black", 100)), xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE,
  xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
)
points(lin_interp(cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)), N = 10))
sdetorus::torusAxis()

plot3d(rbind(simulation_2d$MDS_data, sc),
  col = c(rainbow(n), rep("black", 100)),
  xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), box = FALSE, axes = FALSE,
  xlab = "", ylab = "", zlab = "", pch = 19
)
sphere1.f(x0 = 0, y0 = 0, z0 = 0, r = 1, col = "azure2", lit = FALSE, alpha = 0.5)

### Section 4.2: T^3 Clusters Toy Data

set.seed(100)
n <- 90
d <- 3
Sigma <- 0.1 * diag(c(1, 1, 1))
cluster1 <- mvrnorm(n = n / 3, mu = c(-1, -2, 0), Sigma = Sigma)
cluster2 <- mvrnorm(n = n / 3, mu = c(1, 1, -1), Sigma = Sigma)
cluster3 <- mvrnorm(n = n / 3, mu = c(0, -1, 3), Sigma = Sigma)
clusters <- rbind(cluster1, cluster2, cluster3)
clusters <- (clusters + pi) %% (2 * pi) - pi
col1 <- rep("#009999", (n / 3))
col2 <- rep("#FF0033", (n / 3))
col3 <- rep("#66FF00", (n / 3))
col <- c(col1, col2, col3)
plot3d(x = clusters, col = col)

## Running MDS-PNS
simulation <- MDS_PNS(x = clusters, col = col)

## Prediction with 100 points
d_sphere <- 3
sc <- march_PNS1_shapes(seq(0, simulation$pns_info$radii[d_sphere] * 2 * pi, l = 10),
  PNS = simulation$pns_info
)

pred <- prediction(
  y_new = sc, x_old = clusters, y_old = simulation$MDS_data,
  r_star = simulation$MDS_info$rstar, d_sphere = d_sphere
)

percent_var_explained(simulation, clusters)

## Plotting

all_torus <- rbind(clusters, pred)
all_sphere <- rbind(simulation$MDS_data, sc)
group <- c(rep(1, n), rep(2, 100))

x_unwrapped <- unwrapCircSeries(x = rbind(pred, pred[1, ]))
plot3d(x_unwrapped)
x_interp <- lin_interp(rbind(x_unwrapped), N = 10)
plot3d(x_interp)
lines3d(x_interp, col = 2, type = "points")
points3d(x_unwrapped, pch = 16)
# Wrap the interpolated data
plot3d(pred, axes = FALSE, xlab = "", ylab = "") # Irregular data
sdetorus::torusAxis()
x_interp_wra <- toPiInt(x_interp)
# Nice plot
plot3d(x_interp_wra)
linesTorus(x_interp_wra[, 1], x_interp_wra[, 2], ltyCross = 3, col = 2)
points(x_interp_wra, col = 1)
points(pred, pch = 16)

# Figure 3.a
plot3d(all_torus,
  col = c(col, rep("black", 100)),
  xlim = c(-pi, pi), ylim = c(-pi, pi), zlim = c(-pi, pi), xlab = "",
  ylab = "", zlab = "", box = TRUE,
  type = "p", axes = FALSE
)
box3d(expand = 1)
linesTorus3d(x_interp_wra[, 1], x_interp_wra[, 2], x_interp_wra[, 3], ltyCross = 3, col = 2)
points3d(x_interp_wra, col = 1, pch = 19, size = 6)
axes3d(edges = "x", labels = FALSE, tick = FALSE)
axis3d("x", at = c(-pi, -pi / 2, 0, pi / 2, pi), labels = expression(-pi, -pi / 2, 0, pi / 2, pi))
axis3d("y", at = c(-pi, -pi / 2, 0, pi / 2, pi), labels = expression(-pi, -pi / 2, 0, pi / 2, pi))
axis3d("z", at = c(-pi, -pi / 2, 0, pi / 2, pi), labels = expression(-pi, -pi / 2, 0, pi / 2, pi))

# Figure 3.c
plot(rbind(simulation$pns_data[, 1:2], cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),
  col = c(col, rep("black", 100)), xlim = c(-pi, pi), ylim = c(-pi / 2, pi / 2), axes = FALSE,
  xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
)
points(lin_interp(cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)), N = 10))
sdetorus::torusAxis(sides = 1)
axis(side = 2, at = c(-pi / 2, -pi / 4, 0, pi / 4, pi / 2), labels = expression(-pi / 2, -pi / 4, 0, pi / 4, pi / 2))

plot3d(anglesToSphere(rbind(simulation$pns_data[, 1:2], cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)))),
  col = c(col, rep("black", 100)), xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), box = FALSE, axes = FALSE,
  xlab = "", ylab = "", zlab = "", pch = 19
)
sphere1.f(x0 = 0, y0 = 0, z0 = 0, r = 1, col = "azure2", lit = FALSE, alpha = 0.5)

# all_sphere = anglesToSphere(rbind(simulation$pns_data[,1:2],cbind(seq(from = -pi, to = pi, by = 2*pi/99), rep(0, 100))))

### Section 4.3: T^2 Unisotropic Wrapped Gaussian

set.seed(1)
n <- 50
d <- 2
Sigma <- rbind(c(2, 1.5), c(2, 3))
clusters <- mvrnorm(n = n, mu = c(-1, 0), Sigma)
clusters <- (clusters + pi) %% (2 * pi) - pi
col <- rep("black", n)

plot(clusters, col = col)

simulation_wrapped <- MDS_PNS(x = clusters, col = col, method = "Main")
d_sphere <- 2
percent_var_explained(simulation_wrapped, clusters)

sc <- march_PNS1_shapes(seq(0, simulation_wrapped$pns_info$radii[d_sphere] * 2 * pi, l = 100),
  PNS = simulation_wrapped$pns_info
)

pred <- prediction(
  y_new = sc, x_old = clusters, y_old = simulation_wrapped$MDS_data,
  r_star = simulation_wrapped$MDS_info$rstar, d_sphere = d_sphere
)

all_torus <- rbind(clusters, pred)
all_sphere <- rbind(simulation_wrapped$MDS_data, sc)

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

## Plotting

# Figure 4.a
plot(all_torus,
  col = c(col, rainbow(100)), xlab = "", ylab = "",
  xlim = c(-pi, pi), ylim = c(-pi, pi),
  pch = 20, axes = FALSE, xaxs = "i", yaxs = "i"
)
linesTorus(x_interp_wra[, 1], x_interp_wra[, 2], ltyCross = 3, col = 2)
points(x_interp_wra, col = rainbow(5001))
sdetorus::torusAxis()

# Figure 4.c
plot(rbind(simulation_wrapped$pns_data, cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),
  col = c(rainbow(n), rep("black", 100)), xlim = c(-pi, pi), ylim = c(-pi / 2, pi / 2), axes = FALSE,
  xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
)
points(lin_interp(cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)), 10), col = rainbow(991))
sdetorus::torusAxis(sides = 1)
axis(side = 2, at = c(-pi / 2, -pi / 4, 0, pi / 4, pi / 2), labels = expression(-pi / 2, -pi / 4, 0, pi / 4, pi / 2))

plot3d(rbind(simulation_wrapped$MDS_data, sc),
  col = c(rainbow(n), rep("black", n)),
  xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), box = FALSE, axes = FALSE,
  xlab = "", ylab = "", zlab = "", pch = 19, size = 4
)
sphere1.f(x0 = 0, y0 = 0, z0 = 0, r = 1, col = "azure2", lit = FALSE, alpha = 0.5)
