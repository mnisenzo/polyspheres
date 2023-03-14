### Numerical Experiments of the main paper (Section 4)

source("C:/Users/minise/Downloads/STPCA/Script1_Setup_and_Functions.R")
source("C:/Users/minise/Downloads/samplers.R")

#' # Straight trends on (S^1)^2
set.seed(5)
n <- 20
d_sphere <- 2
samp_1 <- r_path_s1r(n = n, r = 2, k = c(1, 1), angles = TRUE)
plot(samp_1, xlim = c(-pi, pi), ylim = c(-pi, pi), col = rainbow(n),
     axes = FALSE, xlab = "", ylab = "", pch = 16)
sdetorus::torusAxis()
clusters = samp_1
simulation_2d <- MDS_PNS(x = clusters, method = "Main",col=rainbow(n))
percent_var_explained(simulation_2d, x = clusters)

plot(rbind(simulation_2d$pns_data, cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),
     col = c(rainbow(n), rep("black", 100)), xlim = c(-pi, pi), ylim = c(-pi/2, pi/2), axes = FALSE,
     xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
)
points(lin_interp(cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)), N = 10))
sdetorus::torusAxis()


#' # Straight trends on (S^1)^3
set.seed(100)
n <- 20
d_sphere <- 3
samp_2 <- r_path_s1r(n = n, r = 3,k = c(0, 1,2), angles = TRUE)
pairs(samp_2, xlim = c(-pi, pi), ylim = c(-pi, pi), col = rainbow(n),
      pch = 16)
sdetorus::torusAxis()
scatterplot3d::scatterplot3d(
  samp_2, xlim = c(-pi, pi), ylim = c(-pi, pi), zlim = c(-pi, pi),
  xlab = "", ylab = "", zlab = "", color = rainbow(n), pch = 16
)
clusters = samp_2
simulation_2d <- MDS_PNS(x = clusters, method = "Main",col=rainbow(n))
percent_var_explained(simulation_2d, x = clusters)



plot(rbind(simulation_2d$pns_data[, c(1,2)], cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),
     col = c(rainbow(n), rep("black", 100)), xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE,
     xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
)
points(lin_interp(cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)), N = 10))
sdetorus::torusAxis()

sc <- march_PNS1_shapes(seq(0, simulation_2d$pns_info$radii[d_sphere] * 2 * pi, l = 100),
                        PNS = simulation_2d$pns_info
)

plot3d(rbind(simulation_2d$MDS_data, sc),
       col = c(rainbow(n), rep("black", 100)),
       xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), box = FALSE, axes = FALSE,
       xlab = "", ylab = "", zlab = "", pch = 19
)
sphere1.f(x0 = 0, y0 = 0, z0 = 0, r = 1, col = "azure2", lit = FALSE, alpha = 0.5)

#' # Small-circle trends on (S^2)^2
set.seed(20)
d_sphere <- 4
n <- 60
samp_3 <- r_path_s2r(n = n, r = 2, sigma = 0.05)
scatterplot3d::scatterplot3d(
  samp_3[, , 1], xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  xlab = "", ylab = "", zlab = "", color = rainbow(n), pch = 16
)
scatterplot3d::scatterplot3d(
  samp_3[, , 2], xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  xlab = "", ylab = "", zlab = "", color = rainbow(n), pch = 16
)
clusters = cbind(samp_3[, , 1],samp_3[, , 2])
simulation_2d <- MDS_PNS(x = clusters, angular=FALSE,dims=c(2,2),col = rainbow(n), method = "Main")
percent_var_explained(simulation_2d, x = clusters)

sc <- march_PNS1_shapes(seq(0, simulation_2d$pns_info$radii[d_sphere] * 2 * pi, l = 100),
                        PNS = simulation_2d$pns_info
)

plot(rbind(simulation_2d$pns_data[, c(1,3)], cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),
     col = c(rainbow(n), rep("black", 100)), xlim = c(-pi, pi), ylim = c(-pi / 2, pi / 2), axes = FALSE,
     xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
)
points(lin_interp(cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)), N = 10))
sdetorus::torusAxis(sides = 1)
axis(side = 2, at = c(-pi, -pi / 2, 0, pi / 2, pi), labels = expression(-pi / 2, -pi / 4, 0, pi / 4, pi / 2))

plot3d(rbind(simulation_2d$MDS_data, sc)[,1:3],
       col = c(rainbow(n), rep("black", 100)),
       xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), box = FALSE, axes = FALSE,
       xlab = "", ylab = "", zlab = "", pch = 19
)
sphere1.f(x0 = 0, y0 = 0, z0 = 0, r = 1, col = "azure2", lit = FALSE, alpha = 0.5)

#' # Spiral trends on (S^2)^2
n <- 60
d_sphere <- 4

samp_4 <- r_path_s2r(n = n, r = 2, c = 3, spiral = TRUE, sigma = 0.05)

scatterplot3d::scatterplot3d(
  samp_4[, , 1], xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  xlab = "", ylab = "", zlab = "", color = rainbow(n), pch = 16
)
scatterplot3d::scatterplot3d(
  samp_4[, , 2], xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
  xlab = "", ylab = "", zlab = "", color = rainbow(n), pch = 16
)

clusters = cbind(samp_4[, , 1],samp_4[, , 2])
simulation_2d <- MDS_PNS(x = clusters, angular=FALSE,dims=c(2,2),col = rainbow(n), method = "Main")
percent_var_explained(simulation_2d, x = clusters)

plot(rbind(simulation_2d$pns_data[, c(1,2)], cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100))),
     col = c(rainbow(n), rep("black", 100)), xlim = c(-pi, pi), ylim = c(-pi / 2, pi / 2), axes = FALSE,
     xlab = "", ylab = "", pch = 19, xaxs = "i", yaxs = "i"
)
points(lin_interp(cbind(seq(from = -pi, to = pi, by = 2 * pi / 99), rep(0, 100)), N = 10))
sdetorus::torusAxis(sides = 1)
axis(side = 2, at = c(-pi, -pi / 2, 0, pi / 2, pi), labels = expression(-pi / 2, -pi / 4, 0, pi / 4, pi / 2))

plot3d(rbind(simulation_2d$MDS_data, sc),
       col = c(rainbow(n), rep("black", 100)),
       xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), box = FALSE, axes = FALSE,
       xlab = "", ylab = "", zlab = "", pch = 19
)
sphere1.f(x0 = 0, y0 = 0, z0 = 0, r = 1, col = "azure2", lit = FALSE, alpha = 0.5)
