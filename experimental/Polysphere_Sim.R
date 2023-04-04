
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
plot_MDS(clusters,angular=TRUE)


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
plot_MDS(clusters,angular=TRUE)

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
plot_MDS(clusters,c(2,2),angular=TRUE)
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
plot_MDS(clusters,c(2,2),angular=TRUE)
