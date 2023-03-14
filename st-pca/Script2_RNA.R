### Analysis of the Small RNA Dataset
### Section 5.2 of the main Paper

## Load data and preprocessing: as.numeric
## order clusters; convert degrees to radians; cluster membership

source("C:/Users/minise/Downloads/STPCA/Script1_Setup_and_Functions.R")

small_rna <- read.csv("small_rna.csv", header = FALSE, dec = ".")
for (i in 1:190) {
  for (j in 1:9) {
    small_rna[i, j] <- as.numeric(small_rna[i, j])
  }
}
cluster1 <- small_rna[small_rna[, 9] > -60, ]
cluster2 <- small_rna[(small_rna[, 8] > 0) & (small_rna[, 9] < -60), ]
cluster5 <- small_rna[(small_rna[, 8] < 0) & (small_rna[, 9] < -60), ]
clusters_ordered <- rbind(cluster1, cluster2, cluster5)

col1 <- rep("green", 59)
col2 <- rep("red", 88)
col3 <- rep("blue", 43)
col <- c(col1, col2, col3)
plot(clusters_ordered[, 8:9], col = col) # Plot pseudotorsion angles

small_rna_angles <- clusters_ordered[, 1:7]
small_rna_torsion <- clusters_ordered[, 8:9]
small_rna_matrix <- as.matrix(small_rna_angles)
small_rna_matrix <- (small_rna_matrix / 180) * pi # degrees to radians
rna_small_analysis <- MDS_PNS(x = small_rna_matrix, col = col)
plot(rna_small_analysis$pns_data[, 1:2], col = col)
d_sphere <- 7
rna_small_analysis$MDS_info$r.hat # report optimal radius

rna_final <- rna_small_analysis$pns_data # Save scores
membership <- c(rep(1, 59), rep(2, 88), rep(3, 43))
rna_final <- cbind(rna_final, membership) # assign membership to clusters
rna_final[, 1] <- -rna_final[, 1] ## Reordering 1st component to match T-PCA

## Auxiliary functions to run kernel mean shift clustering
eta_hat <- function(x, h, data) {
  h^2 * ks::kdde(
    x = data, h = h, eval.points = x, deriv.order = 1,
    supp = 50, binned = FALSE
  )$estimate /
    ks::kde(
      x = data, h = h, eval.points = x, supp = 50,
      binned = FALSE
    )$estimate
}
euler <- function(x, N = 1e6, verbose = FALSE, epsilon = 1e-5, eta_hat_efic) {
  for (i in 1:N) {
    x_new <- sdetorus::toPiInt(x + eta_hat_efic(x = sdetorus::toPiInt(x)))
    if (verbose) print(x_new)
    if (abs(sdetorus::toPiInt(x_new - x)) < epsilon) {
      if (verbose) cat("Stopped at", i, "\n")
      break
    }
    x <- x_new
  }
  return(x)
}

# Load scores
scores_STPCA <- as.data.frame(rna_final)
names(scores_STPCA)[1] <- "X"
scores_TPCA_SI <- read.csv("tpca_rna_SI.csv", sep = ",", header = TRUE)
scores_TPCA_SO <- read.csv("tpca_rna_SO.csv", sep = ",", header = TRUE)

# Homogenize data: transform to (-pi, pi) flip signs to have better
# comparisons, and rename variables
scores_TPCA_SO[, 1:7] <- -scores_TPCA_SO[, 1:7] / 180 * pi
scores_TPCA_SI[, 1:7] <- -scores_TPCA_SI[, 1:7] / 180 * pi
names(scores_TPCA_SO)[7:8] <- names(scores_TPCA_SI)[7:8] <- c("X", "membership")

# Reparametrize class memberships
scores_TPCA_SO$membership[scores_TPCA_SO$membership == 5] <- 3
scores_TPCA_SI$membership[scores_TPCA_SI$membership == 5] <- 3

# Sanity check: same cluster memberships
table(scores_STPCA$membership)
table(scores_TPCA_SO$membership)
table(scores_TPCA_SI$membership)

# Obtain optimal plug-in bandwidth for density derivative estimation,
# disregarding periodicity
h_STPCA <- ks::hpi(scores_STPCA$X, deriv.order = 1)
h_TPCA_SO <- ks::hpi(scores_TPCA_SO$X, deriv.order = 1)
h_TPCA_SI <- ks::hpi(scores_TPCA_SI$X, deriv.order = 1)

# Hack to have periodicity
samp_STPCA <- c(
  scores_STPCA$X - 2 * pi,
  scores_STPCA$X,
  scores_STPCA$X + 2 * pi
)
samp_TPCA_SO <- c(
  scores_TPCA_SO$X - 2 * pi,
  scores_TPCA_SO$X,
  scores_TPCA_SO$X + 2 * pi
)
samp_TPCA_SI <- c(
  scores_TPCA_SI$X - 2 * pi,
  scores_TPCA_SI$X,
  scores_TPCA_SI$X + 2 * pi
)

# Kdes of data and class memberships
par(mfrow = c(1, 3))
col <- rainbow(3, alpha = 0.75)
for (name in c("STPCA", "TPCA_SO", "TPCA_SI")) {
  samp <- get(paste0("samp_", name))
  scores <- get(paste0("scores_", name))
  h <- get(paste0("h_", name))
  plot(ks::kde(samp, h = h, gridsize = 500),
    xlim = c(-pi, pi), axes = FALSE,
    xlab = "First scores", ylim = c(0, 0.55)
  )
  sdetorus::torusAxis(1)
  axis(2)
  abline(v = c(-pi, pi), lty = 2)
  for (i in unique(scores$membership)) {
    rug(scores$X[scores$membership == i], col = col[i], ticksize = 0.03)
  }
}

# Speedup hack for kernel mean shift clustering
x_spline <- seq(-1.35 * pi, 1.35 * pi, by = 0.005)
y_STPCA <- sapply(x_spline, function(xx) {
  eta_hat(xx, data = samp_STPCA, h = h_STPCA)
})
y_TPCA_SO <- sapply(x_spline, function(xx) {
  eta_hat(xx, data = samp_TPCA_SO, h = h_TPCA_SO)
})
y_TPCA_SI <- sapply(x_spline, function(xx) {
  eta_hat(xx, data = samp_TPCA_SI, h = h_TPCA_SI)
})
eta_hat_spline_STPCA <- splinefun(x = x_spline, y = y_STPCA)
eta_hat_spline_TPCA_SO <- splinefun(x = x_spline, y = y_TPCA_SO)
eta_hat_spline_TPCA_SI <- splinefun(x = x_spline, y = y_TPCA_SI)

# Run kernel mean shift clustering
x_kms <- seq(-1.25 * pi, 1.25 * pi, by = 0.001)
kms_STPCA <- kms_TPCA_SO <- kms_TPCA_SI <- numeric(length(x_kms))
N <- length(x_kms)
pb <- txtProgressBar(style = 3)
for (i in 1:N) {
  kms_STPCA[i] <- euler(x = x_kms[i], eta_hat_efic = eta_hat_spline_STPCA)
  kms_TPCA_SO[i] <- euler(x = x_kms[i], eta_hat_efic = eta_hat_spline_TPCA_SO)
  kms_TPCA_SI[i] <- euler(x = x_kms[i], eta_hat_efic = eta_hat_spline_TPCA_SI)
  setTxtProgressBar(pb, i / N)
}
kms_STPCA <- sdetorus::toPiInt(kms_STPCA)
kms_TPCA_SO <- sdetorus::toPiInt(kms_TPCA_SO)
kms_TPCA_SI <- sdetorus::toPiInt(kms_TPCA_SI)

# Obtain modes for data, and cluster them
unique_modes_STPCA <- sort(unique(round(kms_STPCA, 3)), decreasing = TRUE)
arrival_modes_STPCA <- sapply(scores_STPCA$X, euler,
  eta_hat_efic = eta_hat_spline_STPCA
)
labels_modes_STPCA <- match(
  x = round(arrival_modes_STPCA, 3),
  table = unique_modes_STPCA
)

unique_modes_TPCA_SO <- sort(unique(round(kms_TPCA_SO, 3)), decreasing = TRUE)
arrival_modes_TPCA_SO <- sapply(scores_TPCA_SO$X, euler,
  eta_hat_efic = eta_hat_spline_TPCA_SO
)
labels_modes_TPCA_SO <- match(
  x = round(arrival_modes_TPCA_SO, 3),
  table = unique_modes_TPCA_SO
)
unique_modes_TPCA_SI <- sort(unique(round(kms_TPCA_SI, 3)), decreasing = TRUE)
unique_modes_TPCA_SI <- unique_modes_TPCA_SI[c(3, 1, 2, 5, 4)]
arrival_modes_TPCA_SI <- sapply(scores_TPCA_SI$X, euler,
  eta_hat_efic = eta_hat_spline_TPCA_SI
)
labels_modes_TPCA_SI <- match(
  x = round(arrival_modes_TPCA_SI, 3),
  table = unique_modes_TPCA_SI
)

# Cluster kms
labels_kms_STPCA <- match(x = round(kms_STPCA, 3), table = unique_modes_STPCA)
labels_kms_TPCA_SO <- match(
  x = round(kms_TPCA_SO, 3),
  table = unique_modes_TPCA_SO
)
labels_kms_TPCA_SI <- match(
  x = round(kms_TPCA_SI, 3),
  table = unique_modes_TPCA_SI
)

# Check matching with $membership
tab_STPCA <- table(
  "assigned" = labels_modes_STPCA,
  "true" = scores_STPCA$membership
)
tab_TPCA_SO <- table(
  "assigned" = labels_modes_TPCA_SO,
  "true" = scores_TPCA_SO$membership
)
tab_TPCA_SI <- table(
  "assigned" = labels_modes_TPCA_SI,
  "true" = scores_TPCA_SI$membership
)
tab_STPCA
sum(tab_STPCA) - sum(diag(tab_STPCA))
sum(diag(tab_STPCA)) / sum(tab_STPCA)
tab_TPCA_SO
sum(tab_TPCA_SO) - sum(diag(tab_TPCA_SO))
sum(diag(tab_TPCA_SO)) / sum(tab_TPCA_SO)
tab_TPCA_SI
sum(tab_TPCA_SI) - sum(diag(tab_TPCA_SI))
sum(diag(tab_TPCA_SI)) / sum(tab_TPCA_SI)

plot_scores <- function(name) {

  # Obtain relevant points
  scores1 <- get(paste0("scores_", name))
  scores2 <- switch(name,
    "STPCA" = scores1[, 2],
    "TPCA_SI" = scores1[, 6],
    "TPCA_SO" = scores1[, 6]
  )
  membership <- scores1$membership
  scores1 <- scores1$X
  labels_kms <- get(paste0("labels_kms_", name))
  labels_rle <- rle(labels_kms)
  positions_antimodes <- c(1, cumsum(labels_rle$lengths))

  plot(scores1, scores2,
    pch = 21, xlab = "",
    ylab = "", bg = col[membership], axes = FALSE,
    xlim = c(-pi, pi), ylim = c(-pi / 2, pi / 2)
  )
  sdetorus::torusAxis()
  for (i in unique(membership)) {
    rug(
      x = scores1[membership == i], side = 1,
      col = col[i], ticksize = 0.03
    )
    rug(
      x = scores2[membership == i], side = 2,
      col = col[i], ticksize = 0.03
    )
  }

  # Background info
  for (k in seq_along(labels_rle$values)) {
    begin <- positions_antimodes[k]
    end <- positions_antimodes[k + 1]
    polygon(
      x = c(x_kms[begin], x_kms[begin:end], x_kms[end]),
      y = c(-1.25 * pi, rep(1.25 * pi, length(begin:end)), -1.25 * pi),
      col = rainbow(5, alpha = 0.15)[labels_rle$values[k]],
      border = NA
    )
  }
  abline(v = x_kms[positions_antimodes], lty = 3)
}

par(mfcol = c(2, 3), mai = c(.2, .2, .2, .2))
col <- rainbow(5)
for (name in c("TPCA_SI", "TPCA_SO", "STPCA")) {

  # Plot scores
  plot_scores(name)

  # Obtain relevant points
  samp <- get(paste0("samp_", name))
  scores <- get(paste0("scores_", name))
  h <- get(paste0("h_", name))
  kms <- get(paste0("kms_", name))
  labels_kms <- get(paste0("labels_kms_", name))
  unique_modes <- get(paste0("unique_modes_", name))
  labels_rle <- rle(labels_kms)
  positions_antimodes <- c(1, cumsum(labels_rle$lengths))

  # Plot curve and axes
  plot(ks::kde(samp, h = h, gridsize = 1e3),
    xlim = c(-pi, pi),
    axes = FALSE, xlab = "", ylab = "",
    ylim = c(0, ifelse(name == "TPCA_SI", 0.55, 0.3))
  )
  sdetorus::torusAxis(1)
  axis(2)
  abline(h = 0, lty = 3)

  # Plot modes
  for (i in seq_along(unique_modes)) {
    segments(
      x0 = unique_modes[i], y0 = 0, x1 = unique_modes[i],
      y1 = ks::kde(
        x = samp, h = h, eval.points = unique_modes[i],
        binned = FALSE
      )$estimate,
      col = rainbow(5)[i], lwd = 2
    )
  }

  # Plot domains of attraction
  kde <- ks::kde(x = samp, h = h, eval.points = x_kms, binned = FALSE)
  for (k in seq_along(labels_rle$values)) {
    begin <- positions_antimodes[k]
    end <- positions_antimodes[k + 1]
    polygon(
      x = c(
        kde$eval.points[begin],
        kde$eval.points[begin:end],
        kde$eval.points[end]
      ),
      y = c(0, kde$estimate[begin:end], 0),
      col = rainbow(5, alpha = 0.15)[labels_rle$values[k]],
      border = NA
    )
  }
  lines(kde$eval.points, kde$estimate)
  abline(v = x_kms[positions_antimodes], lty = 3)

  # Plot rug
  for (i in unique(scores$membership)) {
    rug(scores$X[scores$membership == i], col = col[i], ticksize = 0.03)
  }
}
