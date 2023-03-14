### Sunspots Data
### Section 5.1 of the main paper
### Warning: the analysis of sunspots will require multiple hours to run

source("C:/Users/minise/Downloads/STPCA/Script1_Setup_and_Functions.R")

## Loading data
set.seed(10)
data("sunspots_births")
solar_cycle23 <- sunspots_births$theta[sunspots_births$cycle == 23]
n <- length(solar_cycle23)
X23 <- cbind(
  solar_cycle23[-c(n - 1, n)], # Theta_{i}
  solar_cycle23[-c(1, n)], # Theta_{i+1}
  solar_cycle23[-c(1, 2)]
) # Theta_{i+2}
set.seed(1)
Y23 <- X23[sample(1:5371, size = 1000), ]
Y23 <- (Y23 + pi) %% (2 * pi) - pi
n <- 1002
Z23 <- Y23 * 180 / pi
X23 <- (X23 + pi) %% (2 * pi) - pi

## ST-PCA (Steps 1 and 2) for the full sunspots and the sample of 1000 points
sunspots_full <- MDS_PNS(x = X23)
sunspots_1000 <- MDS_PNS(x = Y23)
plot(sunspots_full$pns_data[, 1:2]) ## Figure 1a
## Prediction (step 3) and visualization of data with 1st toroidal component
## used 125 points for predictions
d_sphere <- 3
sc <- march_PNS1_shapes(seq(0, sunspots_1000$pns_info$radii[d_sphere] * 2 * pi,
  l = 125
), PNS = sunspots_1000$pns_info)
pred_sun_full <- prediction(
  y_new = sc, x_old = X23,
  y_old = sunspots_full$MDS_data,
  r_star = sunspots_full$MDS_info$rstar,
  d_sphere = d_sphere
)
all_torus <- rbind(X23, pred_sun_full)
all_sphere <- rbind(sunspots_full$MDS_data, sc)
all_torus <- all_torus * 180 / pi ## Radians --> angles for Python
# Plotting ST-PCA's first component
plot3d(all_torus,
  col = c(rep("#000000", 5371), rep("blue", 125)),
  xlim = c(min(all_torus[, 1]), max(all_torus[, 1])),
  ylim = c(-180, 180), zlim = c(-180, 180),
  pch = c(rep(".", 5371), rep("*", 125))
)
## R^2 for Linear fit of predicted points
## First take the representation of the points on the edges of the cube
## to the straight line
prediction_full_dataframe <-
  as.data.frame(sdetorus::unwrapCircSeries(pred_sun_full))
fit1 <- lm(V3 ~ V1 + V2, data = prediction_full_dataframe)
summary(fit1)$adj.r.squared

## Percent of variance explained for the sample of 1000 points
## This calculations is computationally expensive
sunspots_variance_1000 <- percent_var_explained(MDSPNS = sunspots_1000, Y23)

## Histograms
theta_23 <- subset(sunspots_births, select = theta, cycle == 23)[, 1]
scores_STPCA <- as.data.frame(sunspots_full$pns_data)
scores_TPCA_SI <- read.csv("Sunspotsfull_TPCA_SI.csv",
  sep = ",",
  header = FALSE
)
scores_TPCA_SO <- read.csv("Sunspotsfull_TPCA_SO.csv",
  sep = ",",
  header = FALSE
)
# Homogenize data: transform to (-pi, pi) flip signs to have better
# comparisons, and rename variables
scores_TPCA_SO[, 1:3] <- -scores_TPCA_SO[, 1:3] / 180 * pi
scores_TPCA_SI[, 1:3] <- -scores_TPCA_SI[, 1:3] / 180 * pi
names(scores_TPCA_SO)[3] <- names(scores_TPCA_SI)[3] <- "X"
scores_STPCA$X <- NULL
names(scores_STPCA)[1] <- "X"
# Histograms
par(mfrow = c(1, 3))
col <- rainbow(3, alpha = 0.75)
for (name in c("STPCA", "TPCA_SO", "TPCA_SI")) {
  scores <- get(paste0("scores_", name))
  hist(scores$X,
    breaks = seq(-pi, pi, l = 15), main = name, prob = TRUE,
    ylim = c(0, 0.5), xlab = "First scores", axes = FALSE
  )
  sdetorus::torusAxis(sides = 1)
  axis(2)
  rug(scores$X, col = gray(1, alpha = 0.5))
  abline(h = 1 / (2 * pi), lty = 2, col = 2)
}

## Uniformity tests
type <- c("Gine_Fn", "Watson", "PAD")
set.seed(1)
sphunif::unif_test(data = theta_23, p_value = "asymp", type = type)
sphunif::unif_test(data = scores_STPCA$X, p_value = "asymp", type = type)
sphunif::unif_test(data = scores_TPCA_SO$X, p_value = "asymp", type = type)
sphunif::unif_test(data = scores_TPCA_SI$X, p_value = "asymp", type = type)

## Standard Euclidean PCA in ambient space for visualization

epca <- princomp(x = X23)
summary(epca)
plot(epca$scores[, 1:2]) ## Figure 1b

## Standard Euclidean PCA in ambient space for percent of variance explained
## by the first components, on the subset of 1000 points

epca_subset <- princomp(x = Y23)
summary(epca_subset)

## Dihedral PCA (dPCA) for the subset of 1000 points

dpcadata <- matrix(NA, nrow = nrow(Y23), ncol = 2 * ncol(Y23))

for (i in 1:ncol(Y23)) {
  dpcadata[, 1] <- sapply(Y23[, 1], cos)
  dpcadata[, 2] <- sapply(Y23[, 1], sin)
  dpcadata[, 3] <- sapply(Y23[, 2], cos)
  dpcadata[, 4] <- sapply(Y23[, 2], sin)
  dpcadata[, 5] <- sapply(Y23[, 3], cos)
  dpcadata[, 6] <- sapply(Y23[, 3], sin)
}
dpca <- princomp(x = dpcadata)
summary(dpca)

## Complex dPCA (for subset of 1000 points)
## For formulas and details, see Altis et al.

z <- matrix(NA, nrow = nrow(Y23), ncol = 3)
for (i in 1:nrow(Y23)) {
  for (j in 1:3) {
    z[i, j] <- exp(sqrt(as.complex(-1)) * as.complex(Y23[i, j]))
  }
}

C <- matrix(NA, ncol = 3, nrow = 3)
for (i in 1:3) {
  for (j in 1:3) {
    C[i, j] <- (z[, i] - mean(z[, i])) %*% (Conj(z[, j]) - Conj(mean(z[, j])))
  }
}

eigenvalues <- eigen(C)
variance.explained <- max(eigenvalues$values) / sum(eigenvalues$values)
eigenvalues$values[1] / sum(eigenvalues$values)
eigenvalues$values[2] / sum(eigenvalues$values)
eigenvalues$values[3] / sum(eigenvalues$values)

## Angular PCA of Riccardi et al for subset of 1000 points

S <- rep(NA, ncol(Y23))
C <- rep(NA, ncol(Y23))
M <- rep(NA, ncol(Y23))
apcadata <- matrix(NA, nrow = nrow(Y23), ncol = ncol(Y23))
Z23 <- Z23 + 180 # use the degrees representation
for (j in 1:ncol(Z23)) {
  S[j] <- sum(sin(Z23[, j]))
  C[j] <- sum(cos(Z23[, j]))
  if (C[j] > 0) {
    delta <- 0
  } else if (C[j] < 0) {
    delta <- pi
  } else {
    delta <- 2 * pi
  }
  ## This delta is the correct directional mean for Circular data
  ## Paper by Riccardi et al. has a typo in the delta
  ## One can find the correct formula in the citation provided by Riccardi et al
  ## Page 31 of Fisher's "Statistical Analysis of Circular Data"
  M[j] <- atan2(S[j], C[j]) + delta
  apcadata[, j] <- Y23[, j] - M[j]
  for (i in 1:nrow(Y23)) {
    if (apcadata[i, j] > pi) {
      apcadata[i, j] <- apcadata[i, j] - 2 * pi
    } else if (apcadata[i, j] < -pi) {
      apcadata[i, j] <- apcadata[i, j] + 2 * pi
    } else {
      apcadata[i, j] <- apcadata[i, j]
    }
  }
  ## There is also a typo in the transformed angles, phi, in Riccardi et al.
  ## In the given transformation, the final angles are not in [-pi,pi) as described
  ## We have implemented the transformation that is described in words
  ## According to which this data will be in the required interval
}
apcadata <- ((apcadata - 180) / 180) * pi
apca <- princomp(x = apcadata)
summary(apca)

## dPCA+ of Sitter et al

dpca_plus <- matrix(NA, nrow = nrow(Y23), ncol = ncol(Y23))

# identify lowest density in 1st dimension
histog1 <- hist(Y23[, 1], breaks = seq(from = -pi, to = pi, by = 2 * pi / 10))
translation1 <- (histog1$breaks[which.min(histog1$counts)] +
  histog1$breaks[which.min(histog1$counts) + 1]) / 2
# Bin of lowest density is [1/5* pi, 3/5 *pi], center is 2/5 *pi, translate by -7/5*pi
dpca_plus[, 1] <- Y23[, 1] - translation1

# identify lowest density in 2nd dimension
histog2 <- hist(Y23[, 2], breaks = seq(from = -pi, to = pi, by = 2 * pi / 10))
translation2 <- (histog2$breaks[which.min(histog2$counts)] +
  histog2$breaks[which.min(histog2$counts) + 1]) / 2
# Bin of lowest density is [1/5* pi, 3/5 *pi], center is 2/5 *pi, translate by -7/5*pi
dpca_plus[, 2] <- Y23[, 2] - translation2

## identify lowest density in 3rd dimension
histog3 <- hist(Y23[, 3], breaks = seq(from = -pi, to = pi, by = 2 * pi / 10))
translation3 <- (histog3$breaks[which.min(histog3$counts)] +
  histog3$breaks[which.min(histog3$counts) + 1]) / 2
# Bin of lowest density is [1/5* pi, 3/5 *pi], center is 2/5 *pi, translate by -7/5*pi
dpca_plus[, 3] <- Y23[, 3] - translation3

dpca_plus <- (dpca_plus + pi) %% (2 * pi) - pi ## imposing periodicity
dpca_plus_pca <- princomp(x = dpca_plus)

summary(dpca_plus_pca)
