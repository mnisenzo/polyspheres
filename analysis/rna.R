## RNA torsion angles on T^7 = (S^1)^7 (data from the ST-PCA companion code):
## SP-PCA as scaled torus PCA, against T-PCA (published scores), tangent PCA,
## and CPNS. Clusters are recovered from the first score.
devtools::load_all("sppca",
                   quiet = TRUE)
src <- "st-pca"
fig_dir <- "analysis/figures"
out_dir <- "analysis/results"

# Data and reference clusters (as in the ST-PCA companion script)
rna <- as.matrix(read.csv(file.path(src, "small_rna.csv"), header = FALSE))
cl1 <- rna[, 9] > -60
cl2 <- rna[, 8] > 0 & rna[, 9] < -60
cl3 <- rna[, 8] < 0 & rna[, 9] < -60
ord <- c(which(cl1), which(cl2), which(cl3))
rna <- rna[ord, ]
lab <- rep(1:3, c(sum(cl1), sum(cl2), sum(cl3)))
theta <- rna[, 1:7] * pi / 180
X <- polykde::angles_to_torus(theta)
d <- rep(1, 7)
cat("n =", nrow(X), " clusters:", table(lab), "\n")

# Fits
set.seed(1)
fit <- sppca(X, d = d, p = 7, landmarks = "all")
print(fit)
tp <- tangent_pca(X, d)
cp <- cpns(X, d)
# Published T-PCA scores: already in cluster order; first component is V7
# (degrees), last column is the cluster membership (1, 2, 5)
tpca_si <- read.csv(file.path(src, "tpca_rna_SI.csv"), header = TRUE)
tpca_so <- read.csv(file.path(src, "tpca_rna_SO.csv"), header = TRUE)
stopifnot(all(match(tpca_si[, 8], c(1, 2, 5)) == lab))

# Cluster recovery from the first score: 3-means on the (circular) score,
# best label permutation accuracy and adjusted Rand index
ari <- function(a, b) {
  tab <- table(a, b)
  s <- sum(choose(tab, 2)); s1 <- sum(choose(rowSums(tab), 2))
  s2 <- sum(choose(colSums(tab), 2)); N <- choose(length(a), 2)
  (s - s1 * s2 / N) / ((s1 + s2) / 2 - s1 * s2 / N)
}
perms <- rbind(1:3, c(1, 3, 2), c(2, 1, 3), c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))
acc <- function(a, b) max(apply(perms, 1, function(p) mean(p[a] == b)))
recover <- function(s, circular = FALSE) {
  Z <- if (circular) cbind(cos(s), sin(s)) else cbind(s)
  set.seed(1)
  km <- kmeans(Z, centers = 3, nstart = 50)$cluster
  c(accuracy = acc(km, lab), ARI = ari(km, lab))
}
# The first PNS score is an angle on the last circle (radius 1 in resmat)
s_sp <- fit$scores[, 1] / fit$smds$rho
rec <- rbind(
  "SP-PCA" = recover(s_sp, circular = TRUE),
  "T-PCA (SI)" = recover(tpca_si$V7 * pi / 180, circular = TRUE),
  "T-PCA (SO)" = recover(tpca_so$V7 * pi / 180, circular = TRUE),
  "Tangent PCA" = recover(tp$scores[, 1]),
  "CPNS" = recover(cp$scores[, 1]))
print(round(rec, 3))

# Explained variance with the common leave-one-out inversion
ev <- rbind("SP-PCA" = residual_scores(fit$scores, X, d, Q = 3),
            "Tangent PCA" = residual_scores(tp$scores, X, d, Q = 3),
            "CPNS" = residual_scores(cp$scores, X, d, Q = 3))
print(round(1 - ev / ev[, 1], 3))
saveRDS(list(rec = rec, ev = ev, fit = fit, lab = lab),
        file.path(out_dir, "rna.rds"))

# Figure: first two scores of SP-PCA and tangent PCA, colored by cluster
cols <- c("#1B9E77", "#D95F02", "#7570B3")[lab]
pdf(file.path(fig_dir, "rna_scores.pdf"), width = 7, height = 2.6)
par(mfrow = c(1, 3), mar = c(3, 3, 1.5, 0.5), mgp = c(1.8, 0.5, 0), cex = 0.75)
plot(fit$scores[, 1:2], col = cols, pch = 16, cex = 0.7,
     xlab = "PC1", ylab = "PC2", main = "SP-PCA")
plot(tp$scores[, 1:2], col = cols, pch = 16, cex = 0.7,
     xlab = "PC1", ylab = "PC2", main = "Tangent PCA")
plot(cp$scores[, 1:2], col = cols, pch = 16, cex = 0.7,
     xlab = "PC1", ylab = "PC2", main = "CPNS")
dev.off()
