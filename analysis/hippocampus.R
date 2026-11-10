## Hippocampus s-reps (polykde::hippocampus): SP-PCA vs competitors
devtools::load_all("sppca",
                   quiet = TRUE)
fig_dir <- "analysis/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
out_dir <- "analysis/results"
dir.create(out_dir, showWarnings = FALSE)

data("hippocampus", package = "polykde")
X <- do.call(cbind, lapply(1:168, function(j) hippocampus$dirs[, j, ]))
d <- rep(2, 168)
lab <- hippocampus$ids_labs
col_lab <- ifelse(lab, "#D7191C", "#2C7BB6")  # autism: red; typical: blue

## 1. Stress versus target dimension (cached)
ps <- c(1:10, 12, 15, 20, 30)
st_file <- file.path(out_dir, "hippo_stress.rds")
if (file.exists(st_file)) st <- readRDS(st_file) else st <- t(sapply(ps, function(p) {
  f <- smds_polysph(X, d = d, p = p, landmarks = "all")
  c(p = p, stress = f$stress[[1]], rho = f$rho)
}))
saveRDS(st, file.path(out_dir, "hippo_stress.rds"))
print(st)

## 2. Fits ---------------------------------------------------------------------
p <- 5
set.seed(1)
fit <- sppca(X, d = d, p = p, landmarks = "all")
print(fit)
tp <- tangent_pca(X, d)
cp <- cpns(X, d)
saveRDS(list(fit = fit, tp = tp, cp = cp), file.path(out_dir, "hippo_fits.rds"))

## 3. Leave-one-out explained variance with a common inversion ---------------
res <- rbind("SP-PCA" = residual_scores(fit$scores, X, d, Q = 5),
             "Tangent PCA" = residual_scores(tp$scores, X, d, Q = 5),
             "CPNS" = residual_scores(cp$scores, X, d, Q = 5))
native <- rbind(
  "SP-PCA" = residual_sppca(fit, Q = 5),
  "Tangent PCA" = sapply(0:5, function(q) sum(dist_polysph_pair(
    reconstruct_tangent_pca(tp, q, d), X, d, squared = TRUE))),
  "CPNS" = sapply(0:5, function(q) sum(dist_polysph_pair(
    reconstruct_cpns(cp, q, d), X, d, squared = TRUE))))
saveRDS(list(loo = res, native = native), file.path(out_dir, "hippo_resid.rds"))
print(round(1 - res / res[, 1], 3))
print(round(native, 2))

## 4. Class separation ---------------------------------------------------------
auc <- function(s, y) {
  r <- rank(s); n1 <- sum(y); n0 <- sum(!y)
  a <- (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  max(a, 1 - a)
}
sep <- sapply(list("SP-PCA" = fit$scores, "Tangent PCA" = tp$scores,
                   "CPNS" = cp$scores), function(S) {
  sapply(1:3, function(j) {
    a <- auc(S[, j], lab)
    set.seed(j)
    perm <- replicate(2000, auc(S[, j], sample(lab)))
    c(auc = a, p = mean(perm >= a))
  })
})
print(round(sep, 3))
saveRDS(sep, file.path(out_dir, "hippo_sep.rds"))

## 5. Figures ------------------------------------------------------------------
# (a) stress vs p
pdf(file.path(fig_dir, "hippo_stress.pdf"), width = 3.4, height = 2.6)
par(mar = c(3, 3.2, 0.8, 0.5), mgp = c(1.9, 0.5, 0), cex = 0.8)
plot(st[, "p"], st[, "stress"], log = "y", type = "b", pch = 16,
     xlab = expression("target dimension " * p),
     ylab = "normalized stress")
abline(v = p, lty = 2, col = "grey50")
dev.off()

# (b) scores matrices, OODA style
for (m in c("sppca", "tangent")) {
  S <- if (m == "sppca") fit$scores else tp$scores
  v <- if (m == "sppca") 100 * (1 - res[1, 2:5] / res[1, 1]) -
    c(0, 100 * (1 - res[1, 2:4] / res[1, 1])) else
    100 * (1 - res[2, 2:5] / res[2, 1]) - c(0, 100 * (1 - res[2, 2:4] / res[2, 1]))
  colnames(S) <- paste0("PC", seq_len(ncol(S)))
  pdf(file.path(fig_dir, paste0("hippo_scores_", m, ".pdf")), width = 5.5,
      height = 5.5)
  set.seed(2)
  plot_scores_matrix(S, n_comp = 4, col = col_lab, var = v, cex = 0.55)
  dev.off()
}

# (c) residual curves
pdf(file.path(fig_dir, "hippo_resid.pdf"), width = 3.6, height = 2.8)
par(mar = c(3, 3.2, 0.8, 0.5), mgp = c(1.9, 0.5, 0), cex = 0.8)
cols <- c("SP-PCA" = "#E66101", "Tangent PCA" = "#5E3C99", "CPNS" = "#1B9E77")
matplot(0:5, t(100 * (1 - res / res[, 1])), type = "b", pch = c(16, 17, 15),
        lty = 1, col = cols, xlab = expression("rank " * q),
        ylab = "% variance explained (LOO)", ylim = c(0, 65))
legend("topleft", legend = names(cols), col = cols, pch = c(16, 17, 15),
       lty = 1, bty = "n")
dev.off()

# (d) first mode of variation as hippocampus shapes
base_m <- apply(hippocampus$base, c(2, 3), mean)
rad_m <- colMeans(hippocampus$rads)
sd1 <- sd(fit$scores[, 1])
tt <- c(-2, -1, 0, 1, 2) * sd1
modes <- list("SP-PCA" = mode_sppca(fit, j = 1, t = tt, method = "regression"),
              "Tangent PCA" = reconstruct_tangent_pca(
                list(mu = tp$mu, vbar = tp$vbar, rot = tp$rot,
                     scores = cbind(c(-2, -1, 0, 1, 2) * sd(tp$scores[, 1]))),
                q = 1, d = d))
shape_of <- function(dirs_row) {
  D <- matrix(dirs_row, ncol = 3, byrow = TRUE)
  base_m + rad_m * D
}
pal <- colorRampPalette(c("#2166AC", "#67A9CF", "grey30", "#EF8A62", "#B2182B"))(5)
R <- qr.Q(qr(cbind(c(1, 0.3, 0.2), c(0, 1, 0), c(0, 0, 1))))
proj2 <- function(P, R) P %*% R[, 1:2]
pdf(file.path(fig_dir, "hippo_mode1.pdf"), width = 7, height = 4.2)
par(mfrow = c(2, 2), mar = c(0.5, 0.5, 1.5, 0.5))
pcs <- prcomp(base_m)$rotation
views <- list(pcs[, c(1, 2)], pcs[, c(1, 3)])
for (m in names(modes)) {
  for (v in seq_along(views)) {
    V <- views[[v]]
    allp <- do.call(rbind, lapply(1:5, function(i) shape_of(modes[[m]][i, ]) %*% V))
    plot(allp, type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "",
         main = if (v == 1) paste0(m, ": first mode, view 1") else "view 2",
         cex.main = 0.85)
    B <- base_m %*% V
    for (i in c(1, 5, 3)) {
      P <- shape_of(modes[[m]][i, ]) %*% V
      segments(B[, 1], B[, 2], P[, 1], P[, 2],
               col = adjustcolor(pal[i], 0.35), lwd = 0.6)
      points(P, pch = 16, cex = 0.45, col = pal[i])
    }
    points(B, pch = 16, cex = 0.25, col = "grey60")
    if (v == 1 && m == "SP-PCA") {
      legend("bottomleft", legend = c(expression(-2 * sd), "mean",
                                       expression(+2 * sd)),
             col = pal[c(1, 3, 5)], pch = 16, bty = "n", cex = 0.75)
    }
  }
}
dev.off()

# (e) spoke-wise magnitude of the first mode on the mean skeleton
disp <- sapply(1:168, function(j) {
  ind <- 3 * (j - 1) + 1:3
  a <- modes[["SP-PCA"]][1, ind]; b <- modes[["SP-PCA"]][5, ind]
  acos(pmin(1, sum(a * b)))
})
saveRDS(disp, file.path(out_dir, "hippo_disp.rds"))
pdf(file.path(fig_dir, "hippo_disp.pdf"), width = 6, height = 2.6)
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 1.2, 0.5))
heat <- colorRampPalette(c("#FFFFCC", "#FD8D3C", "#800026"))(100)
ci <- cut(disp, 100, labels = FALSE)
for (v in seq_along(views)) {
  V <- views[[v]]
  P <- shape_of(modes[["SP-PCA"]][3, ]) %*% V
  B <- base_m %*% V
  plot(P, type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "",
       main = if (v == 1) "angular change of spokes along SP-PCA mode 1" else "",
       cex.main = 0.8)
  segments(B[, 1], B[, 2], P[, 1], P[, 2], col = heat[ci], lwd = 1.2)
  points(P, pch = 16, cex = 0.6, col = heat[ci])
}
dev.off()
cat("max angular change (deg):", round(max(disp) * 180 / pi, 1), "\n")
