## Hippocampus figures from cached fits (run hippocampus.R first)
devtools::load_all("sppca",
                   quiet = TRUE)
fig_dir <- "analysis/figures"
out_dir <- "analysis/results"
data("hippocampus", package = "polykde")
X <- do.call(cbind, lapply(1:168, function(j) hippocampus$dirs[, j, ]))
d <- rep(2, 168)
lab <- hippocampus$ids_labs
col_lab <- ifelse(lab, "#D7191C", "#2C7BB6")
fits <- readRDS(file.path(out_dir, "hippo_fits.rds"))
fit <- fits$fit; tp <- fits$tp; cp <- fits$cp
res <- readRDS(file.path(out_dir, "hippo_resid.rds"))$loo
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
mag <- 4
pdf(file.path(fig_dir, "hippo_mode1.pdf"), width = 7.2, height = 4.4)
par(mfrow = c(2, 2), mar = c(0.3, 0.3, 1.4, 0.3))
pcs <- prcomp(base_m)$rotation
views <- list(pcs[, c(1, 2)], pcs[, c(1, 3)])
for (m in names(modes)) {
  M0 <- shape_of(modes[[m]][3, ])
  Mm <- M0 + mag * (shape_of(modes[[m]][1, ]) - M0)
  Mp <- M0 + mag * (shape_of(modes[[m]][5, ]) - M0)
  for (v in seq_along(views)) {
    V <- views[[v]]
    B <- base_m %*% V; P0 <- M0 %*% V; Pm <- Mm %*% V; Pp <- Mp %*% V
    plot(rbind(P0, Pm, Pp, B), type = "n", asp = 1, axes = FALSE, xlab = "",
         ylab = "", main = if (v == 1) paste0(m, ": mode 1 (", mag,
                                              "x magnified)") else "",
         cex.main = 0.85)
    segments(B[, 1], B[, 2], P0[, 1], P0[, 2], col = "grey80", lwd = 0.5)
    points(B, pch = 16, cex = 0.2, col = "grey60")
    arrows(P0[, 1], P0[, 2], Pm[, 1], Pm[, 2], length = 0.03, col = pal[1],
           lwd = 1)
    arrows(P0[, 1], P0[, 2], Pp[, 1], Pp[, 2], length = 0.03, col = pal[5],
           lwd = 1)
    points(P0, pch = 16, cex = 0.35, col = "grey20")
    if (v == 1 && m == "SP-PCA") {
      legend("bottomleft", legend = c(expression(-2 * sd), expression(+2 * sd)),
             col = pal[c(1, 5)], lwd = 1.5, bty = "n", cex = 0.75)
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
