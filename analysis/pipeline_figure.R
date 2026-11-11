## Overview of the SP-PCA pipeline on S^2 x S^2 (four panels)
devtools::load_all("sppca",
                   quiet = TRUE)
fig_dir <- "analysis/figures"

set.seed(11)
n <- 300
X <- r_circles_polysph(n, r = 2, sigma = 0.06, height = 0.5,
                       theta_range = c(0, 1.6 * pi))
th <- attr(X, "theta")
attr(X, "theta") <- NULL
d <- c(2, 2)
col <- hcl.colors(n, "Spectral")[rank(th)]

# Fit with a two-dimensional target sphere so that it can be drawn
fit <- sppca(X, d = d, p = 2, landmarks = "all")
Y <- fit$smds$Y
pn <- fit$pns$PNS
circle <- t(shapes::PNSe2s(rbind(seq(-pi, pi, l = 300), 0), pn))
tt <- seq(min(fit$scores[, 1]), max(fit$scores[, 1]), l = 200)
mode1 <- mode_sppca(fit, j = 1, t = tt, method = "regression")

view_for <- function(P) {
  # Look at the data from the side that shows the fitted circle best
  v <- colMeans(P)
  if (sqrt(sum(v^2)) < 0.2) v <- svd(scale(P, scale = FALSE))$v[, 3]
  v / sqrt(sum(v^2))
}
pdf(file.path(fig_dir, "pipeline.pdf"), width = 7.2, height = 2.1)
par(mfrow = c(1, 4), mar = c(0.2, 0.2, 1.6, 0.2), font.main = 1)
plot_sphere(X[, 1:3], col = col, view = view_for(X[, 1:3]), cex = 0.45,
            main = expression("(a) data, sphere 1"))
plot_sphere(X[, 4:6], col = col, view = view_for(X[, 4:6]), cex = 0.45,
            main = expression("(a) data, sphere 2"))
plot_sphere(Y, col = col, view = view_for(Y), cex = 0.45,
            lines = list(circle), lines_col = "black", lwd = 1.8,
            main = expression("(b) scaled data on " * S[rho]^2 * " + PNS"))
plot_sphere(X[, 1:3], col = adjustcolor(col, 0.35), view = view_for(X[, 1:3]),
            cex = 0.4, lines = list(mode1[, 1:3]), lines_col = "black",
            lwd = 2, main = "(c) mode 1 on sphere 1")
dev.off()
cat("rho =", fit$smds$rho, " var:", round(fit$var_sphere, 1), "\n")

## OODA-style mode-of-variation display: for PC1 and PC2, the score
## distribution (KDE + jitter, colored by score) and the mode on each sphere
## drawn with the same rainbow coloring
pdf(file.path(fig_dir, "ooda_modes.pdf"), width = 7.2, height = 4.4)
par(mfrow = c(2, 3), mar = c(2.6, 0.6, 1.6, 0.6), mgp = c(1.5, 0.4, 0),
    font.main = 1, cex.axis = 0.75, tcl = -0.2)
for (j in 1:2) {
  s <- fit$scores[, j]
  rg <- quantile(s, c(0.02, 0.98))
  tj <- seq(rg[1], rg[2], l = 60)
  mj <- mode_sppca(fit, j = j, t = tj, method = "regression")
  pal <- hcl.colors(60, "Spectral")
  cs <- pal[pmin(60, pmax(1, round(59 * (s - rg[1]) / diff(rg)) + 1))]
  kd <- density(s, bw = "SJ")
  plot(kd$x, kd$y, type = "l", lwd = 1.5, yaxt = "n", ylab = "",
       xlab = paste0("PC", j, " score"),
       ylim = c(-0.35, 1.05) * max(kd$y),
       main = sprintf("PC%d (%.1f%% of scaled variance)", j,
                      fit$var_sphere[j]))
  set.seed(j)
  points(s, -runif(length(s), 0.05, 0.3) * max(kd$y), col = cs, pch = 16,
         cex = 0.5)
  for (sph in 1:2) {
    ind <- 3 * (sph - 1) + 1:3
    segs <- lapply(1:59, function(i) mj[i:(i + 1), ind])
    plot_sphere(X[, ind], col = adjustcolor("grey60", 0.5),
                view = view_for(X[, ind]), cex = 0.3, lines = segs,
                lines_col = pal[1:59], lwd = 3,
                main = paste0("mode ", j, ", sphere ", sph))
  }
}
dev.off()
