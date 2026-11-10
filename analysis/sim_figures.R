## Figures for the simulation study (OODA-style)
devtools::load_all("sppca",
                   quiet = TRUE)
fig_dir <- "analysis/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## 1. Illustration: (S^2)^10, two active spheres + concentrated noise ---------
set.seed(7)
n <- 300
r <- 10
X <- r_circles_polysph(n, r, sigma = 0.05, n_active = 2, kappa = 200)
th <- attr(X, "theta")
attr(X, "theta") <- NULL
d <- rep(2, r)
col <- hcl.colors(n, "Spectral")[rank(th)]   # rainbow along the latent angle
fit <- sppca(X, d = d, p = 6)
tp <- tangent_pca(X, d)
m1 <- mode_sppca(fit, j = 1, t = seq(-pi, pi, l = 200) * fit$smds$rho *
                   fit$pns$PNS$radii[length(fit$pns$PNS$radii)],
                 method = "curve")
tpm <- reconstruct_tangent_pca(list(mu = tp$mu, vbar = tp$vbar, rot = tp$rot,
                                    scores = cbind(seq(min(tp$scores[, 1]),
                                                       max(tp$scores[, 1]),
                                                       l = 200))),
                               q = 1, d = d)
pdf(file.path(fig_dir, "sim_illustration.pdf"), width = 7.2, height = 2.5)
par(mfrow = c(1, 3), mar = c(0.2, 0.2, 1.4, 0.2))
for (j in c(1, 2, 3)) {
  ind <- 3 * (j - 1) + 1:3
  ctr <- colMeans(X[, ind]); ctr <- ctr / sqrt(sum(ctr^2))
  nrm <- svd(scale(X[, ind], scale = FALSE))$v[, 3]
  vw <- if (j <= 2) nrm + 0.6 * ctr else ctr
  plot_sphere(X[, ind], col = col, view = vw, cex = 0.45,
              lines = list(m1[, ind], tpm[, ind]),
              lines_col = c("black", "#5E3C99"), lwd = 1.8,
              main = if (j <= 2) paste0("sphere ", j, " (active)") else
                paste0("sphere ", j, " (vMF noise)"))
}
legend("bottomright", legend = c("SP-PCA mode 1", "tangent PC1"),
       col = c("black", "#5E3C99"), lwd = 1.8, bty = "n", cex = 0.8)
dev.off()

for (m in c("sppca", "tangent")) {
  S <- if (m == "sppca") fit$scores else tp$scores
  colnames(S) <- paste0("PC", 1:ncol(S))
  pdf(file.path(fig_dir, paste0("sim_scores_", m, ".pdf")), width = 4.6,
      height = 4.6)
  set.seed(3)
  plot_scores_matrix(S, n_comp = 3, col = col, cex = 0.5)
  dev.off()
}

## 2. Summary of the simulation study -----------------------------------------
res <- readRDS("analysis/results/simulation_results.rds")
res <- res[!is.na(res$method), ]
cols <- c("SP-PCA" = "#E66101", "Tangent PCA" = "#5E3C99", "CPNS" = "#1B9E77")
designs <- c(great_full = "great circles", small_arc = "small-circle arcs",
             active2 = "2 active + vMF noise", clusters = "3 clusters")
pdf(file.path(fig_dir, "sim_summary.pdf"), width = 7.2, height = 2.6)
par(mfrow = c(1, 4), mar = c(3.2, 2.8, 1.4, 0.3), mgp = c(1.7, 0.45, 0),
    cex.axis = 0.75, tcl = -0.2)
for (ds in names(designs)) {
  sub <- res[res$design == ds, ]
  rs <- sort(unique(sub$r))
  pos <- unlist(lapply(seq_along(rs), function(i) (i - 1) * 4 + 1:3))
  meths <- names(cols)
  vals <- lapply(rs, function(rr) lapply(meths, function(m)
    100 * sub$expl1[sub$r == rr & sub$method == m]))
  vals <- unlist(vals, recursive = FALSE)
  boxplot(vals, at = pos, col = adjustcolor(rep(cols, length(rs)), 0.5),
          border = rep(cols, length(rs)), xaxt = "n", outline = FALSE,
          ylim = c(60, 100), main = designs[ds], cex.main = 0.85,
          ylab = if (ds == "great_full") "% variance explained, rank 1" else "")
  axis(1, at = (seq_along(rs) - 1) * 4 + 2, labels = paste0("r=", rs))
  if (ds == "clusters") {
    ari <- aggregate(ari ~ r + method, data = sub, FUN = mean)
    for (i in seq_along(rs)) for (k in seq_along(meths)) {
      a <- ari$ari[ari$r == rs[i] & ari$method == meths[k]]
      text((i - 1) * 4 + k, 62, sprintf("%.2f", a), cex = 0.5, srt = 90,
           adj = 0, col = cols[k])
    }
  }
}
legend("bottomright", legend = names(cols), fill = adjustcolor(cols, 0.5),
       border = cols, bty = "n", cex = 0.75, inset = c(0, 0.12))
dev.off()

## 3. Runtime table ------------------------------------------------------------
tt <- aggregate(time ~ design + r + method, data = res, FUN = median)
saveRDS(tt, "analysis/results/sim_times.rds")
tab <- aggregate(cbind(expl1, expl2, ari, time) ~ design + r + method,
                 data = res, FUN = mean, na.action = na.pass)
sds <- aggregate(cbind(expl1) ~ design + r + method, data = res, FUN = sd)
tab$expl1_sd <- sds$expl1
saveRDS(tab, "analysis/results/sim_table.rds")
print(tab, digits = 3)
