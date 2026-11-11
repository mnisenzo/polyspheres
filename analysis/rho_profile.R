## Stress profile in the radius: concentrated (hippocampus) vs spread data
devtools::load_all("sppca",
                   quiet = TRUE)
fig_dir <- "analysis/figures"
out_dir <- "analysis/results"

data("hippocampus", package = "polykde")
Xh <- do.call(cbind, lapply(1:168, function(j) hippocampus$dirs[, j, ]))
set.seed(1)
Xc <- r_circles_polysph(n = 177, r = 10, sigma = 0.05)
attr(Xc, "theta") <- NULL
cases <- list(hippocampus = list(X = Xh, d = rep(2, 168), p = 5),
              circles = list(X = Xc, d = rep(2, 10), p = 5))
prof <- lapply(cases, function(cs) {
  f0 <- smds_polysph(cs$X, d = cs$d, p = cs$p, landmarks = "all")
  base <- if (f0$rho > 50) sqrt(length(cs$d)) else f0$rho
  rhos <- base * 2^seq(-2, 6, by = 0.5)
  st <- sapply(rhos, function(rho)
    smds_polysph(cs$X, d = cs$d, p = cs$p, rho = rho, landmarks = "all")$stress[1])
  list(rho = rhos, stress = st, rho_hat = f0$rho, stress_hat = f0$stress[1],
       sqrt_r = sqrt(length(cs$d)))
})
saveRDS(prof, file.path(out_dir, "rho_profile.rds"))

pdf(file.path(fig_dir, "rho_profile.pdf"), width = 6.4, height = 2.6)
par(mfrow = c(1, 2), mar = c(3, 3.2, 1.4, 0.5), mgp = c(1.9, 0.5, 0),
    cex = 0.75)
for (nm in names(prof)) {
  pr <- prof[[nm]]
  plot(pr$rho / pr$sqrt_r, pr$stress, log = "xy", type = "b", pch = 16,
       xlab = expression(rho / sqrt(r)), ylab = "normalized stress",
       main = nm, cex.main = 0.95)
  abline(v = pr$rho_hat / pr$sqrt_r, lty = 2, col = "#E66101")
  abline(v = 1, lty = 3, col = "grey40")
}
dev.off()
print(lapply(prof, function(p) round(c(rho_hat = p$rho_hat, sqrt_r = p$sqrt_r,
                                         min = min(p$stress), max = max(p$stress)), 4)))
