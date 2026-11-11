## Number of landmarks: full-sample stress and runtime (replaces the draft's
## table for circular bands on S^2 x S^2 with n = 400)
devtools::load_all("sppca",
                   quiet = TRUE)
fig_dir <- "analysis/figures"
out_dir <- "analysis/results"

ells <- c(10, 20, 40, 50, 100, 200, 400)
M <- 20
res <- NULL
for (m in seq_len(M)) {
  set.seed(m)
  X <- r_circles_polysph(n = 400, r = 2, sigma = 0.05)
  attr(X, "theta") <- NULL
  for (ell in ells) {
    t0 <- proc.time()[[3]]
    f <- smds_polysph(X, d = c(2, 2), p = 4, n_landmarks = ell,
                      landmarks = if (ell == 400) "all" else "ward")
    tt <- proc.time()[[3]] - t0
    res <- rbind(res, data.frame(rep = m, ell = ell, time = tt,
                                 stress = stress_polysph(X, f$Y, d = c(2, 2)),
                                 rho = f$rho))
  }
}
saveRDS(res, file.path(out_dir, "landmark_study.rds"))
agg <- do.call(rbind, lapply(split(res, res$ell), function(z)
  data.frame(ell = z$ell[1], time = mean(z$time), stress = mean(z$stress),
             sd_stress = sd(z$stress), rho = mean(z$rho))))
print(agg, digits = 3)

## LaTeX table
bs <- function(x) gsub("@", "\\\\", x)
L <- c("@begin{table}[t]", "@centering@small", "@begin{tabular}{rrrr}",
       "@toprule",
       "Landmarks $@ell$ & Time (s) & Mean stress & SD of stress@@",
       "@midrule")
for (i in seq_len(nrow(agg))) {
  L <- c(L, sprintf("%d%s & %.2f & %s & %s@@", agg$ell[i],
                    if (agg$ell[i] == 400) " (full)" else "", agg$time[i],
                    formatC(agg$stress[i], format = "e", digits = 2),
                    formatC(agg$sd_stress[i], format = "e", digits = 1)))
}
L <- c(L, "@bottomrule", "@end{tabular}",
       paste0("@caption{Scaled MDS of $n=400$ points along circular bands on ",
              "$@Sp^2@times@Sp^2$ ($p=4$): mean running time and normalized ",
              "full-sample stress over $20$ replicates, as a function of the ",
              "number of Ward landmarks. The last row is the full scaled MDS.}"),
       "@label{tab:landmarks}", "@end{table}")
writeLines(bs(L), "analysis/tab_landmarks.tex")
