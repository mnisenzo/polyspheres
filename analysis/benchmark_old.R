## Old (GitHub mnisenzo/polyspheres) vs new (sppca) implementation: time & stress
suppressMessages({
  library(polykde); library(sphunif); library(shapes); library(rotasym)
})
old_dir <- "R"
old <- new.env(parent = asNamespace("polykde"))
for (f in c("auxiliary.R", "angles.R", "clust_centers.R", "efficient_mds.R")) {
  src <- readLines(file.path(old_dir, f))
  src <- src[!grepl("^\\s*library\\(", src)]
  eval(parse(text = src), envir = old)
}
devtools::load_all("sppca",
                   quiet = TRUE)

res <- NULL
for (n in c(200, 400, 800, 1600)) {
  for (rep in 1:3) {
    set.seed(rep * 100 + n)
    X <- r_circles_polysph(n = n, r = 4, sigma = 0.1)
    attr(X, "theta") <- NULL
    d <- rep(2, 4)
    ell <- round(0.1 * n)
    # Old: landmark SMDS + per-point optim + PNS (seq.test)
    t_old <- NA; s_old <- NA
    if (n <= 800) {
      t0 <- proc.time()[[3]]
      o <- tryCatch(old$mds_pns(x = X, angular = FALSE, dims = d,
                                landmark_size = 0.1, gradient = FALSE,
                                pc = 1),
                    error = function(e) {message(conditionMessage(e)); NULL})
      t_old <- proc.time()[[3]] - t0
      if (!is.null(o)) s_old <- stress_polysph(X, o$mds_data, d = d)
    }
    # New: same number of landmarks, target S^8 (= sum(d))
    t0 <- proc.time()[[3]]
    f <- sppca(X, d = d, p = 8, n_landmarks = ell)
    t_new <- proc.time()[[3]] - t0
    s_new <- stress_polysph(X, f$smds$Y, d = d)
    t_new_smds <- sum(f$smds$time)
    # New, full SMDS (all points as landmarks) for reference
    t_full <- NA; s_full <- NA
    if (n <= 800) {
      t0 <- proc.time()[[3]]
      g <- smds_polysph(X, d = d, p = 8, landmarks = "all")
      t_full <- proc.time()[[3]] - t0
      s_full <- stress_polysph(X, g$Y, d = d)
    }
    res <- rbind(res, data.frame(n = n, rep = rep, t_old = t_old,
                                 t_new = t_new, t_new_smds = t_new_smds,
                                 t_full = t_full, stress_old = s_old,
                                 stress_new = s_new, stress_full = s_full))
    print(tail(res, 1), digits = 4)
  }
}
saveRDS(res, "analysis/results/benchmark.rds")
print(aggregate(. ~ n, data = res[, -2], FUN = mean, na.action = na.pass),
      digits = 4)
