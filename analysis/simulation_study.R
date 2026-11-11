## Simulation study: SP-PCA vs tangent PCA vs CPNS on structured polyspheres
library(parallel)

settings <- expand.grid(
  design = c("great_full", "small_arc", "active2", "clusters"),
  r = c(2, 10, 50), stringsAsFactors = FALSE)
settings <- subset(settings, !(design == "active2" & r == 2))
M <- 20
n <- 200

one_rep <- function(i, design, r) {
  devtools::load_all("sppca", quiet = TRUE)
  set.seed(1000 * i + r)
  X <- switch(design,
    great_full = r_circles_polysph(n, r, sigma = 0.05),
    small_arc = r_circles_polysph(n, r, sigma = 0.05, height = 0.6,
                                  theta_range = c(0, 1.5 * pi)),
    active2 = r_circles_polysph(n, r, sigma = 0.05, n_active = 2, kappa = 200),
    clusters = r_clusters_polysph(n, r, K = 3, kappa = 30))
  cl <- attr(X, "cluster")
  attr(X, "theta") <- attr(X, "cluster") <- NULL
  d <- rep(2, r)
  out <- list()
  t0 <- proc.time()[[3]]
  f <- sppca(X, d = d, p = min(2 * r, 10))
  t_sp <- proc.time()[[3]] - t0
  t0 <- proc.time()[[3]]
  tp <- tangent_pca(X, d)
  t_tp <- proc.time()[[3]] - t0
  t0 <- proc.time()[[3]]
  cp <- cpns(X, d)
  t_cp <- proc.time()[[3]] - t0
  S <- list("SP-PCA" = f$scores, "Tangent PCA" = tp$scores, "CPNS" = cp$scores)
  res <- lapply(S, function(s) residual_scores(s, X, d, Q = 2))
  ari <- if (!is.null(cl)) {
    sapply(S, function(s) {
      km <- kmeans(s[, 1], centers = 3, nstart = 10)$cluster
      tab <- table(km, cl)
      a <- sum(choose(tab, 2)); b <- sum(choose(rowSums(tab), 2))
      c2 <- sum(choose(colSums(tab), 2)); N <- choose(length(cl), 2)
      (a - b * c2 / N) / ((b + c2) / 2 - b * c2 / N)
    })
  } else rep(NA, 3)
  data.frame(design = design, r = r, rep = i, method = names(S),
             expl1 = sapply(res, function(z) 1 - z[2] / z[1]),
             expl2 = sapply(res, function(z) 1 - z[3] / z[1]),
             ari = ari, time = c(t_sp, t_tp, t_cp),
             stress = c(f$smds$stress[2], NA, NA), rho = c(f$smds$rho, NA, NA))
}

cl <- makeCluster(as.integer(Sys.getenv("SPPCA_CORES", "4")))
clusterExport(cl, c("one_rep", "n"))
jobs <- do.call(rbind, lapply(seq_len(nrow(settings)), function(s)
  data.frame(i = 1:M, design = settings$design[s], r = settings$r[s])))
t0 <- Sys.time()
out <- parLapplyLB(cl, seq_len(nrow(jobs)), function(j, jobs) {
  tryCatch(one_rep(jobs$i[j], jobs$design[j], jobs$r[j]),
           error = function(e) data.frame(design = jobs$design[j], r = jobs$r[j],
                                          rep = jobs$i[j], method = NA,
                                          expl1 = NA, expl2 = NA, ari = NA,
                                          time = NA, stress = NA, rho = NA,
                                          err = conditionMessage(e)))
}, jobs = jobs)
stopCluster(cl)
res <- do.call(rbind, lapply(out, function(z) z[, intersect(names(z), c(
  "design", "r", "rep", "method", "expl1", "expl2", "ari", "time", "stress", "rho"))]))
saveRDS(res, "analysis/results/simulation_results.rds")
print(Sys.time() - t0)
agg <- aggregate(cbind(expl1, expl2, ari, time) ~ design + r + method, data = res,
                 FUN = mean, na.action = na.pass)
print(agg[order(agg$design, agg$r, agg$method), ], digits = 3)
errs <- Filter(function(z) !is.null(z$err), out)
if (length(errs)) print(unique(sapply(errs, function(z) z$err)))
