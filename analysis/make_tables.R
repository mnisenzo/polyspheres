## LaTeX tables for the paper
bs <- function(x) gsub("@", "\\\\", x)  # write "@" for a LaTeX backslash

## Simulation table
tab <- readRDS("analysis/results/sim_table.rds")
lab <- c(great_full = "Great circles", small_arc = "Small-circle arcs",
         active2 = "Two active spheres", clusters = "Clusters")
meths <- c("SP-PCA", "Tangent PCA", "CPNS")
L <- c("@begin{table}[t]", "@centering@small",
       "@begin{tabular}{llrrrrrrrrr}", "@toprule",
       paste0(" & & @multicolumn{3}{c}{$@mathrm{EV}(1)$ (@%)} & ",
              "@multicolumn{3}{c}{$@mathrm{EV}(2)$ (@%)} & ",
              "@multicolumn{3}{c}{Time (s)}@@"),
       "@cmidrule(lr){3-5}@cmidrule(lr){6-8}@cmidrule(lr){9-11}",
       "Design & $r$ & SP & Tan & CPNS & SP & Tan & CPNS & SP & Tan & CPNS@@",
       "@midrule")
for (ds in names(lab)) {
  rs <- sort(unique(tab$r[tab$design == ds]))
  for (r in rs) {
    g <- function(m, v) tab[tab$design == ds & tab$r == r & tab$method == m, v]
    e1 <- sapply(meths, g, v = "expl1") * 100
    e2 <- sapply(meths, g, v = "expl2") * 100
    tt <- sapply(meths, g, v = "time")
    f <- function(x) {
      s <- sprintf("%.1f", x)
      b <- which(round(x, 1) == max(round(x, 1)))
      s[b] <- paste0("@textbf{", s[b], "}")
      s
    }
    L <- c(L, paste0(paste(c(if (r == rs[1]) lab[ds] else "", r, f(e1), f(e2),
                             sprintf("%.2f", tt)), collapse = " & "), "@@"))
  }
}
L <- c(L, "@bottomrule", "@end{tabular}",
       paste0("@caption{Mean leave-one-out explained variance of the first one ",
              "and two scores, and mean fitting time, over $20$ replicates ",
              "($n=200$). SP: SP-PCA; Tan: tangent PCA. Best values in bold.}"),
       "@label{tab:sim}", "@end{table}")
writeLines(bs(L), "analysis/tab_sim.tex")

## Hippocampus explained-variance table
res <- readRDS("analysis/results/hippo_resid.rds")
ev <- 100 * (1 - res$loo / res$loo[, 1])
H <- c("@begin{table}[t]", "@centering@small", "@begin{tabular}{lrrrrr}",
       "@toprule", "Method & $q=1$ & $q=2$ & $q=3$ & $q=4$ & $q=5$@@",
       "@midrule")
for (m in rownames(ev)) {
  x <- ev[m, 2:6]
  best <- apply(ev[, 2:6], 2, max)
  s <- sprintf("%.1f", x)
  s[round(x, 1) == round(best, 1)] <- paste0("@textbf{", s[round(x, 1) == round(best, 1)], "}")
  H <- c(H, paste0(paste(c(m, s), collapse = " & "), "@@"))
}
H <- c(H, "@bottomrule", "@end{tabular}",
       paste0("@caption{Hippocampus s-reps: leave-one-out explained variance ",
              "$@mathrm{EV}(q)$ (@%) of the first $q$ scores, using the common ",
              "local linear inversion~@eqref{eq:ev}.}"),
       "@label{tab:hippo}", "@end{table}")
writeLines(bs(H), "analysis/tab_hippo.tex")
cat(readLines("analysis/tab_sim.tex")[9:11], sep = "\n")
cat(readLines("analysis/tab_hippo.tex")[7:9], sep = "\n")

## Benchmark table (old GitHub implementation vs sppca)
if (file.exists("analysis/results/benchmark.rds")) {
  b <- readRDS("analysis/results/benchmark.rds")
  ag <- do.call(rbind, lapply(split(b, b$n), function(z) data.frame(
    n = z$n[1], t_old = mean(z$t_old), t_new = mean(z$t_new),
    t_smds = mean(z$t_new_smds), t_full = mean(z$t_full),
    s_old = mean(z$stress_old), s_new = mean(z$stress_new),
    s_full = mean(z$stress_full))))
  fmt_t <- function(x) ifelse(is.na(x), "--", sprintf("%.2f", x))
  fmt_s <- function(x) ifelse(is.na(x), "--", formatC(x, format = "e", digits = 1))
  B <- c("@begin{table}[t]", "@centering@small", "@begin{tabular}{rrrrrrrr}",
         "@toprule",
         paste0(" & @multicolumn{4}{c}{Time (s)} & @multicolumn{3}{c}{Stress}@@"),
         "@cmidrule(lr){2-5}@cmidrule(lr){6-8}",
         "$n$ & Original & SP-PCA & Scaling & Full SMDS & Original & SP-PCA & Full SMDS@@",
         "@midrule")
  for (i in seq_len(nrow(ag))) {
    B <- c(B, paste0(paste(c(ag$n[i], fmt_t(ag$t_old[i]), fmt_t(ag$t_new[i]),
                             fmt_t(ag$t_smds[i]), fmt_t(ag$t_full[i]),
                             fmt_s(ag$s_old[i]), fmt_s(ag$s_new[i]),
                             fmt_s(ag$s_full[i])), collapse = " & "), "@@"))
  }
  B <- c(B, "@bottomrule", "@end{tabular}",
         paste0("@caption{Mean running time and normalized full-sample stress ",
                "over $3$ replicates for circles on $(@Sp^2)^4$ ($p=8$, ",
                "$@ell=n/10$ landmarks). Original: our earlier implementation ",
                "(landmark SMDS with pointwise @texttt{optim}, PNS with ",
                "sequential tests). SP-PCA: the whole pipeline in ",
                "@texttt{sppca}, including PNS and bandwidth selection; ",
                "Scaling: its scaled MDS step. Full SMDS: all observations as ",
                "landmarks. A dash means not run (too slow) or not available.}"),
         "@label{tab:benchmark}", "@end{table}")
  writeLines(bs(B), "analysis/tab_benchmark.tex")
}
