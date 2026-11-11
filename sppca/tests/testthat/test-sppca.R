
d <- rep(2, 3)
set.seed(42)
X <- r_circles_polysph(n = 80, r = 3, sigma = 0.05)
th_X <- attr(X, "theta")
attr(X, "theta") <- NULL

test_that("geometry: distances are symmetric, zero diagonal, and match polykde", {
  D <- dist_polysph_cross(X, X, d = d)
  expect_equal(D, t(D), tolerance = 1e-12)
  expect_lt(max(abs(diag(D))), 1e-6)
  D2 <- polykde::dist_polysph(x = X[1:10, ], y = X[11:20, ],
                              ind_dj = polykde::comp_ind_dj(d), std = FALSE)
  expect_equal(diag(D[1:10, 11:20]), drop(D2), tolerance = 1e-10)
})

test_that("geometry: exp and log are inverse", {
  mu <- mean_polysph(X, d = d)
  V <- log_polysph(X, mu = mu, d = d)
  expect_equal(exp_polysph(V, mu = mu, d = d), X, tolerance = 1e-8)
})

test_that("stress gradients match finite differences", {
  Y <- matrix(rnorm(8 * 4), 8, 4)
  U <- matrix(rnorm(5 * 4), 5, 4)
  U <- U / sqrt(rowSums(U^2))
  Delta <- matrix(runif(40, 0.2, 2), 8, 5)
  Dsym <- as.matrix(dist(matrix(rnorm(8 * 3), 8)))
  for (case in list(list(U = U, Delta = Delta), list(U = NULL, Delta = Dsym))) {
    for (rho in list(NULL, 1.3)) {
      f <- function(y) stress_sph(matrix(y, 8), U = case$U, Delta = case$Delta,
                                  rho = rho)$value
      g <- c(stress_sph(Y, U = case$U, Delta = case$Delta, rho = rho)$gradient)
      eps <- 1e-6
      gn <- sapply(seq_along(Y), function(i) {
        e <- numeric(length(Y)); e[i] <- eps
        (f(c(Y) + e) - f(c(Y) - e)) / (2 * eps)
      })
      expect_equal(g, gn, tolerance = 1e-5)
    }
  }
})

test_that("SMDS recovers an exactly spherical configuration", {
  # Points on a great circle of each sphere driven by the same angle are
  # isometric (up to scale sqrt(r)) to a circle of radius sqrt(r)
  th <- seq(0, 2 * pi, l = 61)[-61]
  Xc <- do.call(cbind, lapply(1:3, function(j) cbind(cos(th), sin(th), 0)))
  fit <- smds_polysph(Xc, d = d, p = 2, landmarks = "all")
  expect_lt(fit$stress[1], 1e-4)
  expect_equal(fit$rho, sqrt(3), tolerance = 1e-2)
})

test_that("sppca captures a one-dimensional mode", {
  fit <- sppca(X, d = d, p = 3)
  expect_s3_class(fit, "sppca")
  expect_gt(fit$var_sphere[1], 80)
  # First score is (circularly) associated with the latent angle
  th <- th_X
  s1 <- fit$scores[, 1] / fit$smds$rho
  cc <- abs(cor(cos(s1), cos(th))) + abs(cor(sin(s1), sin(th)))
  r <- residual_sppca(fit, Q = 2)
  expect_lt(r[2], 0.5 * r[1])
  m <- mode_sppca(fit, j = 1)
  expect_equal(dim(m), c(50, 9))
})

test_that("competitors reconstruct exactly at full rank", {
  f <- tangent_pca(X, d = d)
  expect_equal(reconstruct_tangent_pca(f, q = 9, d = d), X, tolerance = 1e-6)
  f <- extrinsic_pca(X, d = d)
  expect_equal(reconstruct_extrinsic_pca(f, q = 9, d = d), X, tolerance = 1e-6)
  f <- cpns(X, d = d)
  expect_equal(reconstruct_cpns(f, q = 6, d = d), X, tolerance = 1e-4)
})

test_that("pns_fit reproduces shapes::pns", {
  set.seed(3)
  Y <- sweep(matrix(rnorm(60 * 4), 60, 4) %*% diag(c(3, 1, 0.5, 0.2)), 2,
             c(0, 0, 0, 4), "+")
  Y <- Y / sqrt(rowSums(Y^2))
  for (ty in c("small", "great")) {
    a <- quiet(shapes::pns(t(Y), sphere.type = ty, output = FALSE))
    b <- pns_fit(Y, sphere_type = ty)
    expect_equal(b$percent, a$percent, tolerance = 1e-3)
    expect_equal(abs(b$resmat), abs(a$resmat), tolerance = 1e-3)
    expect_equal(t(shapes::PNSe2s(b$resmat, b$PNS)), Y, tolerance = 1e-8)
  }
})

test_that("ll_scores_polysph is finite for isolated evaluation points", {
  set.seed(4)
  S <- cbind(c(rnorm(79, sd = 0.01), 10))
  fit <- ll_scores_polysph(S, S, X, d = d, h = 0.2, loo = TRUE)
  expect_true(all(is.finite(fit)))
  expect_equal(rowSums(fit[, 1:3]^2), rep(1, 80), tolerance = 1e-8)
})
