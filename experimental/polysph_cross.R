library(polykde)
library(sphunif)


set.seed(5)
n <- 50
tol <- 5e-3
landmark_tol <- 5e-3
nspheres <- 2
dims <- c(1, 1)
r_star <- sqrt(length(dims))
ind_dj <- comp_ind_dj(dims)
landmark_size <- .2
ls <- n * landmark_size

samp_1 <- r_path_s1r(n = n, r = nspheres, k = rep(1, nspheres), angles = TRUE)
cart_x <- polyAnglesToSphere(samp_1, dims)
dist <- dist_polysph_matrix(cart_x, ind_dj = ind_dj, std = FALSE)
d_ij_polysphere <- dist[upper.tri(dist)]
landmarks <- get_center(x = cart_x, dims = dims, angular = angular, landmark_size = landmark_size)
y <-
  cmdscale(
    d = dist,
    k = d_sphere + 1,
    eig = FALSE,
    x.ret = FALSE
  )
y <- y / sqrt(rowSums(y^2))
y_initial <- y
r_star_initial <- r_star
not_landmarks <- c(1:nrow(samp_1))
not_landmarks <- which(!not_landmarks %in% landmarks)

# Objective function
stress <- function(y, r_star) {
  r_star <- y[n * (d_sphere + 1) + 1] # r_star in optimization
  y <-
    matrix(y[1:(n * (d_sphere + 1))], nrow = n, ncol = d_sphere + 1)
  y <- y / sqrt(rowSums(y^2)) # Normalization of y

  dim(y) <- c(dim(y), 1) # As an array
  # Spherical distance
  ind_tri <- sphunif::upper_tri_ind(n)
  d_ij_sphere <- sphunif::Psi_mat(
    data = y,
    use_ind_tri = TRUE,
    ind_tri = ind_tri
  )

  d_ij_sphere <- drop(d_ij_sphere)
  # Stress function (only sum the distinct distances and disregard diagonal)
  return(sum((d_ij_polysphere - r_star * d_ij_sphere)^2) / length(d_ij_polysphere))
}

landmark_stress <- function(y, ls, d_sphere, dist, landmarks, r_star) {
  r_star <- y[ls * (d_sphere + 1) + 1] # r_star in optimization
  y <-
    matrix(y[1:(ls * (d_sphere + 1))], nrow = ls, ncol = d_sphere + 1)
  y <- y / sqrt(rowSums(y^2)) # Normalization of y

  dim(y) <- c(dim(y), 1) # As an array
  # Spherical distance
  ind_tri <- sphunif::upper_tri_ind(ls)
  d_ij_sphere <- sphunif::Psi_mat(
    data = y,
    use_ind_tri = TRUE,
    ind_tri = ind_tri
  )
  d_ij_sphere <- drop(d_ij_sphere)
  dist_sub <- dist[landmarks, landmarks]
  landmark_polysphere <- dist_sub[upper.tri(dist_sub)]
  # Stress function (only sum the distinct distances and disregard diagonal)
  return(sum((landmark_polysphere - r_star * d_ij_sphere)^2) / length(landmark_polysphere))
}
# Solves MDS for the non-landmark coordinates yn given landmark coordinates yl
remaining_stress <- function(yn, yl, d_sphere, n, ls, landmarks, not_landmarks, dist) {
  # landmark points defined outside the function
  yn <-
    matrix(yn[1:((n - ls) * (d_sphere + 1))], nrow = n - ls, ncol = d_sphere + 1)
  y1 <- rbind(yn, yl)
  y1 <-
    matrix(y1[1:(n * (d_sphere + 1))], nrow = n, ncol = d_sphere + 1)
  y1 <- y1 / sqrt(rowSums(y1^2)) # Normalization of y

  dim(y1) <- c(dim(y1), 1) # As an array
  # Spherical distance
  ind_tri <- sphunif::upper_tri_ind(n)
  d_ij_sphere <- sphunif::Psi_mat( # coerce to dist class
    data = y1,
    use_ind_tri = TRUE,
    ind_tri = ind_tri
  )
  sph_mat <- matrix(0, n, n)
  sph_mat[upper.tri(sph_mat, diag = FALSE)] <- d_ij_sphere
  # polysphere-kde
  sph_mat <- drop(sph_mat)[1:(n - ls), (n - ls + 1):n]
  ind_dj <- comp_ind_dj(d_sphere)
  s1 <- y[1:(n - ls), ]
  s2 <- y[(n - ls + 1):n, ]
  sph_mat1 <- dist_polysph_cross(x = s1, y = s2, ind_dj = ind_dj, std = FALSE)
  # print(max(sph_mat- sph_mat1))
  ind_dj <- comp_ind_dj(d_sphere)
  # sph_mat <- dist_polysph_cross(x=y[1:(n-ls),],y=y[(n - ls + 1):n,],ind_dj=ind_dj,std=FALSE)
  landmark_polysphere <- dist[not_landmarks, landmarks]
  # Stress function (only sum the distinct distances and disregard diagonal)
  return(sum((landmark_polysphere - r_star * sph_mat)^2) / length(landmark_polysphere))
}

# Solve landmark MDS first
res_landmark <- optim(
  fn = landmark_stress,
  par = c(y[landmarks, ], r_star), ls = ls, d_sphere = d_sphere, dist = dist, landmarks = landmarks,
  method = "BFGS",
  control = list(reltol = landmark_tol)
)
r_star <- res_landmark$par[ls * (d_sphere + 1) + 1] # r_star in optimization
yl <-
  matrix(res_landmark$par[1:(ls * (d_sphere + 1))], nrow = ls, ncol = d_sphere + 1)
y[landmarks, ] <- yl
# If there are non-landmarks left, run the analysis on them
if (n - ls > 0) {
  res <- optim(
    fn = remaining_stress,
    par = c(y[not_landmarks, ]), yl = yl, d_sphere = d_sphere, n = n, ls = ls, landmarks = landmarks, not_landmarks = not_landmarks, dist = dist,
    method = "BFGS",
    control = list(reltol = tol)
  )
  yn <-
    matrix(res$par, nrow = n - ls, ncol = d_sphere + 1)
  y[not_landmarks, ] <- yn
}

set.seed(1)
r <- 2
d <- c(2, 3)
n1 <- 5
n2 <- 2
X1 <- r_unif_polysph(n = n1, d = d)
X2 <- r_unif_polysph(n = n2, d = d)
ind_dj <- comp_ind_dj(d)
dist_polysph_cross(x = X1, y = X2, ind_dj = ind_dj, std = FALSE)
x1 <- X1 / sqrt(rowSums(X1^2))
x11 <- x1
dim(x11) <- c(dim(x11), 1) # As an array
# Spherical distance
ind_tri <- sphunif::upper_tri_ind(n1)
d_ij_sphere <- sphunif::Psi_mat( # coerce to dist class
  data = x11,
  use_ind_tri = TRUE,
  ind_tri = ind_tri
)

sph_mat <- matrix(0, n1, n1)
sph_mat[upper.tri(sph_mat, diag = FALSE)] <- d_ij_sphere
sph_mat2 <- drop(sph_mat)[1:3, 4:5]
ind_dj <- comp_ind_dj(6)
s1 <- x1[1:3, ]
s2 <- x1[4:5, ]
sph_mat1 <- dist_polysph_cross(x = s1, y = s2, ind_dj = ind_dj, std = FALSE)
# print(max(sph_mat- sph_mat1))
ind_dj <- comp_ind_dj(d_sphere)



set.seed(5)
n <- 50
tol <- 5e-3
landmark_tol <- 5e-3
nspheres <- 2
dims <- c(1, 1)
r_star <- sqrt(length(dims))
ind_dj <- comp_ind_dj(dims)
landmark_size <- .2
ls <- n * landmark_size




samp_1 <- r_path_s1r(n = n, r = nspheres, k = rep(1, nspheres), angles = TRUE)
cart_x <- polyAnglesToSphere(samp_1, dims)
dist <- dist_polysph_matrix(cart_x, ind_dj = ind_dj, std = FALSE)
d_ij_polysphere <- dist[upper.tri(dist)]
landmarks <- get_center(x = cart_x, dims = dims, angular = angular, landmark_size = landmark_size)
y <-
  cmdscale(
    d = dist,
    k = d_sphere + 1,
    eig = FALSE,
    x.ret = FALSE
  )
y <- y / sqrt(rowSums(y^2))
y_initial <- y
r_star_initial <- r_star


landmarks <- get_centroids(x = x, dims = dims, angular = angular, landmark_size = landmark_size)
landmarks <- landmarks / sqrt(rowSums(landmarks^2))
# project the landmarks onto the sphere
dist_ls <- dist_polysph_matrix(landmarks, ind_dj = ind_dj, std = FALSE)
d_ij_polysphere_ls <- dist_ls[upper.tri(dist_ls)]
y_ls <-
  cmdscale(
    d = dist_ls,
    k = d_sphere + 1,
    eig = FALSE,
    x.ret = FALSE
  )
d_sphere <- ncol(y_ls) - 1
y_ls <- y_ls / sqrt(rowSums(y_ls^2))
stress <- function(y, r_star, n, d_sphere, d_ij_polysphere) {
  r_star <- y[n * (d_sphere + 1) + 1] # r_star in optimization
  y <-
    matrix(y[1:(n * (d_sphere + 1))], nrow = n, ncol = d_sphere + 1)
  y <- y / sqrt(rowSums(y^2)) # Normalization of y

  dim(y) <- c(dim(y), 1) # As an array
  # Spherical distance
  ind_tri <- sphunif::upper_tri_ind(n)
  d_ij_sphere <- sphunif::Psi_mat(
    data = y,
    use_ind_tri = TRUE,
    ind_tri = ind_tri
  )

  d_ij_sphere <- drop(d_ij_sphere)
  # Stress function (only sum the distinct distances and disregard diagonal)
  return(sum((d_ij_polysphere - r_star * d_ij_sphere)^2) / length(d_ij_polysphere))
}

landmark_stress <- function(y, d_sphere, r_star) {
  r_star <- y[ls * (d_sphere + 1) + 1] # r_star in optimization
  y <-
    matrix(y[1:(ls * (d_sphere + 1))], nrow = ls, ncol = d_sphere + 1)
  y <- y / sqrt(rowSums(y^2)) # Normalization of y

  dim(y) <- c(dim(y), 1) # As an array
  # Spherical distance
  ind_tri <- sphunif::upper_tri_ind(ls)
  d_ij_sphere <- sphunif::Psi_mat(
    data = y,
    use_ind_tri = TRUE,
    ind_tri = ind_tri
  )
  d_ij_sphere <- drop(d_ij_sphere)
  # Stress function (only sum the distinct distances and disregard diagonal)
  return(sum((d_ij_polysphere_ls - r_star * d_ij_sphere)^2) / length(d_ij_polysphere_ls))
}

# Solves MDS for y given landmark centroids ls
remaining_stress <- function(y, d_sphere, dims, n, yl, r_star, dist) {
  # landmark points defined outside the function
  y1 <-
    matrix(y[1:(n * (d_sphere + 1))], nrow = n, ncol = d_sphere + 1)
  y1 <- y1 / sqrt(rowSums(y1^2)) # Normalization of y
  y2 <- y1
  dim(y2) <- c(dim(y2), 1) # As an array
  # Spherical distance
  ind_tri <- sphunif::upper_tri_ind(n)
  # polysphere-kde
  ind_dj <- comp_ind_dj(d_sphere)
  sph_mat <- dist_polysph_cross(x = y1, y = yl, ind_dj = ind_dj, norm_x = FALSE, norm_y = FALSE, std = FALSE)
  landmark_polysphere <- dist
  # Stress function (only sum the distinct distances and disregard diagonal)
  return(sum((landmark_polysphere - r_star * sph_mat)^2) / length(landmark_polysphere))
}
# Solve landmark MDS first
res_landmark <- optim(
  fn = landmark_stress,
  par = c(y_ls, r_star), d_sphere = d_sphere,
  method = "BFGS",
  control = list(reltol = landmark_tol)
)
r_star <- res_landmark$par[ls * (d_sphere + 1) + 1] # r_star in optimization
yl <-
  matrix(res_landmark$par[1:(ls * (d_sphere + 1))], nrow = ls, ncol = d_sphere + 1)
yl <- yl / sqrt(rowSums(yl^2)) # Normalization of y

ind_dj_poly <- comp_ind_dj(dims)
# distances between original points and landmarks on the polysphere
dist_initial <- dist_polysph_cross(x = cart_x, y = landmarks, norm_x = TRUE, norm_y = TRUE, ind_dj = ind_dj_poly, std = FALSE)
y <-
  cmdscale(
    d = dist,
    k = d_sphere + 1,
    eig = FALSE,
    x.ret = FALSE
  )
y <- y / sqrt(rowSums(y^2))
y_initial <- y
# Compare landmarks to non-landmarks
res1 <- optim(
  fn = remaining_stress,
  par = c(y), r_star = r_star, yl = yl, d_sphere = d_sphere, n = n, dist = dist_initial,
  method = "BFGS",
  control = list(reltol = 5e-5)
)
# y <- matrix(res$par, nrow = n, ncol = d_sphere + 1)
