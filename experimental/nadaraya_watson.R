# Nadaraya-Watson estimator from the sphere to the polysphere
mNW <- function(x, X, Y, dims = rep(1,ncol(X)), h=.5, K = dnorm) {
  
  # Arguments
  # x: matrix (n rows) of evaluation points in R^d
  # X: matrix (n rows) with points in S^d
  # Y: matrix (n rows) with points in S^(d_1,...d_r)
  # dims: (d_1,...,d_r) are the dimensions of the polysphere
  # h: bandwidth
  # K: kernel
  
  # Matrix of size n x length(x)
  n <- nrow(X)
  Kx <- matrix(0, nrow(x),n)
  for(i in 1:nrow(x)){
    for(j in 1:n){
      Kx[i,j] <- K(dist_polysph_cross(t(matrix(x[i,])),t(matrix(X[j,])), ind_dj = comp_ind_dj(sum(dims))) / h) / h
    }
  }

  # Weights
  W <- Kx / rowSums(Kx) # Column recycling!
  # Means at x ("drop" to drop the matrix attributes)
  res<- drop(W %*% Y)
  for(i in 1:nrow(x)){
    res[i,]<- normalize_poly(res[i,],dims)
  }
  return(res)
}

# Cross validation for Nadaraya-Watson estimator from sphere to polysphere
cvNW <- function(X, Y, dims = dims, h, K = dnorm) {
  n <- nrow(X)
  Kx <- matrix(0, n,n)
  for(i in 1:n){
    for(j in 1:n){
      Kx[i,j] <- K(dist_polysph_cross(t(matrix(X[i,])),t(matrix(X[j,])), ind_dj = comp_ind_dj(sum(dims))) / h) / h
    }
  }
  results <- rep(0,n)
  for(i in 1:n){
    # Renormalize weights
    weights <- Kx[i,]/(1- Kx[i,i])
    weights <- weights[-i]
    guess <- weights %*% Y[-i,]
    guess <- normalize_poly(guess,dims)
    results[i] <- dist_polysph_cross(t(matrix(guess)),t(matrix(Y[i,])), ind_dj = comp_ind_dj(dims))
    
  }
  return(sum(results^2))

}


# Nadaraya-Watson estimator from the sphere to the polysphere with optimized bandwidth
opt_mNW <- function(x, X, Y, dims = rep(1,ncol(X)), K = dnorm,tol = 5e-2) {
  opt_h <- optim(
    fn = cvNW,
    par = c(.5), X= X, Y = Y, dims = dims, K=K,
    method = "BFGS",
    control = list(reltol = tol)
  )
  result <- mNW(x=x, X=X, Y=Y, dims = dims, h=opt_h$par, K = dnorm) 

}

# Find optimum CV bandwidth, with sensible grid
bw.cv.grid <- function(X, Y, dims = rep(1,ncol(X)),
                       h.grid = diff(range(X)) * (seq(0.1, 1, l = 200))^2,
                       K = dnorm, plot.cv = FALSE) {
  
  obj <- sapply(h.grid, function(h) cvNW(X = X, Y = Y, dims = dims, h = h, K = K))
  h <- h.grid[which.min(obj)]
  if (plot.cv) {
    plot(h.grid, obj, type = "o")
    rug(h.grid)
    abline(v = h, col = 2, lwd = 2)
  }
  h
  
}


