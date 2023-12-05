# Nadaraya-Watson estimator from the sphere to the polysphere
mNW <- function(x, X, Y, radius=1, dims = rep(1,ncol(X)), h=.5, dist_x, K = dnorm) {
  
  # Arguments
  # x: matrix (n rows) of evaluation points in R^d
  # X: matrix (n rows) with points in S^d
  # Y: matrix (n rows) with points in S^(d_1,...d_r)
  # dims: (d_1,...,d_r) are the dimensions of the polysphere
  # h: bandwidth
  # K: kernel

  # Set bandwidth to positive
  n = nrow(X)
  h <- abs(h)
  Kx <- matrix(K(radius*dist_x / h) / h, ncol = n)
  
  # Weights
  if(nrow(Kx) > 1){
    W <- Kx / rowSums(Kx) # Column recycling!
  } else{
    W <- Kx / sum(Kx) # Column recycling!
  }

  # Account for na entries in weights
  W[is.na(W)] <-0

  # Means at x
  res<- W %*% Y
  res <- normalize_poly(res,dims)
  return(res)
}

# Cross validation for Nadaraya-Watson estimator from sphere to polysphere
cvNW <- function(X, Y, radius = 1, dims = dims, h, dist_X, K = dnorm) {
  n <- nrow(X)
  results <- rep(0,n)
  # Set bandwidth to positive
  h <- abs(h)
  #print(h)
  
  Kx <- K(radius*dist_X / h) / h
  guesses <- matrix(0, nrow = n, ncol = sum(dims)+ length(dims))
  for(i in 1:n){
    # Renormalize weights
     guess <- mNW(x, X[-i,], Y[-i,], radius=radius, dims = dims, h=h, dist_x = dist_X[i,-i], K = K)
     guesses[i,] <- guess
  }
  guesses <- normalize_poly(guesses,dims)
  for(i in 1:n){
    results[i] <- dist_polysph_cross(t(matrix(guesses[i,])),t(matrix(Y[i,])), ind_dj = comp_ind_dj(dims))
  }

  # Weigh CV distances
  out <- rep(0,n)
  weights <- matrix(0,nrow = 1, ncol = n)
  #for(i in 1:n){
     # out[i] <- results[i]
    #weight_i <- K(0)/sum(Kx[i,])
   # weight_i <- K(0)/ colSums(K(outer(guesses, guesses, "-") / h))
    #out[i] <- results[i]/(1 - weight_i)
  #s}
  
  out <- results/(1 - K(0)/ colSums(K(outer(guesses, guesses, "-") / h)))
  # Penalty to minimize average distances
  return(sum(out^2))
}


# Nadaraya-Watson estimator from the sphere to the polysphere with optimized bandwidth
opt_mNW <- function(x, X, Y, radius = 1, dims = rep(1,ncol(X)), init_bw=.5,K = dnorm,tol = 5e-2) {
  # Distance matrix between spherical coordinates
  dist_X <- dist_polysph_cross(X,X, ind_dj = comp_ind_dj(sum(dims)))
  # Distance matrix between principal curve and spherical coordinates
  dist_x <- dist_polysph_cross(x,X, ind_dj = comp_ind_dj(sum(dims)))
  
  # Weight constant based on norm in cross validation
  opt_h <- optim(
    fn = cvNW,
    par = c(h=init_bw), X= X, Y = Y,radius = radius, dist_X= dist_X, dims = dims, K=K,
    method = "BFGS",
    control = list(reltol = tol)
  )
  opt_h <- opt_h$par
  print(opt_h)
  # opt_h <- bw.cv.grid(X=X,Y=Y,dist_X=dist_X, dims = dims,h.grid = diff(range(X)) * (seq(0.1, 1, l = 200))^2,
  #                     K=K)
  # print(opt_h)
  
  result <- mNW(x=x, X=X, Y=Y, dims = dims, dist_x = dist_x, h=opt_h, K = dnorm) 
  return(result)
}

# Find optimum CV bandwidth, with sensible grid
bw.cv.grid <- function(X, Y, dist_X, dims = rep(1,ncol(X)),
                       h.grid = diff(range(X)) * (seq(0.05, 1, l = 200))^2,
                       K = dnorm, plot.cv = FALSE) {
  
  obj <- sapply(h.grid, function(h) cvNW(X = X, Y = Y,dist_X=dist_X, dims = dims, h = h, K = K))
  h <- h.grid[which.min(obj)]
  if (plot.cv) {
    plot(h.grid, obj, type = "o")
    rug(h.grid)
    abline(v = h, col = 2, lwd = 2)
  }
  return(h)
  
}


