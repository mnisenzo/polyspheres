## Auxiliary functions

# Plot a sphere
sphere1.f <- function(x0 = 0,
                      y0 = 0,
                      z0 = 0,
                      r = 1,
                      n = 101,
                      ...) {
  
  f <- function(s, t) {
    
    cbind(
      r * cos(t) * cos(s) + x0,
      r * sin(s) + y0,
      r * sin(t) * cos(s) + z0
    )
    
  }
  
  persp3d(
    f,
    slim = c(-pi / 2, pi / 2),
    tlim = c(0, 2 * pi),
    n = n,
    add = TRUE,
    ...
  )
  
}

# Linear interpolation
lin_interp <- function(X, N = 10) {
  
  n_X <- nrow(X)
  X_new <- matrix(nrow = N * (n_X - 1) + 1, ncol = ncol(X))
  t0 <- seq(0, 1, l = N + 1)[-(N + 1)]
  t1 <- 1 - t0
  for (i in 1:(n_X - 1)) {
    X_new[N * (i - 1) + 1:N, ] <- t1 %x% X[i, , drop = FALSE] +
      t0 %x% X[i + 1, , drop = FALSE]
  }
  X_new[N * (n_X - 1) + 1, ] <- X[n_X, ]
  return(X_new)
  
}

# Returns Euclidean norm of a vector
norm_vec <- function(x) {
  
  return(sqrt(sum(x^2)))
  
}

# Returns indices of rows in first array that are not contained in second array
unique_indices <- function(arr1, arr2) {
  
  inds <- 1:nrow(arr1)
  diff <- intersect(arr1, arr2)
  for (i in inds) {
    if (length(intersect(arr1[i, ], diff)) == length(arr1[i, ])) {
      inds[i] <- 0
      
    }
    
  }
  
  return(inds[inds != 0])
  
}

# Returns indices of rows in first array that are contained in second array
subset_indices <- function(arr1, arr2) {
  
  inds <- 1:nrow(arr1)
  diff <- intersect(arr1, arr2)
  for (i in inds) {
    
    if (length(intersect(arr1[i, ], diff)) == 0) {
      
      inds[i] <- 0
      
    }
    
  }
  
  return(inds[inds != 0])
  
}

# Compare predicted 1D data with a given function using the formula RSS/TSS up to a flip
modified_R2 <- function(y, f) {
  
  n <- length(y)
  flip_y <- y[c(n:1)]
  result <- min(sum((y - f)^2), sum((flip_y - f)^2)) / (n * var(f))
  return(result)
  
}

# Repeat a set of objects n times in order
multiple_reps <- function(obj_list, n){
  
  l = length(obj_list)
  out <- rep(0, n*l)
  for(i in 1:l){
    out[((i-1)*n+1):(i*n)] = rep(obj_list[i],n)
  }
  return(out)
}

# Create a gradient for n objects
col_grad <- function(n, start=0,stop =.65){
  out <- rep(0,n)
  rgbs <- lin_interp(matrix(rep(c(start,stop),3),nrow=2), n)
  cols <- apply(rgbs,1,{function(x) rgb(x[1],x[2],x[3])})
  return(cols)
}
