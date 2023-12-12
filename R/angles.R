## Angles to Spheres transformation
angles_to_sph <- function(theta) {
  
  if (!is.matrix(theta)) {
    
    theta <- matrix(theta, nrow = 1)
    
  }
  
  d <- ncol(theta)
  x <- matrix(0, nrow = nrow(theta), ncol = d + 1)
  cosTheta <- cos(theta)
  sinTheta <- sin(theta)
  if (d == 1) {
    
    x <- cbind(cosTheta, sinTheta)
    
  } else {
    
    prodSinTheta <- t(apply(sinTheta, 1, cumprod))
    x <- cbind(
      cosTheta[, 1, drop = FALSE], prodSinTheta[, -d, drop = FALSE] *
        cosTheta[, -1, drop = FALSE],
      prodSinTheta[, d, drop = FALSE]
    )[, , drop = TRUE]
    
  }
  
  return(x)
  
}

## Rectangular to polysphere coordinate transformation
polysph_to_angles <- function(x, dims) {
  if (!is.matrix(x)) {
    
    x <- matrix(x, nrow = 1)
    
  }
  theta <-
    matrix(0, nrow = nrow(x), ncol = sum(dims))
  ind_p <- 0
  ind_r <- 0
  for (d in dims) {
    
    input <- as.matrix(x[, (ind_p + 1):(ind_p + d + 1)])
    if (nrow(x) == 1) {
      
      input <- t(input)
    }
    
    theta[, (ind_r + 1):(ind_r + d)] <- sphere_to_angles(input)
    ind_p <- ind_p + d + 1
    ind_r <- ind_r + d
    
  }
  
  return(theta)
  
}

## Turn degrees into radians
deg2rad <- function(deg) {
  
  (deg * pi) / (180)
  
}

## Polysphere to rectangular coordinate transformation
angles_to_polysph <- function(theta, dims) {
  
  if (!is.matrix(theta)) {
    
    theta <- matrix(theta, nrow = 1)
    
  }
  
  x <-
    matrix(0, nrow = nrow(theta), ncol = length(dims) + sum(dims))
  ind_p <- 0
  ind_r <- 0
  for (d in dims) {
    
    input <- as.matrix(theta[, (ind_p + 1):(ind_p + d)])
    if (nrow(theta) == 1) {
      
      input <- t(input)
      
    }
    
    x[, (ind_r + 1):(ind_r + d + 1)] <- angles_to_sphere(input)
    ind_p <- ind_p + d
    ind_r <- ind_r + d + 1
    
  }
  
  return(x)
  
}

## Toroidal distance
torus_dist <- function(a, b) {
  
  dist.coord <- rep(NA, length(a))
  for (i in seq_along(a)) {
    dist1 <- abs(a[i] - b[i])
    dist2 <- 2 * pi - abs(a[i] - b[i])
    dist.coord[i] <- min(dist1, dist2)^2
    
  }
  
  dist.torus <- sqrt(sum(dist.coord))
  return(dist.torus)
  
}

## Compute distance between two points on polysphere
polysph_dist <-
  function(a,
           b,
           angular = TRUE,
           dims = rep(1, length(a))) {
    
    if (angular) {
      
      a <- angles_to_polysph(a, dims)
      b <- angles_to_polysph(b, dims)
      
    }
    
    dist.coord <- rep(NA, length(dims))
    ind <- 0
    for (i in seq_along(dims)) {
      
      d <- dims[i]
      a1 <- a[(ind + 1):(ind + d + 1)]
      b1 <- b[(ind + 1):(ind + d + 1)]
      ind <- ind + d + 1
      dist.coord[i] <- acos(a1 %*% b1)
      
    }

    dist.sphere <- sqrt(sum(dist.coord**2))
    return(dist.sphere)
    
  }
