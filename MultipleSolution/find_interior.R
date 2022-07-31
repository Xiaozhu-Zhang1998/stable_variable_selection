# Find the interior points

source("beta_trans.R")

find_interior = function(X, beta, tol2) {
  # define A and b
  d = ncol(X)
  A = rbind( cbind(X, -X), rep(1, 2 * d))
  b = rbind(X %*% beta, norm(matrix(beta), '1'))
  n = nrow(A)
  
  # find the y
  y_par = c(rep(1, n-1), - sum(b[1:(n-1)]) / b[n])
  V = pracma::nullspace(t(b))
  
  A_ = - t(A) %*% V
  b_ = t(A) %*% y_par
  P = volesti::Hpolytope(A = A_, b = as.vector(b_) + abs(rnorm(2 * d, sd = 1e-10)))
  random_walk = list(
    walk = 'BiW',
    walk_length = 1
  )
  alpha = volesti::sample_points(P, n = 100, random_walk = random_walk)
  y = y_par + V %*% alpha[, 100]
  
  # find the V
  z = t(A) %*% y
  sz = z > tol2
  V = diag(rep(1, 2 * d))[, !sz]
  
  # find the x_par
  Pr = t(pracma::orth(A %*% V))
  v = pracma::pinv(Pr %*% A %*% V) %*% (Pr %*% b)
  x_par = beta_tight(V %*% v)
  
  return(x_par)
}
