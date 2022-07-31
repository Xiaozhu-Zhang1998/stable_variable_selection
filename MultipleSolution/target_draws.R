# make "target draws"

source("beta_trans.R")

target_draws = function(X, orthant, x_par, npoints = 10000, walk = 'BiW', walk_length = 100) {
  # preparation
  n = nrow(X)
  d = ncol(X)
  A = rbind(X %*% diag(orthant), 
            rep(1, d))
  V = pracma::nullspace(A)
  # draw points
  P = volesti::Hpolytope(A = -V, b = as.vector(abs(x_par) ))
  random_walk = list(
    walk = walk,
    walk_length = walk_length
  )
  alpha = volesti::sample_points(P, n = npoints, random_walk = random_walk)
  x = abs(x_par) %*%  matrix(1, nrow = 1, ncol = npoints) + V %*% alpha
  soln = t(x) %*% diag(orthant)
}
