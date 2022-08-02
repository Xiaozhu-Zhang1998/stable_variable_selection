# make "target draws"

source("beta_trans.R")

target_draws = function(A, FR, beta, npoints = 10000, walk = 'BiW', walk_length = 100) {
  # preparation
  A = A %*% FR
  V = pracma::nullspace(A)
  x_par = t(FR) %*% beta_slack(beta)
  # draw points
  P = volesti::Hpolytope(A = -V, b = as.vector(x_par))
  random_walk = list(
    walk = 'BiW',
    walk_length = 100
  )
  alpha = volesti::sample_points(P, n = npoints, random_walk = random_walk)
  soln = t(sapply(1:npoints, function(i) {
    beta_tight(FR %*% (x_par + V %*% alpha[, i]))
  }))
  
  return(soln)
}
