# facial reduction

facial_reduction = function(X, beta, tol2) {
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
  
  flag = FALSE
  while(TRUE) {
    # try-catch
    tryCatch({
      P = volesti::Hpolytope(A = A_, b = as.vector(b_) + abs(rnorm(2 * d, sd = 1e-10)))
      random_walk = list(
        walk = 'BiW',
        walk_length = 1
      )
      alpha = volesti::sample_points(P, n = 100, random_walk = random_walk)
    },
    error = function(e) {
      flag <<- TRUE
    })
    # flag
    if(!flag) {
      break
    }
  }
  
  y = y_par + V %*% alpha[, 100]
  
  # find the FR
  z = t(A) %*% y
  sz = z > tol2
  FR = diag(rep(1, 2 * d))[, !sz]
  
  return(list(A = A, b = b, FR = FR))
}
