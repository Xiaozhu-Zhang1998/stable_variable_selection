# facial reduction

facial_reduction = function(X, beta, tol2) {
  # define A and b
  d = ncol(X)
  A = rbind( cbind(X, -X), rep(1, 2 * d))
  b = rbind(X %*% beta, norm(matrix(beta), '1'))
  n = nrow(A)
  
  # find the y
  for(i in rev(1:n)) {
    # try yi = 1
    e = rep(0, n)
    e[i] = 1
    f.obj = rep(1, nrow(A))
    f.con = rbind(t(A), t(b), t(e))
    f.dir = c(rep(">=", ncol(A)), "==", "==")
    f.rhs = c(rep(0, ncol(A) + 1), 1)
    op = lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs)
    y = op$solution
    
    # find the FR
    z = t(A) %*% y
    sz = z > tol2
    if(sum(sz) != 0) {
      break
    }
    
    # try yi = -1
    e[i] = -1
    f.con = rbind(t(A), t(b), t(e))
    op = lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs)
    y = op$solution
    
    # find the FR
    z = t(A) %*% y
    sz = z > tol2
    if(sum(sz) != 0) {
      break
    }
  }
  
  FR = diag(rep(1, 2 * d))[, !sz]
  
  return(list(A = A, b = b, FR = FR))
}
