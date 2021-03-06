# find one beta_star

op_beta_star = function(X, beta) {
  p = length(beta)
  f.obj = rep(1, 2 * p)
  f.con = cbind(X, -X)
  f.dir = rep("=", nrow(X))
  f.rhs = X %*% beta
  op = lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs)
  slack = op$solution
  return(slack[1:p] - slack[(p + 1):(2 * p)])
}
