# Beta transformation

beta_slack = function(beta){
  beta_pos = sapply(beta, function(i){
    max(i, 0)
  })
  beta_neg = sapply(beta, function(i){
    -min(i, 0)
  })
  return( matrix(c(beta_pos, beta_neg), ncol = 1))
}

beta_tight = function(beta){
  len = length(beta)
  return( beta[1:(len/2)] - beta[(len/2 + 1):len] )
}