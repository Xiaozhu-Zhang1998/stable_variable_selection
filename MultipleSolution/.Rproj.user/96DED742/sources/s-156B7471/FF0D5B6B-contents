# sample beta

library(walkr)
library(pracma)

source("beta_trans.R")
source("equi_corr.R")


sample_beta = function(type, X, y, lambda, beta, points, method = 'dikin', ntry = 2000, ...) {
  # check the nullspace
  Nullspace = nullspace(X)
  if(is.null(Nullspace)) {
    stop("The design matrix is of full column rank")
  }
  
  # find the equicorrelation set
  p = length(beta)
  beta = matrix(as.vector(beta), ncol = 1)
  if(type == "hat") {
    ind = equi_index(X, y, lambda, beta)
  }
  else {
    ind = gen_inv_index(X, beta)
  }
  
  X = X[, ind]
  beta = beta[ind]
  
  # sample solutions on equi-set
  l1norm = sum(abs(beta))
  A = cbind(X, -X)
  Rowspace = t(orth(t(A)))
  b = as.vector(Rowspace %*% beta_slack(beta)) / l1norm
  
  for(j in 1:ntry){
    error_occur = FALSE
    tryCatch({
      sampled_points <- walkr(A = Rowspace, b = b + rnorm(n = length(b), sd = 0.0001), points = points, 
                              method = "dikin", ...)
    },
    error = function(e) {error_occur <<- TRUE}
    )
    if(error_occur == FALSE) {
      break
    }
  }
  
  # deal with samples
  if(j == ntry & error_occur == TRUE) {
    print("Either the solution is unique, or the sampler is broken.")
  }
  else {
    beta = t(apply(sampled_points, 2, beta_tight)) * l1norm
    # convert beta to the full length
    result = matrix(0, nrow = points, ncol = p)
    result[, ind] = beta 
    return(result)
  }
}
