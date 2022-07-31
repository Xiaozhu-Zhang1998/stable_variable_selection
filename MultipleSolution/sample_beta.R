# sample beta

source("beta_trans.R")
source("equi_corr.R")
source("find_orthant.R")
source("find_x_par.R")
source("find_interior.R")
source("target_draws.R")

sample_beta = function(type, X, beta, y = NULL, lambda = NULL, npoints = 1000, walk = 'BiW', walk_length = 100, tol1 = 1e-2, tol2 = 1e-7) {
  
  # save the beta
  beta_store = beta
  
  # check
  if(type == "hat") {
    if(is.null(y) | is.null(lambda)) {
      stop("y or lambda is missing!")
    }
    ind = equi_index(X, y, lambda, beta, tol = tol1)
  } else if(type == "star") {
    ind = gen_inv_index(X, beta, tol = tol1)
  } else {
    stop("The type is invalid!")
  }
  
  # check whether all elements of ind are 0
  if(length(ind) == 0) {
    return(matrix(1, nrow = npoints, ncol = 1) %*% t(as.matrix(beta_store)))
  }
  
  # subset according to equi-correlation / gen_inv index
  X = as.matrix(X[, ind])
  beta = beta[ind]
  
  # first find an interior point
  x_par = find_interior(X, beta, tol2) 

  # find the orthant the interior point
  orthant = find_orthant(x_par)
  
  # make "target draws"
  tryCatch({
    soln = target_draws(X, orthant, x_par, npoints, walk, walk_length = walk_length)
  },
  error = function(e) {
    soln <<- matrix(1, nrow = npoints, ncol = 1) %*% t(as.matrix(beta))
    print("The LASSO solution is unique.")
  }
  )
  
  # add back the fixed columns
  result = matrix(1, nrow = npoints, ncol = 1) %*% t(as.matrix(beta_store))
  result[, ind] = soln
  
  return(result)
}
