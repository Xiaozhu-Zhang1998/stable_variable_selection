# sample beta

source("beta_trans.R")
source("equi_corr.R")
source("find_orthant.R")
source("find_x_par.R")
source("pilot_draws.R")
source("target_draws.R")

sample_beta = function(type, X, beta, y = NULL, lambda = NULL, npoints = 1000, walk = 'BiW', walk_length = 100, tol = 1e-2) {

  # save the beta
  beta_store = beta
  
  # check
  if(type == "hat") {
    if(is.null(y) | is.null(lambda)) {
      stop("y or lambda is missing!")
    }
    ind = equi_index(X, y, lambda, beta, tol = tol)
  }
  else if(type == "star") {
    ind = gen_inv_index(X, beta, tol = tol)
  }
  else {
    stop("The type is invalid!")
  }
  
  # subset according to equi-correlation / gen_inv index
  X = X[, ind]
  beta = beta[ind]
  
  # first make several pilot draws
  tryCatch({
    pilot_soln = pilot_draws(X, beta, npoints = 1000)
    },
  error = function(e) {
    print("The LASSO solution is unique.")}
  )
  
  # find the orthant the draws
  orthant = find_orthant(pilot_soln)
  
  # find a x_par that in the orthant
  x_par = find_x_par(pilot_soln, orthant)
  
  # make "target draws"
  tryCatch({
    soln = target_draws(X, beta, orthant, x_par, npoints, walk, walk_length)
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

