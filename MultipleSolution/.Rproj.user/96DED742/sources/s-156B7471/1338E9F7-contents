# find the orthant the draws

find_orthant = function(soln) {
  prop = apply(soln, 2, function(x) {
    mean(x > 0)
  })
  orthant = ifelse(prop > 0.5, 1, -1)
  return(orthant)
}