# find the orthant the draws

find_orthant = function(x_par) {
  # prop = apply(soln, 2, function(x) {
  #   mean(x > 0)
  # })
  # orthant = ifelse(prop > 0.5, 1, -1)
  orthant = sign(x_par)
  orthant[orthant == 0] = 1
  return(orthant)
}