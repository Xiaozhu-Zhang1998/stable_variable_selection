# find a x_par that in the orthant

find_x_par = function(pilot_soln, orthant) {
  n = nrow(pilot_soln)
  while(1) {
    id = sample(x = n, size = 1)
    if(all.equal(sign(pilot_soln[id, ]), orthant) == TRUE) {
      return(pilot_soln[id, ])
    }
  }
}
