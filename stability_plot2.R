# user-specified rate
# a function that finds prob given rate
rate_2_prob = function(gamma, new_Theta, d) {
  prob = sapply(1:d, function(j) {
    mean(abs(new_Theta[,j]) > gamma)
  })
  return(prob)
}

# begin!
n = 500
d = 100
rho = -0.75
Sigma = matrix(0, nrow = d, ncol = d)
for(i in 1:d) {
  for(j in 1:d) {
    Sigma[i,j] = rho^abs(i-j)
  }
}
X = rmvnorm(n = n, mean = rep(0, d), sigma = Sigma)
X_new = X[,1] + rnorm(n, mean = 0, sd = 0.1)
X = cbind(X_new, X)
beta = c(0, rep(c(10, 0), 20), rep(0, 60))
epsilon = rnorm(n, sd = 0.1)
y = X %*% beta + epsilon

n = nrow(X)
d = ncol(X)

Lambda = seq(from = 0.01, to = 0.1, by = 0.005)
Store = list()
for(l in 15:19) {
  # page 1
  lambda = Lambda[l]
  model = glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, intercept = FALSE, standardize = TRUE)
  
  loss = sum((y - X %*% as.vector(model$beta))^2) / (2 * n) + lambda * sum(abs(model$beta))

  # page 2
  thres = loss + lambda * sum(abs(model$beta)) 
  # store
  N = 5000
  results = matrix(0, nrow = N, ncol = d)
  # start
  chain = 0
  pointer = 1
  xi = as.vector(model$beta)
  while(pointer <= N) {
    # store
    if(chain > 1000 & chain %% 300 == 0) {
      results[pointer, ] = xi
      pointer = pointer + 1
    }
    # randomly choose a coordinate
    j = sample(x = 1:d, size = 1)
    coord = xi[j]
    # first find p
    ubd = find_upper_bound(j, xi, type = 'p', tol2)
    p = binary_search(xi, coord, ubd, thres, tol)
    # then find q
    ubd = find_upper_bound(j, xi, type = 'q', tol2)
    q = binary_search(xi, coord, ubd, thres, tol)
    # rearrange p and q
    p_ = min(p, q)
    q_ = max(p, q)
    # uniformly draw
    coord = runif(1, min = p_, max = q_)
    xi[j] = coord
    # increment
    chain = chain + 1
  }
  new_loss = sapply(1:N, function(x){ loss_f(results[x,]) })
  new_Theta = data.frame(cbind(results, new_loss))
  # store
  Store[[l]] = new_Theta
  saveRDS(new_Theta, paste0(l, ".rds"))
}

# page 3
saveRDS(Store, "case.rds")

# user-specified
gamma = 0.2 * max(abs(as.vector(model$beta)))
Prob = matrix(0, nrow = length(Lambda), ncol = d)
for(l in seq_along(Lambda)) {
  #new_Theta = readRDS(paste0(l, ".rds"))
  new_Theta = Store[[l]]
  Prob[l,] = rate_2_prob(gamma, new_Theta, d)
}


ind = rep(as.vector(beta)!= 0, length(Lambda))
Prob = Prob %>%
  data.frame() %>%
  mutate(lambda = Lambda) %>%
  pivot_longer(!lambda, names_to = "feature", values_to = "prob") %>%
  mutate(ind = ind)

ggplot(Prob, aes(x = lambda, y = prob, type = feature, col = ind, lty = ind)) +
  geom_line(alpha = 0.3) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed")) +
  labs(col = "LASSO", title = "case 7") 

# importance
Prob %>%
  group_by(feature) %>%
  summarise(max_prob = max(prob),  avg_prob = mean(prob), med_prob = median(prob), min_prob = min(prob)) %>%
  arrange(desc(max_prob)) %>%
  View()


