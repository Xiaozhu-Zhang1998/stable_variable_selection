# different cases

# case 1: large n / beta-type 1 / small rho ----
n = 400
d = 50
rho = -0.1
Sigma = matrix(0, nrow = d, ncol = d)
for(i in 1:d) {
  for(j in 1:d) {
    Sigma[i,j] = rho^abs(i-j)
  }
}
X = rmvnorm(n = n, mean = rep(0, d), sigma = Sigma)
beta = c(rep(10, 10), rep(0, 40))
beta[2] = 0
#beta = c(rep(c(10, 0), 20), rep(0, 60))
epsilon = rnorm(n, sd = 0.1)
y = X %*% beta + epsilon
lambda = 2 * max(abs(t(X) %*% epsilon)) / n
model = glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE,
               maxit = 1e8, thresh = 1e-20)
as.vector(model$beta)

loss = sum((y - X %*% model$beta)^2) / (2 * n) + lambda * sum(abs(model$beta))
sum((y - X %*% beta)^2) / (2 * n) + lambda * sum(abs(beta)) - loss

s = sum(model$beta != 0)
ind = (1:d)[as.vector(model$beta != 0)]
c_min = min(svd(t(X[, ind]) %*% X[, ind] / n)$d) 



# case 2: large n / beta-type 1 / large rho ----
n = 400
d = 50
rho = -0.99
Sigma = matrix(0, nrow = d, ncol = d)
for(i in 1:d) {
  for(j in 1:d) {
    Sigma[i,j] = rho^abs(i-j)
  }
}
X = rmvnorm(n = n, mean = rep(0, d), sigma = Sigma)
beta = c(rep(10, 10), rep(0, 40))
beta[2] = 0
#beta = c(rep(c(10, 0), 20), rep(0, 60))
epsilon = rnorm(n, sd = 0.1)
y = X %*% beta + epsilon
lambda = 2 * max(abs(t(X) %*% epsilon)) / n
model = glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE,
               maxit = 1e8, thresh = 1e-20)
as.vector(model$beta)

loss = sum((y - X %*% model$beta)^2) / (2 * n) + lambda * sum(abs(model$beta))
sum((y - X %*% beta)^2) / (2 * n) + lambda * sum(abs(beta)) - loss

s = sum(model$beta != 0)
ind = (1:d)[as.vector(model$beta != 0)]
c_min = min(svd(t(X[, ind]) %*% X[, ind] / n)$d) 



# case 3: large n / beta-type 2 / small rho ----
n = 400
d = 50
rho = -0.1
Sigma = matrix(0, nrow = d, ncol = d)
for(i in 1:d) {
  for(j in 1:d) {
    Sigma[i,j] = rho^abs(i-j)
  }
}
X = rmvnorm(n = n, mean = rep(0, d), sigma = Sigma)
beta = c(rep(c(10, 0), 5), rep(0, 40))
epsilon = rnorm(n, sd = 0.1)
y = X %*% beta + epsilon
lambda = 2 * max(abs(t(X) %*% epsilon)) / n
model = glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE,
               maxit = 1e8, thresh = 1e-20)
as.vector(model$beta)

loss = sum((y - X %*% model$beta)^2) / (2 * n) + lambda * sum(abs(model$beta))
sum((y - X %*% beta)^2) / (2 * n) + lambda * sum(abs(beta)) - loss

s = sum(model$beta != 0)
ind = (1:d)[as.vector(model$beta != 0)]
c_min = min(svd(t(X[, ind]) %*% X[, ind] / n)$d) 



# case 4: large n / beta-type 2 / large rho ----
n = 400
d = 50
rho = -0.99
Sigma = matrix(0, nrow = d, ncol = d)
for(i in 1:d) {
  for(j in 1:d) {
    Sigma[i,j] = rho^abs(i-j)
  }
}
X = rmvnorm(n = n, mean = rep(0, d), sigma = Sigma)
beta = c(rep(c(10, 0), 5), rep(0, 40))
epsilon = rnorm(n, sd = 0.1)
y = X %*% beta + epsilon
lambda = 2 * max(abs(t(X) %*% epsilon)) / n
model = glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE,
               maxit = 1e8, thresh = 1e-20)
as.vector(model$beta)

loss = sum((y - X %*% model$beta)^2) / (2 * n) + lambda * sum(abs(model$beta))
sum((y - X %*% beta)^2) / (2 * n) + lambda * sum(abs(beta)) - loss

s = sum(model$beta != 0)
ind = (1:d)[as.vector(model$beta != 0)]
c_min = min(svd(t(X[, ind]) %*% X[, ind] / n)$d) 



# case 5: small n / beta-type 1 / small rho ----
n = 50
d = 50
rho = -0.1
Sigma = matrix(0, nrow = d, ncol = d)
for(i in 1:d) {
  for(j in 1:d) {
    Sigma[i,j] = rho^abs(i-j)
  }
}
X = rmvnorm(n = n, mean = rep(0, d), sigma = Sigma)
beta = c(rep(10, 10), rep(0, 40))
beta[2] = 0
#beta = c(rep(c(10, 0), 20), rep(0, 60))
epsilon = rnorm(n, sd = 0.1)
y = X %*% beta + epsilon
lambda = 2 * max(abs(t(X) %*% epsilon)) / n
model = glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE,
               maxit = 1e8, thresh = 1e-20)
as.vector(model$beta)

loss = sum((y - X %*% model$beta)^2) / (2 * n) + lambda * sum(abs(model$beta))
sum((y - X %*% beta)^2) / (2 * n) + lambda * sum(abs(beta)) - loss

s = sum(model$beta != 0)
ind = (1:d)[as.vector(model$beta != 0)]
c_min = min(svd(t(X[, ind]) %*% X[, ind] / n)$d)  



# case 6: small n / beta-type 1 / large rho ----
n = 50
d = 50
rho = -0.99
Sigma = matrix(0, nrow = d, ncol = d)
for(i in 1:d) {
  for(j in 1:d) {
    Sigma[i,j] = rho^abs(i-j)
  }
}
X = rmvnorm(n = n, mean = rep(0, d), sigma = Sigma)
beta = c(rep(10, 10), rep(0, 40))
beta[2] = 0
#beta = c(rep(c(10, 0), 20), rep(0, 60))
epsilon = rnorm(n, sd = 0.1)
y = X %*% beta + epsilon
lambda = 2 * max(abs(t(X) %*% epsilon)) / n
model = glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE,
               maxit = 1e8, thresh = 1e-20)
as.vector(model$beta)

loss = sum((y - X %*% model$beta)^2) / (2 * n) + lambda * sum(abs(model$beta))
sum((y - X %*% beta)^2) / (2 * n) + lambda * sum(abs(beta)) - loss

s = sum(model$beta != 0)
ind = (1:d)[as.vector(model$beta != 0)]
c_min = min(svd(t(X[, ind]) %*% X[, ind] / n)$d) 



# case 7: small n / beta-type 2 / small rho ----
n = 50
d = 50
rho = -0.1
Sigma = matrix(0, nrow = d, ncol = d)
for(i in 1:d) {
  for(j in 1:d) {
    Sigma[i,j] = rho^abs(i-j)
  }
}
X = rmvnorm(n = n, mean = rep(0, d), sigma = Sigma)
beta = c(rep(c(10, 0), 5), rep(0, 40))
epsilon = rnorm(n, sd = 0.1)
y = X %*% beta + epsilon
lambda = 2 * max(abs(t(X) %*% epsilon)) / n
model = glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE,
               maxit = 1e8, thresh = 1e-20)
as.vector(model$beta)

loss = sum((y - X %*% model$beta)^2) / (2 * n) + lambda * sum(abs(model$beta))
sum((y - X %*% beta)^2) / (2 * n) + lambda * sum(abs(beta)) - loss

s = sum(model$beta != 0)
ind = (1:d)[as.vector(model$beta != 0)]
c_min = min(svd(t(X[, ind]) %*% X[, ind] / n)$d) 



# case 8: small n / beta-type 2 / large rho ----
n = 50
d = 50
rho = -0.99
Sigma = matrix(0, nrow = d, ncol = d)
for(i in 1:d) {
  for(j in 1:d) {
    Sigma[i,j] = rho^abs(i-j)
  }
}
X = rmvnorm(n = n, mean = rep(0, d), sigma = Sigma)
beta = c(rep(c(10, 0), 5), rep(0, 40))
epsilon = rnorm(n, sd = 0.1)
y = X %*% beta + epsilon
lambda = 2 * max(abs(t(X) %*% epsilon)) / n
model = glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE,
               maxit = 1e8, thresh = 1e-20)
as.vector(model$beta)

loss = sum((y - X %*% model$beta)^2) / (2 * n) + lambda * sum(abs(model$beta))
sum((y - X %*% beta)^2) / (2 * n) + lambda * sum(abs(beta)) - loss

s = sum(model$beta != 0)
ind = (1:d)[as.vector(model$beta != 0)]
c_min = min(svd(t(X[, ind]) %*% X[, ind] / n)$d) 



# case 9: large n / one column highly correlated ----
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
lambda = 2 * max(abs(t(X) %*% epsilon)) / n
model = glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE,
               maxit = 1e8, thresh = 1e-20)
as.vector(model$beta)

loss = sum((y - X %*% model$beta)^2) / (2 * n) + lambda * sum(abs(model$beta))
sum((y - X %*% beta)^2) / (2 * n) + lambda * sum(abs(beta)) - loss

s = sum(model$beta != 0)
ind = (1:d)[as.vector(model$beta != 0)]
c_min = min(svd(t(X[, ind]) %*% X[, ind] / n)$d) 
d = d + 1


# case 10: small n / one column highly correlated ----
n = 100
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
lambda = 2 * max(abs(t(X) %*% epsilon)) / n
model = glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE,
               maxit = 1e8, thresh = 1e-20)
as.vector(model$beta)

loss = sum((y - X %*% model$beta)^2) / (2 * n) + lambda * sum(abs(model$beta))
sum((y - X %*% beta)^2) / (2 * n) + lambda * sum(abs(beta)) - loss

s = sum(model$beta != 0)
ind = (1:d)[as.vector(model$beta != 0)]
c_min = min(svd(t(X[, ind]) %*% X[, ind] / n)$d) 
d = d + 1


# case 11: small n / complicated beta-type 1 / large rho ----
n = 50
d = 50
rho = -0.99
Sigma = matrix(0, nrow = d, ncol = d)
for(i in 1:d) {
  for(j in 1:d) {
    Sigma[i,j] = rho^abs(i-j)
  }
}
X = rmvnorm(n = n, mean = rep(0, d), sigma = Sigma)
beta = c(rep(c(10, 1), 5), rep(c(2, 1), 5), rep(0, 30))
epsilon = rnorm(n, sd = 0.1)
y = X %*% beta + epsilon
lambda = 2 * max(abs(t(X) %*% epsilon)) / n
model = glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, intercept = FALSE, standardize = FALSE,
               maxit = 1e8, thresh = 1e-20)
as.vector(model$beta)

loss = sum((y - X %*% model$beta)^2) / (2 * n) + lambda * sum(abs(model$beta))
sum((y - X %*% beta)^2) / (2 * n) + lambda * sum(abs(beta)) - loss

s = sum(model$beta != 0)
ind = (1:d)[as.vector(model$beta != 0)]
c_min = min(svd(t(X[, ind]) %*% X[, ind] / n)$d) 



# real dataset
library(hdi)
data(riboflavin)

X = riboflavin$x
y = riboflavin$y
n = nrow(X)
d = ncol(X)
lambda = 0.01
model = glmnet(X, y, family = "gaussian", alpha = 1, lambda = lambda, intercept = FALSE, standardize = TRUE)

loss = sum((y - X %*% as.vector(model$beta))^2) / (2 * n) + lambda * sum(abs(model$beta))
s = sum(model$beta != 0)
ind = (1:d)[as.vector(model$beta != 0)]
c_min = min(svd(t(X[, ind]) %*% X[, ind] / n)$d) 
