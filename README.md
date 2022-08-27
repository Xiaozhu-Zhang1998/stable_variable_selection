# Stable Variable Selection

This repository provides a series of `R` functions to sample uniformly on the polytope formed by non-unique LASSO solutions. The dependent packages include `lpsolve` and `volesti`.


## Preliminaries
Consider the linear model $y = X\theta^* + w$, where the design matrix is $X\in\mathbb{R}^{n\times d}$, the true theta is $\theta^* \in\mathbb{R}^d$ with support $\mathcal{A}(\theta^*) \subseteq \lbrace 1,\dots, d \rbrace$, the noise is $w\in\mathbb{R}^n$, and the response is $y\in\mathbb{R}^n$.

We define the true theta set as

$$ 
\Theta^* := \lbrace \theta \in \mathbb{R}^d \mid \theta \in \arg\min_{ \vartheta \in \mathbb{R}^d} \\| \vartheta \\|_1 \text{ s.t. } X \vartheta = X\theta^*  \rbrace, 
$$

and define the predicted theta set as

$$
\hat{\Theta} := \lbrace \theta \in \mathbb{R}^d \mid \theta \in \arg\min_{\vartheta \in \mathbb{R}^d} \frac{1}{2n} \\| X \vartheta - y \\|_2^2 + \lambda \\| \vartheta \\|_1  \rbrace.
$$

If the set ![equation](https://latex.codecogs.com/svg.image?\Theta^*) (or the set ![equation](https://latex.codecogs.com/svg.image?\hat\Theta)) contains more than one element, then the set forms a polytope geometrically. 


## Find one solution

- True theta: The following function returns one ![equation](https://latex.codecogs.com/svg.image?\theta^*\in\Theta^*)

```
op_beta_star(X, beta)
```
where `X` denotes the design matrix, and `beta` denotes one vector such that `X` * `beta` achieves the desired prediction value. 


- Predicted theta: The popular function `glmnet` in the package `glmnet` is able to output one LASSO solution $\hat\theta$ given the design matrix $X$, the response $y$, the tuning parameter $\lambda$, and the parameter $\alpha = 1$.


## Find other solutions

- True thetas: 
The following function returns `n` uniformly distributed samples on the polytope $\Theta^*$:
```
sample_beta(type = "star", X = X, beta = beta, npoints = n)
```
where `X` is the design matrix, and `beta` is any one ![equation](https://latex.codecogs.com/svg.image?\theta^*\in\Theta^*)


- Predicted thetas:
The following function returns `n` uniformly distributed samples on the polytope $\Theta^*$:
```
sample_beta(type = "hat", X = X, beta = beta, y = y, lambda = lambda, npoints = n)
```
where `X` is the design matrix, `y` is the response, `lambda` is the tuning parameter, and `beta` is any one ![equation](https://latex.codecogs.com/svg.image?\hat\theta\in\hat\Theta).


### MCMC algorithm 
The MCMC procedure of sampling uniformly on the polytope ![equation](https://latex.codecogs.com/svg.image?\Theta) takes advantage of the function `volesti::sample_points()`, from which the function `sample_beta()` inherits several arguments regarding the MCMC setting such as `walk` and `walk_length`. The default arguments are
```
sample_beta(type, X, beta, y = NULL, lambda = NULL, npoints = 1000, walk = 'BiW', walk_length = 100, tol = 1e-2)
```
The argument `tol` controls the precision of determining the equicorrelation set, which depends on the precision of the given `beta`.


