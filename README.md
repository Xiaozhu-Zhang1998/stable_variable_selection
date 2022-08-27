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

If the set $\Theta^*$ (or the set $\hat\Theta$) contains more than one element, then the set forms a polytope geometrically.


## Find one solution
- The function `op_beta_star` inputs the design matrix `X` and one `beta` such that `X``beta` achieves the desired prediction value. It outputs one corresponding $\theta \in \Theta^*$.

- 
