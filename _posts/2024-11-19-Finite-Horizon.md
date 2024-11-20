<a href="https://github.com/chiangnicholas/chiangnicholas.github.io/blob/main/_posts/2024-11-19-Finite-Horizon.md">View as markdown</a>
<a id="post-top"></a>
# A finite-horizon cake-eating problem

Consider a **finite** cake-eating problem of the following form:
```math
\displaylines{V(W) = \max_{W' \in [0,W]} log(W-W') + \beta V(W'),\\
\beta = 0.95, ~W \in [0.1,10]}
```

Here's a simple solution using a recursive procedure in R:

```R
#assign parameters
beta <- 0.95
t <- c(20)
s <- length(t)

#select state space [wlb wub]
wlb <- 0.1
wub <- 10

# Discretize state space
n <- 1001
inc <- (wub - wlb) / (n - 1)
wvec <- seq(wlb, wub, by = inc)

# Initialize the value function matrix
wm <- matrix(rep(wvec, n), ncol = n, byrow = TRUE)
cm <- wm - t(wm)
cm[cm == 0] <- 1e-6
um <- matrix(0, nrow = n, ncol = n)
um[cm < 0] <- -Inf
um[cm > 0] <- log(cm[cm > 0])
v <- matrix(0, nrow = n, ncol = s)

# Compute the value function
for (j in 1:s) {
  for (i in 2:t[j]) {
    vm <- um + beta * v[, j]
    p <- apply(vm, 2, which.max)
    v[, j] <- vm[cbind(p, seq_along(p))]
  }
}

# Define the true value function
a <- (log(1-beta) + (beta/(1-beta))* log(beta))/(1-beta)
b <- 1 / (1 - beta)
```

# Results
The plots for the value function with T=20 (LHS) and T=50 (RHS) are in the image below:

![image](https://github.com/user-attachments/assets/0c38c29f-c402-4f5b-bc77-4c721c31982a)
