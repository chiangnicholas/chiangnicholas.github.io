<a href="https://github.com/chiangnicholas/chiangnicholas.github.io/blob/main/_posts/2024-11-17-Value-Function-Iteration.md">View as markdown</a>

# A basic cake-eating problem

Consider a basic cake-eating problem of the following form:
```math
V(W) = \max_{W' \in [0,W]} log(W-W') + \beta V(W'),
~\beta = 0.95, ~W \in [0.1,10]
```

This problem has a closed form solution to the value function and policy function. 

Here's a simple solution using discrete-state value function iteration in R:

```
#assign parameters
beta <- 0.95

#select state space [wlb wub]
wlb <- 0.1
wub <- 10

#define vectors and matrices 
#discretise the state space 
N <- 4001;
inc <- (wub-wlb)/(N-1)

wvec <- as.vector(seq(wlb, wub, by = inc))

#initial guess for value function 
V <- wvec; 
VNEW <- rep(0, N)

# Define matrices to record current utility
wM <- matrix(wvec, nrow = N, ncol = N, byrow = TRUE)
cM <- wM - wvec
cM[cM == 0] <- .Machine$double.eps
UM <- matrix(0, nrow = N, ncol = N)
UM[cM <= 0] <- -Inf
UM[cM > 0] <- log(cM[cM > 0])

# Beginning value function iteration
Vtol <- 10
epsilon <- 1e-6

while(Vtol > epsilon) {
  VpM <- V %*% matrix(1, ncol = N)
  VM <- UM + beta * VpM
  result <- apply(VM, 2, function(x) { c(max_value = max(x), max_index = which.max(x)) })
  VNEW <- result[1,]
  I <- result[2,]
  #Vtol <- max(abs(VNEW - V)) / (1 + max(abs(V)))
  diff_V <- abs(VNEW - V)
  diff_V <- diff_V[is.finite(diff_V)]
  Vtol <- max(diff_V) / (1 + max(abs(V)))
  V <- VNEW
}

wp <- wvec[I] # Optimal w'

# Define the true value function
a <- (log(1-beta) + (beta/(1-beta))* log(beta))/(1-beta)
b <- 1 / (1 - beta)
```

# Results
VFI gives the following results:

![image](https://github.com/user-attachments/assets/276e27ce-0d7e-4646-8bac-cb6184f7c899)