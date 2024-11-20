<a href="https://github.com/chiangnicholas/chiangnicholas.github.io/blob/main/_posts/2024-11-17-Value-Function-Iteration.md">View as markdown</a>
<a id="post-top"></a>
# A basic cake-eating problem

Consider a basic cake-eating problem of the following form:
```math
\displaylines{V(W) = \max_{W' \in [0,W]} log(W-W') + \beta V(W'),\\
\beta = 0.95, ~W \in [0.1,10]}
```

This problem has a [closed form solution](#analytical-solution) to the value function and policy function. 

Here's a simple solution using discrete-state value function iteration in R:

```R
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

#Define matrices to record current utility
wM <- matrix(wvec, nrow = N, ncol = N, byrow = TRUE)
cM <- wM - wvec
cM[cM == 0] <- .Machine$double.eps
UM <- matrix(0, nrow = N, ncol = N)
UM[cM <= 0] <- -Inf
UM[cM > 0] <- log(cM[cM > 0])

#Beginning value function iteration
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

#Define the true value function
a <- (log(1-beta) + (beta/(1-beta))* log(beta))/(1-beta)
b <- 1 / (1 - beta)
```

# Results
VFI gives the following results:

![image](https://github.com/user-attachments/assets/276e27ce-0d7e-4646-8bac-cb6184f7c899)

# Analytical solution

We conjecture that the solution to $V(W)$ takes the form $A + B log(W)$.
```math
V(W) = \max_{W' \in [0,W]} \{ {log(W-W')+\beta (A + B(log (W'))} \}
```
FOC for the maximisation problem is given by:
```math
\displaylines{\frac{-1}{W-W'}+ \frac{\beta B}{W'}=0 \\
W' = \frac{\beta B}{1+ \beta B} W}
```
Substituting the optimal W' into V(W), we can drop the max operator:
```math
\displaylines{A + B log(W) &= log(W-\frac{\beta B}{1+ \beta B} W) + \beta(A + B log(\frac{\beta B W}{1 + \beta B})) \\
        &= log(W(\frac{1+\beta B - \beta B}{1 + \beta B}) + \beta(A + B log(\frac{\beta B W}{1 + \beta B})) \\
        &= log(\frac{W}{1+\beta B})+\beta(A+Blog(\frac{\beta B W}{1 +\beta B})) \\
        &= (1 + \beta B) log(W) - log(1 + \beta B) + \beta A + \beta B log(\frac{\beta B}{1 + \beta B})}
```
By comparing coefficients, $B = 1 + \beta B$, and hence $B = \frac{1}{1-\beta}$. 

Substituting B back into the FOC, we get the optimal consumption path:
```math
\displaylines{W' = \frac{\beta(\frac{1}{1-\beta})}{1+\beta \frac{1}{1-\beta}}W = \beta W \\
        C = W - W' = (1-\beta) W}
```
    
Substituting B back into the value function, we get:
```math
\displaylines{A + B log(W) &= (1 + \beta B) log(W) - log(1 + \beta \frac{1}{1-\beta}) + \beta A + \beta \frac{1}{1-\beta} log(\frac{\beta \frac{1}{1-\beta}}{1 + \beta \frac{1}{1-\beta}}) \\
        & = B log(W) - log(\frac{1}{1-\beta}) + \beta A + \frac{\beta}{1-\beta}log(\beta)   }
```
    
By comparing coefficients for A, we get:
```math
\displaylines{A - \beta A = - log(\frac{1}{1-\beta}) + \frac{\beta}{1-\beta}log(\beta) \\
        A - \beta A = log(1-\beta) + \frac{\beta}{1-\beta}log(\beta) \\
        A = \frac{log(1-\beta) + \frac{\beta}{1-\beta}log(\beta)}{1-\beta} }
```

<p align="right">(<a href="#post-top">back to top</a>)</p>
