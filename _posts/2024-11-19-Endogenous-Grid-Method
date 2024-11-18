<a href="https://github.com/chiangnicholas/chiangnicholas.github.io/blob/main/_posts/2024-11-17-Value-Function-Iteration.md">View as markdown</a>
<a id="post-top"></a>
# A deterministic growth model

Consider the following deterministic growth model:
```math
\displaylines{v(k) = \max_{k' \in [0,k^\alpha]} \{ log(k ^\alpha -k') + \beta v(k')\}\\
\beta = 0.9932, \alpha = 0.36.}
```

We can solve the model using the endogenous grid method.

We note that in the steady state, 
```math
k_{ss} = \frac{\alpha \beta}{1/(1-\alpha)}
```
since at the steady state, 
```math
1 = \beta f'(k_{ss})
```

Here's a simple solution in R:

```R
#assign parameters
beta <- 0.9932
alpha <- 0.36
#compute kss
kss <- (alpha * beta)^(1 / (1 - alpha))
#select state space [klb kub]
klb <- 0.001
kub <- 1.5 * kss
#define vectors and matrices
#discretize the state space
N <- 101
inc <- (kub - klb) / (N - 1)
kpgrid <- seq(klb, kub, by = inc) #construct grid for k’
cpgrid <- kpgrid ^ alpha #initial guess for the policy function
# Endogenous grid method
tolerance <- 1e-6
diff <- 100
cpgrid_1 <- 100
while (diff > tolerance) {
cgrid <- cpgrid/(alpha*beta) * kpgrid ^ (1-alpha) #from the euler equation
kgrid <- (cgrid + kpgrid)^(1/alpha) #from the resource constraint
for (i in 1:N) {
cpgrid[i] <- approx(x = kgrid, y = cgrid, xout = kpgrid[i])$y #updated function for c=h(k)
}
diff <- max(abs(cpgrid_1 - cpgrid)) #convergence is achieved if implied c’ converges
cpgrid_1 <- cpgrid
}
```

# Results
VFI gives the following results:

![image](https://github.com/user-attachments/assets/276e27ce-0d7e-4646-8bac-cb6184f7c899)
