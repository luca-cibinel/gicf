# An example with a banded covariance matrix
library(mvtnorm)

set.seed(1234)

p <- 10
n <- 500

# Create banded covariance matrix with three bands
band1 <- cbind(1:(p - 1), 2:p)
band2 <- cbind(1:(p - 2), 3:p)
band3 <- cbind(1:(p - 3), 4:p)
idxs <- rbind(band1, band2, band3)

Sigma <- matrix(0, p, p)
Sigma[idxs] <- 0.5
Sigma <- Sigma + t(Sigma)
diag(Sigma) <- 2

# Generate data
data <- rmvnorm(n, sigma = Sigma)
S <- cov(data) * (n - 1)/n

# Fix a value of lambda and compute k_max
lambda <- 0.07
k.max <- kappamax(S, lambda = lambda)

# Check that fit is diagonal
fit <- gicf(S = S, n = n, lambda = lambda, kappa = k.max)
image(fit$sigma != 0)


# Fix a value of kappa and compute l_max
kappa <- 1.15
l.max <- lambdamax(S, kappa = kappa)

# Check that fit is diagonal
fit <- gicf(S = S, n = n, lambda = l.max, kappa = kappa)
image(fit$sigma != 0)

# Repeat steps above, but with correct sparsity pattern specified
lambda <- 0.07
k.max <- kappamax(S, lambda = lambda, adj = Sigma)
fit <- gicf(S = S, n = n, lambda = lambda, kappa = k.max, adj = Sigma)
image(fit$sigma != 0)

kappa <- 1.15
l.max <- lambdamax(S, kappa = kappa, adj = Sigma)
fit <- gicf(S = S, n = n, lambda = l.max, kappa = kappa, adj = Sigma)
image(fit$sigma != 0)
