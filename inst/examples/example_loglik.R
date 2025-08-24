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

# Fix a value of lambda and kappa
lambda <- 0.07
kappa <- 0.5

# Gaussian loglikelihood
print("Gaussian loglikelihood:")
print(gcgmloglik(Sigma, S, n))

# Penalised Gaussian loglikelihood
print("Penalised Gaussian loglikelihood:")
print(gcgmloglik(Sigma, S, n, lambda, kappa))
