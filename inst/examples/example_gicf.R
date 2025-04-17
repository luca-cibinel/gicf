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

# Fit a path of estimates
lambdas <- seq(0, 0.15, 0.01)
fit <- gicf(data, lambda = lambdas, kappa = 0.1)

# Explore one particular estimate
onefit <- fit[[5]]
image(onefit$sigma != 0)

# Redo the fit, but this time fix the correct sparsity pattern
fit2 <- gicf(data, lambda = lambdas, kappa = 0.1, adj = Sigma)

onefit2 <- fit2[[5]]
image(onefit2$sigma != 0)
