# HEADER ====
rm(list = ls()) # clear environment
source("../gicf/gicf_source.R")

set.seed(1234)

d.seq <- c(1, 3, 5) # Sequence of 10*(true graph densities)
n.seq <- c(45, 75, 100, 250, 500, 1000) # Sequence of dataset size
p <- 50 # Model size

N.sim <- 10 # Number of simulations
N.folds <- 5 # Number of CV folds
N.kappa <- 30 # Number of values of \kappa to sample
N.d <- length(d.seq)
N.n <- length(n.seq)

zero <- 1e-4 # Tolerance: absolute values below this threshold are treated as 0

# Parameters to access simulated data (DO NOT CHANGE)
model <- "RB" # Model name: Random entries with Banded structure
n.key <- 2000 # Key to access simulated data

# METRICS ====
metrics.names <- c(
  "kappa", # optimal \kappa
  "condnum", # condition number
  "EL", # Entropy Loss
  "RMSE" # Root Mean Square Error
)

metrics <- array(zero, c(length(metrics.names), N.d, N.sim, N.n, 2)) # [metric name, density, simulation, n, MLE/RIDGE]
dimnames(metrics)[[1]] <- metrics.names
dimnames(metrics)[[5]] <- c("MLE", "RIDGE")

true.condnum <- rep(0, N.d) # Store true condition number

condnum <- function(A){
  eig <- eigen(A)$values
  
  return(max(eig)/min(eig))
}

entropy.loss <- function(A, Sigma, detSigma, Theta){
  return(sum(A*Theta) + log(detSigma/det(A)) - nrow(A))
}

rmse <- function(A, Sigma){
  return(norm(A - Sigma, "F")/(nrow(A)^2))
}

# MODEL SELECTION ====
model.selection.cv <- function(y, # data
                                 N.folds, # num. of CV folds
                                 N.k, # num. of candidates
                                 k.max, # max value of \kappa to be explored
                                 adj = 1 - diag(ncol(y)) # prespecified adj. matrix
                                 ){
  n <- nrow(y)
  p <- ncol(y)
  
  S <- cov(y) * (n - 1)/n
  
  folded.data <- fold(data.frame(y), k = N.folds)
  
  vals <- rep(0, N.k)
  
  if(n > p)
    k.seq <- seq(0, k.max, length.out = N.k)
  else
    k.seq <- seq(0, k.max, length.out = N.k + 1)[-1]
  
  for(fold in 1:N.folds){
    train <- folded.data[folded.data$.folds != fold, ]
    test <- folded.data[folded.data$.folds == fold, ]
    
    train <- unname(as.matrix(
      train[, names(train) != ".folds"]
    ))
    
    test <- unname(as.matrix(
      test[, names(test) != ".folds"]
    ))
    
    n.train <- nrow(train)
    n.test <- nrow(test)
    
    S.train <- cov(train) * (n.train - 1)/n.train
    S.test <- cov(test) * (n.test - 1)/n.test
    
    for(K in 1:N.k){
      D <- diag(k.seq[K], p)
      fit <- gicf.run(S = S.train + D, n = n.train, lambda = 0,
                      adj = adj)$sigma
      
      vals[K] <- vals[K] + gicf.likelihood(fit, S.test, n.test, 0)
    }
  }
  
  return(k.seq[which.max(vals)])
}

# SIMULATION ====
for(d in 1:N.d){ # For each graph density
  density <- d.seq[d]/10
  
  print("")
  print("")
  print(paste("DENSITY:", density))
  
  Sigma <- as.matrix( # Read true Sigma from simulated data
    read.table(
      paste0("data/sigma_mod_", model, "_d_", d.seq[d], "_p_", p, "_n_", n.key, ".dat")
    )
  )
  
  Theta <- solve(Sigma) # \Sigma^{-1}
  detsigma <- det(Sigma) # det(\Sigma)
  
  adj <- matrix(0, p, p) # Compute adjacency matrix
  adj[abs(Sigma) > zero] <- 1
  diag(adj) <- 0
  
  true.condnum[d] <- condnum(Sigma) # Store true condition number
  
  data <- as.matrix( # Read simulated data
    read.table(
      paste0("data/simul_mod_", model, "_d_", d.seq[d], "_p_", p, "_n_", n.key, ".dat")
    )
  )
  
  for(s in 1:N.sim){ # Repeat for each simulation:
    print(paste("Simulation:", s))
    
    for(m in 1:N.n){ # For each desired dataset size:
      cat(paste("n =", n.seq[m], "-"))
      n.loc <- n.seq[m] # Local dataset size
      
      y <- data[(s - 1)*n.seq[N.n] + 1:n.loc, ] # Extract local dataset
      
      S <- cov(y) * (n.loc - 1)/n.loc # Compute local sample covariance matrix
      
      # Compute MLE (if possible)
      if(n.loc > p){
        sigma.mle <- gicf.run(S = S, n = n.loc, adj = adj)$sigma # MLE
        
        metrics["kappa", d, s, m, "MLE"] <- 0
        metrics["condnum", d, s, m, "MLE"] <- condnum(sigma.mle)
        metrics["EL", d, s, m, "MLE"] <- entropy.loss(sigma.mle, Sigma, detsigma, Theta)
        metrics["RMSE", d, s, m, "MLE"] <- rmse(sigma.mle, Sigma)
      }
      
      # Compute shrinked estimate
      kappa.sel <- model.selection.cv(y, N.folds, N.kappa, k.max = 3, adj = adj)
      sigma.ridge <- gicf.run(S = S + diag(kappa.sel, p), n = n.loc, adj = adj)$sigma
      
      metrics["kappa", d, s, m, "RIDGE"] <- kappa.sel
      metrics["condnum", d, s, m, "RIDGE"] <- condnum(sigma.ridge)
      metrics["EL", d, s, m, "RIDGE"] <- entropy.loss(sigma.ridge, Sigma, detsigma, Theta)
      metrics["RMSE", d, s, m, "RIDGE"] <- rmse(sigma.ridge, Sigma)
    }
    
    print("")
  }
}

# OUTPUT ====
library(latex2exp)

or.par <- par(mar = c(5, 4, 4, 2) + 0.5)
clrs <- palette.colors(palette = "Okabe-Ito", n = N.d)

# Titles of the plots
titles <- c(
  kappa = TeX("$\\hat{\\kappa}$"),
  condnum = "Condition Number",
  EL = "Entropy Loss",
  RMSE = "RMSE"
)

# Y-labels of the plots
ylabs <- c(
  kappa = TeX("$\\hat{\\kappa}$"),
  condnum = TeX('cond(\\hat{\\Sigma})'),
  EL = "Entropy Loss",
  RMSE = "RMSE"
)

compute.lim <- metrics.names # For which metrics should the plot limits be computed?
plot.mle <- setdiff(metrics.names, c("kappa")) # For which metrics should the mle-related
                                                  # metrics be plotted?

for(nm in metrics.names){
  if(nm %in% compute.lim){
    tempM <- 0
    tempm <- Inf
    
    for(d in 1:N.d){
      temp1 <- colMeans(metrics[nm, d, , , "MLE"])[-1]
      temp2 <- colMeans(metrics[nm, d, , , "RIDGE"])
      
      tempM <- max(tempM, max(temp1, temp2))
      tempm <- min(tempm, min(temp1, temp2))
    }
  }else{
    tempM <- 1
    tempm <- 0
  }
  
  plot(n.seq, colMeans(metrics[nm, 1, , , "RIDGE"]),
       ylim = c(tempm, tempM),
       col = clrs[1], type = "l", log = "x",
       ylab = ylabs[nm],
       xlab = "n",
       main = titles[nm],
       lwd = 1)
  points(n.seq, colMeans(metrics[nm, 1, , , "RIDGE"]),
         col = clrs[1], pch = 20 + 1, bg = clrs[1] )
  
  if(nm %in% plot.mle){
    lines(n.seq[-1], colMeans(metrics[nm, 1, , , "MLE"])[-1],
          col = clrs[1], type = "l", lty = 2)
    points(n.seq[-1], colMeans(metrics[nm, 1, , , "MLE"])[-1],
           col = clrs[1], pch = 20 + 1 )
  }
  
  for(d in 2:N.d){
    lines(n.seq, colMeans(metrics[nm, d, , , "RIDGE"]),
          col = clrs[d], type = "l", lwd = 1)
    points(n.seq, colMeans(metrics[nm, d, , , "RIDGE"]),
           col = clrs[d], pch = 20 + d, bg = clrs[d] )
    
    if(nm %in% plot.mle){
      lines(n.seq[-1], colMeans(metrics[nm, d, , , "MLE"])[-1],
            col = clrs[d], type = "l", lty = 2)
      points(n.seq[-1], colMeans(metrics[nm, d, , , "MLE"])[-1],
             col = clrs[d], pch = 20 + d )
    }
  }
  
  legend("topright", col = clrs, pch = (20 + d.seq), lty = 1,
         legend = paste0(d.seq*10, "%"), title = "Density",
         pt.bg = clrs, lwd = 2)
}

par <- or.par