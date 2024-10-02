# HEADER ====
rm(list = ls()) # clear environment

library(covglasso)
library(mvtnorm)
source("../gicf/gicf_for_time_sim/gicf_source_time.R")
set.seed(1234)

d.seq <- c(0.15, 0.65) # Sequence of 10*(true graph densities), correspond to 3 and 10 bands
p.seq <- c(25, 50, 100, 200, 300) # Sequence of models size
n <- 2000 # Dataset size

p.max <- max(p.seq) # Maximum model size

N.sim <- 10 # Number of simulations
N.d <- length(d.seq)
N.p <- length(p.seq)

zero <- 1e-4 # Tolerance: absolute values below this threshold are treated as 0

# Parameters to access simulated data (DO NOT CHANGE)
model <- "RB" # Model name: Random entries with Banded structure
n.key <- 2000 # Key to access simulated data

# METRICS ====
metrics.names <- c(
  "time", # computational time
  "loglik", # Optimised loglikelihood
  "condnum", # condition number
  "EL", # Entropy Loss
  "RMSE", # Root Mean Square Error
  "bRMSE", # Root Mean Square Error (between the two estimates)
  "iter" # Number of iterations
)

metrics <- array(zero, c(length(metrics.names), N.d, N.sim, N.p, 2)) # [metric name, density, simulation, n, HYB/CVGL]
dimnames(metrics)[[1]] <- metrics.names
dimnames(metrics)[[5]] <- c("GICF", "CVGL")

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

# SIMULATION ====
for(d in 1:N.d){ # For each graph density
  density <- d.seq[d]/10
  
  print("")
  print("")
  print(paste("DENSITY:", density))
  
  Sigma <- as.matrix( # Read true Sigma from simulated data
    read.table(
      paste0("data/sigma_mod_", model, "_d_", d.seq[d], "_p_", p.max, "_n_", n.key, ".dat")
    )
  )
  
  Theta <- solve(Sigma) # \Sigma^{-1}
  detsigma <- det(Sigma) # det(\Sigma)
  
  adj <- matrix(0, p.max, p.max) # Compute adjacency matrix
  adj[abs(Sigma) > zero] <- 1
  diag(adj) <- 0
  
  true.condnum[d] <- condnum(Sigma) # Store true condition number
  
  for(s in 1:N.sim){ # Repeat for each simulation:
    print(paste("Simulation:", s))
    
    for(q in 1:N.p){ # For each desired model size:
      cat(paste("p =", p.seq[q], "-"))
      p.loc <- p.seq[q] # Local model size
      
      Sigma.loc <- Sigma[1:p.loc, 1:p.loc] # Local Sigma
      adj.loc <- adj[1:p.loc, 1:p.loc] # Local adjacency
      
      y <- rmvnorm(n, sigma = Sigma.loc)
      
      S <- cov(y) * (n - 1)/n # Compute local sample covariance matrix
      
      # Compute MLE via GICF
      tic <- proc.time()
      fit.hyb <- gicf.run(S = S, n = n, adj = adj.loc)
      toc <- proc.time() - tic
      toc <- unname(toc["elapsed"])
      sigma.hyb <- fit.hyb$sigma
      
      metrics["time", d, s, q, "GICF"] <- toc
      metrics["condnum", d, s, q, "GICF"] <- condnum(sigma.hyb)
      metrics["EL", d, s, q, "GICF"] <- entropy.loss(sigma.hyb, Sigma.loc, det(Sigma.loc), solve(Sigma.loc))
      metrics["RMSE", d, s, q, "GICF"] <- rmse(A = sigma.hyb, Sigma = Sigma.loc)
      metrics["loglik", d, s, q, "GICF"] <- gicf.likelihood(sigma.hyb, S, n, 0)
      metrics["iter", d, s, q, "GICF"] <- fit.hyb$it
      
      
      # Compute MLE via covglasso algorithm
      L <- matrix(0, p.loc, p.loc)
      L[adj.loc == 0] <- 1e6
      diag(L) <- 0
      temp <- array(0, c(p.loc, p.loc, 1))
      temp[,,1] <- L
      L <- temp
      
      tic <- proc.time()
      fit.cvgl <- covglasso(S = S, n = n, lambda = L, path = T)$path[[1]] # path needed to extract n. of iterations
      toc <- proc.time() - tic
      toc <- unname(toc["elapsed"])
      
      sigma.cvgl <- fit.cvgl$sigma
      
      metrics["time", d, s, q, "CVGL"] <- toc
      metrics["condnum", d, s, q, "CVGL"] <- condnum(sigma.cvgl)
      metrics["EL", d, s, q, "CVGL"] <- entropy.loss(sigma.cvgl, Sigma.loc, det(Sigma.loc), solve(Sigma.loc))
      metrics["RMSE", d, s, q, "CVGL"] <- rmse(A = sigma.cvgl, Sigma = Sigma.loc)
      metrics["loglik", d, s, q, "CVGL"] <- gicf.likelihood(sigma.cvgl, S, n, 0)
      metrics["iter", d, s, q, "CVGL"] <- fit.cvgl$it
      
      # Compute error between estimates
      metrics["bRMSE", d, s, q, "GICF"] <- rmse(sigma.cvgl, sigma.hyb)
      
    }
    
    print("")
  }
}

# OUTPUT ====

## Plotting ====
library(latex2exp)

or.par <- par(mar = c(5, 4, 4, 2) + 0.5)
clrs <- palette.colors(palette = "Okabe-Ito", n = N.d)

# Plot titles
titles <- c(
  time = "Computational time",
  loglik = "Log-likelihood",
  condnum = "Condition Number",
  EL = "Entropy Loss",
  RMSE = "RMSE",
  bRMSE = "Error between estimates",
  iter = "Number of iterations"
)

# Plot y-labs
ylabs <- c(
  time = "Time (s)",
  loglik = TeX("$\\ell$"),
  condnum = TeX('cond(\\hat{\\Sigma})'),
  EL = "Entropy Loss",
  RMSE = "RMSE",
  bRMSE = "RMSE",
  iter = "Iterations"
)

compute.lim <- metrics.names # For which metrics should the plot limits be computed?
plot.cvgl <- setdiff(metrics.names, c("bRMSE")) # For which metrics should the covglasso-related
                                                  # metrics be plotted?

for(nm in metrics.names){
  if(nm %in% compute.lim){
    tempM <- 0
    tempm <- Inf
    
    for(d in 1:N.d){
      temp1 <- colMeans(metrics[nm, d, , , "CVGL"])
      temp2 <- colMeans(metrics[nm, d, , , "GICF"])
      
      tempM <- max(tempM, max(temp1, temp2))
      tempm <- min(tempm, min(temp1, temp2))
    }
  }else{
    tempM <- 1
    tempm <- 0
  }
  
  plot(p.seq, colMeans(metrics[nm, 1, , , "GICF"]),
       ylim = c(tempm, tempM),
       col = clrs[1], type = "l", log = "x",
       ylab = ylabs[nm],
       xlab = "p",
       main = titles[nm],
       lwd = 1)
  points(p.seq, colMeans(metrics[nm, 1, , , "GICF"]),
         col = clrs[1], pch = 20 + 1, bg = clrs[1] )
  
  if(nm %in% plot.cvgl){
    lines(p.seq, colMeans(metrics[nm, 1, , , "CVGL"]),
          col = clrs[1], type = "l", lty = 2)
    points(p.seq, colMeans(metrics[nm, 1, , , "CVGL"]),
           col = clrs[1], pch = 20 + 1 )
  }
  
  for(d in 2:N.d){
    lines(p.seq, colMeans(metrics[nm, d, , , "GICF"]),
          col = clrs[d], type = "l", lwd = 1)
    points(p.seq, colMeans(metrics[nm, d, , , "GICF"]),
           col = clrs[d], pch = 20 + d, bg = clrs[d] )
    
    if(nm %in% plot.cvgl){
      lines(p.seq, colMeans(metrics[nm, d, , , "CVGL"]),
            col = clrs[d], type = "l", lty = 2)
      points(p.seq, colMeans(metrics[nm, d, , , "CVGL"]),
             col = clrs[d], pch = 20 + d )
    }
  }
  
  legend("topright", col = clrs, pch = (20 + 1:N.d), lty = 1,
         legend = c("3", "10"), title = "N. of bands",
         pt.bg = clrs, lwd = 2)
}

par <- or.par

## Table ====

bands <- c(3, 10)
for(d in 1:length(d.seq)){
  print(paste(bands[d], "BANDS"))
  
  print(paste("p___", "GICF___", "CVGL___"))
  
  for(q in 1:length(p.seq)){
    print(paste(
      stringr::str_pad(p.seq[q], 4),
      
      stringr::str_pad(round(
        mean(metrics["time", d, , q, "GICF"]), 4
      ), 7),
      
      stringr::str_pad(round(
        mean(metrics["time", d, , q, "CVGL"]), 4
      ), 7)
    ))
  }
}