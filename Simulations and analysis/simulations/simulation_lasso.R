# HEADER ====
rm(list = ls()) # clear environment
source("../gicf/gicf_source.R")

set.seed(1234)

d.seq <- c(1, 3, 5) # Sequence of 10*(true graph densities)
n.seq <- c(45, 75, 100, 250, 500, 1000) # Sequence of dataset size
p <- 50 # Model size

kappa.max <- 10 # Maximum value of \kappa to be examined

N.sim <- 10 # Number of simulations
N.folds <- 5 # Number of CV folds
N.kappa <- 18 # Number of values of \kappa to sample
N.lambda <- 18 # Number of values of \lambda to sample
N.d <- length(d.seq)
N.n <- length(n.seq)

zero <- 1e-4 # Tolerance: absolute values below this threshold are treated as 0

# Parameters to access simulated data (DO NOT CHANGE)
model <- "RB" # Model name: Random entries with Banded structure
n.key <- 2000 # Key to access simulated data

# METRICS ====
metrics.names <- c(
  "lambda", # optimal \lambda
  "lambdar0", # (optimal \lambda)/(\lambda_max(0))
  "kappa", # optimal \kappa
  "condnum", # condition number
  "EL", # Entropy Loss
  "RMSE", # Root Mean Square Error
  "F1", # F_1 score
  "eTPR", # True positive rate
  "eTNR", # True negative rate
  "ePPV", # Positive predicted value
  "d" #estimated density
)

metrics <- array(zero, c(length(metrics.names), N.d, N.sim, N.n, 2)) # [metric name, density, simulation, n, LASSO/LRIDGE]
dimnames(metrics)[[1]] <- metrics.names
dimnames(metrics)[[5]] <- c("LASSO", "LRIDGE")

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

eTNR <- function(A, Sigma, toll = zero){
  true.negative <- (abs(A) <= zero) & (abs(Sigma) <= zero) 
  negative <- (abs(Sigma) <= zero)
  
  if(sum(negative) == 0)
    return(1)
  
  return(sum(true.negative)/sum(negative))
}

eTPR <- function(A, Sigma, toll = zero){
  p <- dim(Sigma)[1]
  
  true.positive <- (abs(A) > toll) & (abs(Sigma) > toll)
  positive <- (abs(Sigma) > toll)
  
  if(sum(positive) - p == 0)
    return(1)
  
  return( (sum(true.positive) - p)/(sum(positive) - p) )
}

ePPV <- function(A, Sigma, toll = zero){
  p <- dim(Sigma)[1]
  
  true.positive <- (abs(A) > toll) & (abs(Sigma) > toll)
  predicted.positive <- (abs(A) > toll)
  
  if(sum(predicted.positive) - p == 0)
    return(1)
  
  return( (sum(true.positive) - p)/(sum(predicted.positive) - p) )
}

f1 <- function(A, Sigma, toll = zero){
  prec <- ePPV(A, Sigma, toll)
  rec <- eTPR(A, Sigma, toll)
  
  f <- 2*(prec*rec)/(prec + rec)
  
  return( ifelse(is.na(f), 0 , f) )
}

# MODEL SELECTION ====
model.selection.cv <- function(y, # data
                               N.folds, # num. of CV folds
                               N.l, # num. of lambda candidates
                               N.k, # num. of kappa candidates
                               l.min = 0, # min value of \lambda to be explored
                               l.max = Inf, # max value of \lambda to be explored
                               k.max = 10 # max value of \kappa to be explored
                              ){
  n <- nrow(y)
  p <- ncol(y)
  
  S <- cov(y) * (n - 1)/n
  
  folded.data <- fold(data.frame(y), k = N.folds)
  
  # (kappa, lambda, val)
  model.opt.lridge <- rep(-Inf, 3) # lasso + ridge regularisation
  model.opt.lasso <- rep(-Inf, 3) # lasso regularisation
  
  lambda.max <- min( max.lambda(S), l.max )
  lambda.seq <- seq(l.min, lambda.max, length.out = N.l)
  
  kappa.max <- min( max.kappa(S, l.min), k.max )
  
  for(L in 1:N.l){
    lambda.loc <- lambda.seq[L]
    
    kappa.max.loc <- min( max.kappa(S, lambda.loc), k.max )
    seq.length.loc <- max( ceiling(N.k * kappa.max.loc/kappa.max) , 3 ) # Keep "density" of points approx. constant
    
    if(n > p)
      kappa.seq <- seq(0, kappa.max.loc, length.out = seq.length.loc)
    else
      kappa.seq <- seq(0, kappa.max.loc, length.out = seq.length.loc + 1)[-1]
    
    for(K in 1:seq.length.loc){
      kappa.loc <- kappa.seq[K]
      val <- 0
      
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
        
        D <- diag(kappa.loc, p)
        fit <- gicf.run(S = S.train + D, n = n.train, lambda = lambda.loc)$sigma
        
        val <- val + gicf.likelihood(fit, S.test, n.test, 0)
      }
      
      if(n > p && K == 1){
        if(val > model.opt.lasso[3])
          model.opt.lasso <- c(0, lambda.loc, val)
      }
      
      if(val > model.opt.lridge[3])
        model.opt.lridge <- c(kappa.loc, lambda.loc, val)
    }
  }
  
  return(list(lasso = model.opt.lasso, lridge = model.opt.lridge))
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
      
      selected.models <- model.selection.cv(y, 
                                            N.folds, 
                                            N.lambda, 
                                            N.kappa,
                                            0.01*max.lambda(S),
                                            0.95*max.lambda(S),
                                            kappa.max)
      
      # Compute LASSO (if possible)
      if(n.loc > p){
        model.lasso <- selected.models$lasso
        sigma.lasso <- gicf.run(S = S, n = n.loc, lambda = model.lasso[2])$sigma # MLE
        
        metrics["lambda", d, s, m, "LASSO"] <- model.lasso[2]
        metrics["lambdar0", d, s, m, "LASSO"] <- model.lasso[2]/max.lambda(S)
        metrics["kappa", d, s, m, "LASSO"] <- 0
        metrics["condnum", d, s, m, "LASSO"] <- condnum(sigma.lasso)
        metrics["EL", d, s, m, "LASSO"] <- entropy.loss(sigma.lasso, Sigma, detsigma, Theta)
        metrics["RMSE", d, s, m, "LASSO"] <- rmse(sigma.lasso, Sigma)
        metrics["eTPR", d, s, m, "LASSO"] <- eTPR(sigma.lasso, Sigma)
        metrics["eTNR", d, s, m, "LASSO"] <- eTNR(sigma.lasso, Sigma)
        metrics["ePPV", d, s, m, "LASSO"] <- ePPV(sigma.lasso, Sigma)
        metrics["F1", d, s, m, "LASSO"] <- f1(sigma.lasso, Sigma)
        metrics["d", d, s, m, "LASSO"] <- (sum(abs(sigma.lasso) > zero) - p)/(p*(p-1))
      }
      
      # Compute shrinked estimate
      model.lridge <- selected.models$lridge
      D <- diag(model.lridge[1], p)
      sigma.lridge <- gicf.run(S = S + D, n = n.loc, lambda = model.lridge[2])$sigma
      
      metrics["lambda", d, s, m, "LRIDGE"] <- model.lridge[2]
      metrics["lambdar0", d, s, m, "LRIDGE"] <- model.lridge[2]/max.lambda(S)
      metrics["kappa", d, s, m, "LRIDGE"] <- model.lridge[1]
      metrics["condnum", d, s, m, "LRIDGE"] <- condnum(sigma.lridge)
      metrics["EL", d, s, m, "LRIDGE"] <- entropy.loss(sigma.lridge, Sigma, detsigma, Theta)
      metrics["RMSE", d, s, m, "LRIDGE"] <- rmse(sigma.lridge, Sigma)
      metrics["eTPR", d, s, m, "LRIDGE"] <- eTPR(sigma.lridge, Sigma)
      metrics["eTNR", d, s, m, "LRIDGE"] <- eTNR(sigma.lridge, Sigma)
      metrics["ePPV", d, s, m, "LRIDGE"] <- ePPV(sigma.lridge, Sigma)
      metrics["F1", d, s, m, "LRIDGE"] <- f1(sigma.lridge, Sigma)
      metrics["d", d, s, m, "LRIDGE"] <- (sum(abs(sigma.lridge) > zero) - p)/(p*(p - 1))
    }
    
    print("")
  }
}

# OUTPUT ====
library(latex2exp)

or.par <- par(mar = c(5, 4, 4, 2) + 0.5)
clrs <- palette.colors(palette = "Okabe-Ito", n = N.d)

# plot titles
titles <- c(
  lambda = TeX("$\\hat{\\lambda}$"),
  lambdar0 = "Lasso Parameter",
  kappa = TeX("\\kappa"),
  condnum = "Condition Number Ratio",
  EL = "Entropy Loss",
  RMSE = "RMSE",
  F1 = TeX("$F_{1}$"),
  eTPR = "eTPR",
  eTNR = "eTNR",
  ePPV = "ePPV",
  d = "Estimated Density"
)

# plot y-labels
ylabs <- c(
  lambda = TeX("$\\hat{\\lambda}$"),
  lambdar0 = TeX("$\\hat{\\lambda}/\\lambda_{MAX}(0)$"),
  kappa = TeX("\\kappa"),
  condnum = TeX('cond(\\hat{\\Sigma})/cond(\\Sigma)'),
  EL = "Entropy Loss",
  RMSE = "RMSE",
  F1 = TeX("$F_{1}$ Score"),
  eTPR = "eTPR",
  eTNR = "eTNR",
  ePPV = "ePPV",
  d = "Est. density"
)

# Placement of legends
placements <- c(
  lambda = "bottomright",
  lambdar0 = "bottomright",
  kappa = "topright",
  condnum = "topleft",
  EL = "bottomright",
  RMSE = "bottomright",
  F1 = "bottomright",
  eTPR = "bottomright",
  eTNR = "bottomleft",
  ePPV = "bottomright",
  d = "bottomright"
)

compute.lim <- setdiff(metrics.names, c("eTPR", "eTNR", "ePPV", "F1", "d")) # For which metrics should the plot limits be computed?
highlight.level <- c( # Which plots need to highlight a level? 
  condnum = max(true.condnum[1]) 
)
plot.mle <- setdiff(metrics.names, c("kappa"))# For which metrics should the mle-related
                                               # metrics be plotted?

for(nm in metrics.names){
  if(nm %in% compute.lim){
    tempM <- 0
    tempm <- Inf
    
    for(d in 1:N.d){
      temp1 <- colMeans(metrics[nm, d, , , "LASSO"])[-1]
      temp2 <- colMeans(metrics[nm, d, , , "LRIDGE"])
      
      tempM <- max(tempM, max(temp1, temp2))
      tempm <- min(tempm, min(temp1, temp2))
    }
  }else{
    tempM <- 1
    tempm <- 0
  }
  
  if(nm %in% names(highlight.level))
      tempM <- max(tempM, highlight.level[nm])

  plot(n.seq, colMeans(metrics[nm, 1, , , "LRIDGE"]),
       ylim = c(tempm, tempM),
       col = clrs[1], type = "l", log = "x",
       ylab = ylabs[nm],
       xlab = "n",
       main = titles[nm],
       lwd = 1)
  points(n.seq, colMeans(metrics[nm, 1, , , "LRIDGE"]),
         col = clrs[1], pch = 20 + 1, bg = clrs[1] )
  
  if(nm %in% plot.mle){
    lines(n.seq[-1], colMeans(metrics[nm, 1, , , "LASSO"])[-1],
          col = clrs[1], type = "l", lty = 2)
    points(n.seq[-1], colMeans(metrics[nm, 1, , , "LASSO"])[-1],
           col = clrs[1], pch = 20 + 1 )
  }
  
  for(d in 2:N.d){
    lines(n.seq, colMeans(metrics[nm, d, , , "LRIDGE"]),
          col = clrs[d], type = "l", lwd = 1)
    points(n.seq, colMeans(metrics[nm, d, , , "LRIDGE"]),
           col = clrs[d], pch = 20 + d, bg = clrs[d] )
    
    if(nm %in% plot.mle){
      lines(n.seq[-1], colMeans(metrics[nm, d, , , "LASSO"])[-1],
            col = clrs[d], type = "l", lty = 2)
      points(n.seq[-1], colMeans(metrics[nm, d, , , "LASSO"])[-1],
             col = clrs[d], pch = 20 + d )
    }
  }
  
  if(nm %in% names(highlight.level))
    abline(h = highlight.level[nm], lwd = 2, col = "gray")
  
  legend(placements[nm], col = clrs, pch = (20 + d.seq), lty = 1,
         legend = paste0(d.seq*10, "%"), title = "Density",
         pt.bg = clrs, lwd = 2)
}

par <- or.par