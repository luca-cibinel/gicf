rm(list = ls()) # clear environment

path.join <- function(...){
  sep <- .Platform$file.sep
  pat <- paste0("\\", sep, "$")
  
  paths <- lapply(list(...), function(x){gsub(pat, "", x)})
  
  return(paste0(paths, collapse = sep))
}

# HEADER ====
source("../gicf/gicf.R")

set.seed(1234)

b.seq <- c(3, 10) # Sequence of n. of bands
n <- 2000 # Dataset size
p.seq <- c(25, 50, 100, 200, 300) # Sequence of model sizes

kappa.max <- 5 # Maximum value of \kappa to be examined

simulation.batch <- 1:10 # Desired simulations
N.sim <- length(simulation.batch)
N.folds <- 5 # Number of CV folds
N.kappa <- 30 # Number of values of \kappa to sample
N.lambda <- 20 # Number of values of \lambda to sample
N.b <- length(b.seq)
N.p <- length(p.seq)

zero <- 1e-4 # Tolerance: absolute values below this threshold are treated as 0

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
  "d", #estimated density
  "iters", # N. of iterations
  "time.user.self", # Computational time (fit only) "user.self" entry
  "time.sys.self", # Computational time (fit only) "sys.self" entry
  "time.elapsed", # Computational time (fit only) "elapsed" entry
  "time.user.child", # Computational time (fit only) "user.child" entry
  "time.sys.child", # Computational time (fit only) "sys.child" entry
  "cv.iters", # N. of iterations (CV, average)
  "cv.time.user.self", # Computational time (CV only, average) "user.self" entry
  "cv.time.sys.self", # Computational time (CV only, average) "sys.self" entry
  "cv.time.elapsed", # Computational time (CV only, average) "elapsed" entry
  "cv.time.user.child", # Computational time (CV only, total) "user.child" entry
  "cv.time.sys.child", # Computational time (CV only, total) "sys.child" entry
  "cv.all.time.user.self", # Computational time (CV only, total) "user.self" entry
  "cv.all.time.sys.self", # Computational time (CV only, total) "sys.self" entry
  "cv.all.time.elapsed", # Computational time (CV only, total) "elapsed" entry
  "cv.all.time.user.child", # Computational time (CV only, total) "user.child" entry
  "cv.all.time.sys.child", # Computational time (CV only, total) "sys.child" entry
  "n.cv.par" # N. of CV candidates
)

metrics <- array(zero, c(length(metrics.names), N.b, N.sim, N.p, 2)) # [metric name, density, simulation, n, LASSO/LRIDGE]
dimnames(metrics)[[1]] <- metrics.names
dimnames(metrics)[[2]] <- sapply(b.seq, toString)
dimnames(metrics)[[3]] <- sapply(simulation.batch, toString)
dimnames(metrics)[[4]] <- sapply(p.seq, toString)
dimnames(metrics)[[5]] <- c("LASSO", "LRIDGE")

true.condnum <- rep(0, N.b) # Store true condition number

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
model.selection.cv.gicf <- function(y, # data
                               adj, # adjacency matrix
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
  
  lambda.max <- min( lambdamax(S), l.max )
  lambda.seq <- seq(l.min, lambda.max, length.out = N.l)
  
  kappa.max <- min( kappamax(S, l.min), k.max )
  
  n.of.pars <- 0
  times <- NULL
  iters <- NULL
  
  for(L in 1:N.l){
    lambda.loc <- lambda.seq[L]
    
    kappa.max.loc <- min( kappamax(S, lambda.loc), k.max )
    seq.length.loc <- max( ceiling(N.k * kappa.max.loc/kappa.max) , 3 ) # Keep "density" of points approx. constant
    
    if(n > p)
      kappa.seq <- seq(0, log(kappa.max.loc + 1), length.out = seq.length.loc)
    else
      kappa.seq <- seq(0, log(kappa.max.loc + 1), length.out = seq.length.loc + 1)[-1]
    
    kappa.seq <- exp(kappa.seq) - 1
    
    for(K in 1:seq.length.loc){
      n.of.pars <- n.of.pars + 1
      
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
        
        loc.time <- system.time(
          fit <- gicf(S = S.train, n = n.train, adj = adj, lambda = lambda.loc, kappa = kappa.loc)
        )
        val <- val + gcgmloglik(fit$sigma, S.test, n.test)
        
        times <- rbind(times, c(kappa.loc, lambda.loc, loc.time))
        iters <- rbind(iters, c(kappa.loc, lambda.loc, fit$it))
      }
      
      if(val > model.opt.lridge[3])
        model.opt.lridge <- c(kappa.loc, lambda.loc, val)
    }
  }
  
  return(list(lridge = model.opt.lridge, n.of.pars = n.of.pars, cv.time = times, cv.iters = iters))
}

model.selection.cv.lasso <- function(y, # data
                                    adj, # adjacecny matrix
                                    N.folds, # num. of CV folds
                                    N.l, # num. of lambda candidates
                                    # N.k, # num. of kappa candidates
                                    l.min = 0, # min value of \lambda to be explored
                                    l.max = Inf, # max value of \lambda to be explored
                                    k.max = 10 # max value of \kappa to be explored
){
  n <- nrow(y)
  p <- ncol(y)
  
  L <- matrix(1, p, p)
  L[adj == 0] <- Inf
  
  S <- cov(y) * (n - 1)/n
  
  folded.data <- fold(data.frame(y), k = N.folds)
  
  # (kappa, lambda, val)
  model.opt.lasso <- rep(-Inf, 3) # lasso regularisation
  
  lambda.max <- min( lambdamax(S), l.max )
  lambda.seq <- seq(l.min, lambda.max, length.out = N.l)
  
  n.of.pars <- 0
  times <- NULL
  iters <- NULL
  
  for(L in 1:N.l){
    lambda.loc <- lambda.seq[L]*L
    
    n.of.pars <- n.of.pars + 1
    
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
      
      loc.time <- system.time(
        fit <- gicf(S = S.train, n = n.train, lambda = lambda.loc)
      )
      val <- val + gcgmloglik(fit$sigma, S.test, n.test)
      
      times <- rbind(times, c(0, lambda.loc, loc.time))
      iters <- rbind(iters, c(0, lambda.loc, fit$it))
    }
    
    if(val > model.opt.lasso[3])
      model.opt.lasso <- c(0, lambda.loc, val)
  }
  
  return(list(lasso = model.opt.lasso, n.of.pars = n.of.pars, cv.time = times, cv.iters = iters))
}

# SIMULATION ====
for(b in 1:N.b){ # For each n. of bands
  n.bands <- b.seq[b]
  
  print("")
  print("")
  print(paste("N. OF BANDS:", n.bands))
  
  
  for(q in 1:N.p){ # Repeat for each simulation:
    p <- p.seq[q]
    
    print(paste("p =", p))
    
    Sigma <- as.matrix( # Read true Sigma from simulated data
      read.table(
        path.join("data", paste0("p", p), paste0("b", n.bands), "sigma.dat")
      )
    )
    
    Theta <- solve(Sigma) # \Sigma^{-1}
    detsigma <- det(Sigma) # det(\Sigma)
    
    true.condnum[b] <- condnum(Sigma) # Store true condition number
    
    adj <- matrix(0, p, p)
    adj[Sigma != 0] <- 1
    
    adj.lasso <- adj
    adj.lasso[adj == 0] <- Inf
    
    for(s in simulation.batch){ # For each desired dataset size:
      print(paste("Simulation:", s))
      
      data <- as.matrix( # Read simulated data
        read.table(
          path.join("data", paste0("p", p), paste0("b", n.bands), paste0("sim", simulation.batch[s], ".dat"))
        )
      )
      
      y <- data[1:n, ] # Extract local dataset
      
      S <- cov(y) * (n - 1)/n # Compute local sample covariance matrix
      
      # Select models
      selected.models.lridge <- model.selection.cv.gicf(y,
                                            adj,
                                            N.folds, 
                                            N.lambda, 
                                            N.kappa,
                                            0*0.01*lambdamax(S),
                                            0.95*lambdamax(S),
                                            kappa.max)
      
      cv.time.lridge <- colMeans(selected.models.lridge$cv.time[,-c(1,2)])
      cv.iters.lridge <- mean(selected.models.lridge$cv.iters[,-c(1,2)])
      cv.all.time.lridge <- colMeans(selected.models.lridge$cv.time[,-c(1,2)])
      
      selected.models.lasso <- model.selection.cv.lasso(y,
                                                        adj,
                                                        N.folds, 
                                                        N.lambda,
                                                        0*0.01*lambdamax(S),
                                                        0.95*lambdamax(S),
                                                        kappa.max)
      cv.time.lasso <- colMeans(selected.models.lasso$cv.time[, -c(1,2)])
      cv.iters.lasso <- mean(selected.models.lasso$cv.time[, -c(1,2)])
      cv.all.time.lasso <- colSums(selected.models.lasso$cv.time[, -c(1,2)])
      
      # Compute LASSO (if possible)
      model.lasso <- selected.models.lasso$lasso
      
      time <- system.time(
        sigma.lasso <- gicf(S = S, n = n, lambda = model.lasso[2]*adj.lasso)
      )
      iters.lasso <- sigma.lasso$it
      sigma.lasso <- sigma.lasso$sigma
      
      metrics["lambda", b, s, q, "LASSO"] <- model.lasso[2]
      metrics["lambdar0", b, s, q, "LASSO"] <- model.lasso[2]/lambdamax(S)
      metrics["kappa", b, s, q, "LASSO"] <- 0
      metrics["condnum", b, s, q, "LASSO"] <- condnum(sigma.lasso)
      metrics["EL", b, s, q, "LASSO"] <- entropy.loss(sigma.lasso, Sigma, detsigma, Theta)
      metrics["RMSE", b, s, q, "LASSO"] <- rmse(sigma.lasso, Sigma)
      metrics["eTPR", b, s, q, "LASSO"] <- eTPR(sigma.lasso, Sigma)
      metrics["eTNR", b, s, q, "LASSO"] <- eTNR(sigma.lasso, Sigma)
      metrics["ePPV", b, s, q, "LASSO"] <- ePPV(sigma.lasso, Sigma)
      metrics["F1", b, s, q, "LASSO"] <- f1(sigma.lasso, Sigma)
      metrics["d", b, s, q, "LASSO"] <- (sum(abs(sigma.lasso) > zero) - p)/(p*(p-1))
      metrics["n.cv.par", b, s, q, "LASSO"] <- N.lambda
      
      metrics["iters", b, s, q, "LASSO"] <- iters.lasso
      metrics["cv.iters", b, s, q, "LASSO"] <- cv.iters.lasso
      
      for(tname in names(time)){
        metrics[paste0("time.", tname), b, s, q, "LASSO"] <- time[tname]
        metrics[paste0("cv.time.", tname), b, s, q, "LASSO"] <- cv.time.lasso[tname]
        metrics[paste0("cv.all.time.", tname), b, s, q, "LASSO"] <- cv.all.time.lasso[tname]
      }
      
      
      # Compute shrinked estimate
      model.lridge <- selected.models.lridge$lridge
      
      time <- system.time(
        sigma.lridge <- gicf(S = S, n = n, lambda = model.lridge[2], kappa = model.lridge[1], adj = adj)
      )
      iters.lridge <- sigma.lridge$it
      sigma.lridge <- sigma.lridge$sigma
      
      metrics["lambda", b, s, q, "LRIDGE"] <- model.lridge[2]
      metrics["lambdar0", b, s, q, "LRIDGE"] <- model.lridge[2]/lambdamax(S)
      metrics["kappa", b, s, q, "LRIDGE"] <- model.lridge[1]
      metrics["condnum", b, s, q, "LRIDGE"] <- condnum(sigma.lridge)
      metrics["EL", b, s, q, "LRIDGE"] <- entropy.loss(sigma.lridge, Sigma, detsigma, Theta)
      metrics["RMSE", b, s, q, "LRIDGE"] <- rmse(sigma.lridge, Sigma)
      metrics["eTPR", b, s, q, "LRIDGE"] <- eTPR(sigma.lridge, Sigma)
      metrics["eTNR", b, s, q, "LRIDGE"] <- eTNR(sigma.lridge, Sigma)
      metrics["ePPV", b, s, q, "LRIDGE"] <- ePPV(sigma.lridge, Sigma)
      metrics["F1", b, s, q, "LRIDGE"] <- f1(sigma.lridge, Sigma)
      metrics["d", b, s, q, "LRIDGE"] <- (sum(abs(sigma.lridge) > zero) - p)/(p*(p - 1))
      metrics["n.cv.par", b, s, q, "LRIDGE"] <- selected.models.lridge$n.of.pars
      
      metrics["iters", b, s, q, "LRIDGE"] <- iters.lridge
      metrics["cv.iters", b, s, q, "LRIDGE"] <- cv.iters.lridge
      
      for(tname in names(time)){
        metrics[paste0("time.", tname), b, s, q, "LRIDGE"] <- time[tname]
        metrics[paste0("cv.time.", tname), b, s, q, "LRIDGE"] <- cv.time.lridge[tname]
        metrics[paste0("cv.all.time.", tname), b, s, q, "LRIDGE"] <- cv.all.time.lridge[tname]
      }
    }
    
    print("")
  }
}

df <- array2DF(metrics)
colnames(df) <- c("metric", "n_bands", "simulation", "p", "method", "value")
write.table(df, path.join("results", "simulation_time.csv"), row.names = F)

# OUTPUT ====

# Read results (designed for independent plotting)

results <- read.table(path.join("results", "simulation_lasso.csv"), header = 1)

dims <- sapply(1:(ncol(results) - 1), function(i){length(unique(results[,i]))})
dimnms <- list()

for(i in 1:(ncol(results) - 1)){
  dimnms[[i]] <- unique(results[,i])
}

results <- array(results$value, dims)
dimnames(results) <- dimnms
metric.names <- dimnms[[1]]

# Start plotting
library(latex2exp)

or.par <- par(mar = c(5, 4, 4, 2) + 0.5)
clrs <- palette.colors(palette = "Okabe-Ito", n = N.b)
names(clrs) <- dimnms[[2]]
band.color.palette <- scale_color_manual(name = "Var2", values = clrs)

ltys <- 1:2
names(ltys) <- c("LRIDGE", "LASSO")
method.linetype <- scale_linetype_manual(name = "Var3", values = ltys)

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
  d = "Estimated Density",
  time.user.self = "Time (user.self)",
  time.sys.self = "Time (sys.self)",
  time.elapsed = "Time (elapsed)",
  time.user.child = "Time (sys.child)",
  time.sys.child = "Time (user.self)",
  cv.time.user.self = "CV-Time (user.self)",
  cv.time.sys.self = "CV-Time (sys.self)",
  cv.time.elapsed = "CV-Time (elapsed)",
  cv.time.user.child = "CV-Time (sys.child)",
  cv.time.sys.child = "CV-Time (user.self)"
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
  d = "Est. density",
  time.user.self = "Time (s)",
  time.sys.self = "Time (s)",
  time.elapsed = "Time (s)",
  time.user.child = "Time (s)",
  time.sys.child = "Time (s)",
  cv.time.user.self = "CV-Time (s)",
  cv.time.sys.self = "CV-Time (s)",
  cv.time.elapsed = "CV-Time (s)",
  cv.time.user.child = "CV-Time (s)",
  cv.time.sys.child = "CV-Time (s)"
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
    
    for(d in 1:N.b){
      temp1 <- colMeans(results[nm, d, , , "LASSO"])[-1]
      temp2 <- colMeans(results[nm, d, , , "LRIDGE"])
      
      tempM <- max(tempM, max(temp1, temp2))
      tempm <- min(tempm, min(temp1, temp2))
    }
  }else{
    tempM <- 1
    tempm <- 0
  }
  
  if(nm %in% names(highlight.level))
    tempM <- max(tempM, highlight.level[nm])
  
  plot.data.lasso <-array2DF(apply(results[nm,,,-1,"LASSO"], 1, colMeans))
  plot.data.lasso$Var1 <- sapply(plot.data.lasso$Var1, as.integer)
  plot.data.lridge <-array2DF(apply(results[nm,,,,"LRIDGE"], 1, colMeans))
  plot.data.lridge$Var1 <- sapply(plot.data.lridge$Var1, as.integer)
  
  g <- ggplot() +
    scale_x_log10() + 
    band.color.palette + 
    method.linetype +
    geom_line(data = plot.data.lasso, aes(x = Var1, y = Value, color = Var2, linetype = "LASSO")) + 
    geom_point(data = plot.data.lasso, aes(x = Var1, y = Value, pch = Var2, color = Var2), size = 2.1) +
    geom_line(data = plot.data.lridge, aes(x = Var1, y = Value, color = Var2, linetype = "LRIDGE")) + 
    geom_point(data = plot.data.lridge, aes(x = Var1, y = Value, pch = Var2, color = Var2), size = 2.1)
  print(g)
  
  plot(n.seq, colMeans(results[nm, 1, , , "LRIDGE"]),
       ylim = c(tempm, tempM),
       col = clrs[1], type = "l", log = "x",
       ylab = ylabs[nm],
       xlab = "n",
       main = titles[nm],
       lwd = 1)
  points(n.seq, colMeans(results[nm, 1, , , "LRIDGE"]),
         col = clrs[1], pch = 20 + 1, bg = clrs[1] )
  
  if(nm %in% plot.mle){
    lines(n.seq[-1], colMeans(results[nm, 1, , , "LASSO"])[-1],
          col = clrs[1], type = "l", lty = 2)
    points(n.seq[-1], colMeans(results[nm, 1, , , "LASSO"])[-1],
           col = clrs[1], pch = 20 + 1 )
  }
  
  for(d in 2:N.b){
    lines(n.seq, colMeans(results[nm, d, , , "LRIDGE"]),
          col = clrs[d], type = "l", lwd = 1)
    points(n.seq, colMeans(results[nm, d, , , "LRIDGE"]),
           col = clrs[d], pch = 20 + d, bg = clrs[d] )
    
    if(nm %in% plot.mle){
      lines(n.seq[-1], colMeans(results[nm, d, , , "LASSO"])[-1],
            col = clrs[d], type = "l", lty = 2)
      points(n.seq[-1], colMeans(results[nm, d, , , "LASSO"])[-1],
             col = clrs[d], pch = 20 + d )
    }
  }
  
  if(nm %in% names(highlight.level))
    abline(h = highlight.level[nm], lwd = 2, col = "gray")
  
  # legend(placements[nm], col = clrs, pch = (20 + b.seq), lty = 1,
  #        legend = paste0(b.seq*10, "%"), title = "Density",
  #        pt.bg = clrs, lwd = 2)
}

par <- or.par