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

b.seq <- c(14) # Sequence of n. of bands
n.seq <- c(45, 75, 100, 250, 500, 1000) # Sequence of dataset size
p <- 50 # Model size

kappa.max <- 5 # Maximum value of \kappa to be examined

simulation.batch <- 1:20 # Desired simulations
N.sim <- length(simulation.batch)
N.folds <- 5 # Number of CV folds
N.kappa <- 30 # Number of values of \kappa to sample
N.lambda <- 20 # Number of values of \lambda to sample
N.b <- length(b.seq)
N.n <- length(n.seq)

zero <- 1e-4 # Tolerance: absolute values below this threshold are treated as 0

backup.period <- 2 # After how many simulations should a partial result be stored for backup?

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
  "cv.iters", # N. of iterations (CV)
  "cv.time.user.self", # Computational time (CV only) "user.self" entry
  "cv.time.sys.self", # Computational time (CV only) "sys.self" entry
  "cv.time.elapsed", # Computational time (CV only) "elapsed" entry
  "cv.time.user.child", # Computational time (CV only) "user.child" entry
  "cv.time.sys.child", # Computational time (CV only) "sys.child" entry
  "n.cv.par" # N. of CV candidates
)

metrics <- array(zero, c(length(metrics.names), N.b, N.sim, N.n, 2)) # [metric name, density, simulation, n, LASSO/LRIDGE]
dimnames(metrics)[[1]] <- metrics.names
dimnames(metrics)[[2]] <- sapply(b.seq, toString)
dimnames(metrics)[[3]] <- sapply(simulation.batch, toString)
dimnames(metrics)[[4]] <- sapply(n.seq, toString)
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
          fit <- gicf(S = S.train, n = n.train, lambda = lambda.loc, kappa = kappa.loc)
        )
        val <- val + gcgmloglik(fit$sigma, S.test, n.test)
        
        times <- rbind(times, c(kappa.loc, lambda.loc, loc.time))
        iters <- rbind(iters, c(kappa.loc, lambda.loc, fit$it))
      }
      
      if(n > p && K == 1){
        if(val > model.opt.lasso[3])
          model.opt.lasso <- c(0, lambda.loc, val)
      }
      
      if(val > model.opt.lridge[3])
        model.opt.lridge <- c(kappa.loc, lambda.loc, val)
    }
  }
  
  return(list(lasso = model.opt.lasso, lridge = model.opt.lridge, n.of.pars = n.of.pars, cv.time = times, cv.iters = iters))
}

# SIMULATION ====
for(b in 1:N.b){ # For each n. of bands
  n.bands <- b.seq[b]
  
  print("")
  print("")
  print(paste("N. OF BANDS:", n.bands))
  
  Sigma <- as.matrix( # Read true Sigma from simulated data
    read.table(
      path.join("data", paste0("p", p), paste0("b", n.bands), "sigma.dat")
    )
  )
  
  Theta <- solve(Sigma) # \Sigma^{-1}
  detsigma <- det(Sigma) # det(\Sigma)
  
  true.condnum[b] <- condnum(Sigma) # Store true condition number
  
  for(s in 1:N.sim){ # Repeat for each simulation:
    print(paste("Simulation:", s))
    
    data <- as.matrix( # Read simulated data
      read.table(
        path.join("data", paste0("p", p), paste0("b", n.bands), paste0("sim", simulation.batch[s], ".dat"))
      )
    )
    
    for(m in 1:N.n){ # For each desired dataset size:
      cat(paste("n =", n.seq[m], "-"))
      n.loc <- n.seq[m] # Local dataset size
      
      y <- data[1:n.loc, ] # Extract local dataset
      
      S <- cov(y) * (n.loc - 1)/n.loc # Compute local sample covariance matrix
      
      selected.models <- model.selection.cv(y, 
                                            N.folds, 
                                            N.lambda, 
                                            N.kappa,
                                            0*0.01*lambdamax(S),
                                            0.95*lambdamax(S),
                                            kappa.max)
      
      cv.time <- colMeans(selected.models$cv.time[,-c(1,2)])
      cv.iters <- mean(selected.models$cv.iters[,-c(1,2)])
      cv.time.lasso <- colMeans(selected.models$cv.time[selected.models$cv.time[,1] == 0, -c(1,2)])
      cv.iters.lasso <- mean(selected.models$cv.time[selected.models$cv.time[,1] == 0, -c(1,2)])
      
      # Compute LASSO (if possible)
      if(n.loc > p){
        model.lasso <- selected.models$lasso
        
        time <- system.time(
          sigma.lasso <- gicf(S = S, n = n.loc, lambda = model.lasso[2])
        )
        iters.lasso <- sigma.lasso$it
        sigma.lasso <- sigma.lasso$sigma
        
        metrics["lambda", b, s, m, "LASSO"] <- model.lasso[2]
        metrics["lambdar0", b, s, m, "LASSO"] <- model.lasso[2]/lambdamax(S)
        metrics["kappa", b, s, m, "LASSO"] <- 0
        metrics["condnum", b, s, m, "LASSO"] <- condnum(sigma.lasso)
        metrics["EL", b, s, m, "LASSO"] <- entropy.loss(sigma.lasso, Sigma, detsigma, Theta)
        metrics["RMSE", b, s, m, "LASSO"] <- rmse(sigma.lasso, Sigma)
        metrics["eTPR", b, s, m, "LASSO"] <- eTPR(sigma.lasso, Sigma)
        metrics["eTNR", b, s, m, "LASSO"] <- eTNR(sigma.lasso, Sigma)
        metrics["ePPV", b, s, m, "LASSO"] <- ePPV(sigma.lasso, Sigma)
        metrics["F1", b, s, m, "LASSO"] <- f1(sigma.lasso, Sigma)
        metrics["d", b, s, m, "LASSO"] <- (sum(abs(sigma.lasso) > zero) - p)/(p*(p-1))
        metrics["n.cv.par", b, s, m, "LASSO"] <- N.lambda
        
        metrics["iters", b, s, m, "LASSO"] <- iters.lasso
        metrics["cv.iters", b, s, m, "LASSO"] <- cv.iters.lasso
        
        for(tname in names(time)){
          metrics[paste0("time.", tname), b, s, m, "LASSO"] <- time[tname]
          metrics[paste0("cv.time.", tname), b, s, m, "LASSO"] <- cv.time.lasso[tname]
        }
      }
      
      # Compute shrinked estimate
      model.lridge <- selected.models$lridge
      
      time <- system.time(
        sigma.lridge <- gicf(S = S, n = n.loc, lambda = model.lridge[2], kappa = model.lridge[1])
      )
      iters.lridge <- sigma.lridge$it
      sigma.lridge <- sigma.lridge$sigma
      
      metrics["lambda", b, s, m, "LRIDGE"] <- model.lridge[2]
      metrics["lambdar0", b, s, m, "LRIDGE"] <- model.lridge[2]/lambdamax(S)
      metrics["kappa", b, s, m, "LRIDGE"] <- model.lridge[1]
      metrics["condnum", b, s, m, "LRIDGE"] <- condnum(sigma.lridge)
      metrics["EL", b, s, m, "LRIDGE"] <- entropy.loss(sigma.lridge, Sigma, detsigma, Theta)
      metrics["RMSE", b, s, m, "LRIDGE"] <- rmse(sigma.lridge, Sigma)
      metrics["eTPR", b, s, m, "LRIDGE"] <- eTPR(sigma.lridge, Sigma)
      metrics["eTNR", b, s, m, "LRIDGE"] <- eTNR(sigma.lridge, Sigma)
      metrics["ePPV", b, s, m, "LRIDGE"] <- ePPV(sigma.lridge, Sigma)
      metrics["F1", b, s, m, "LRIDGE"] <- f1(sigma.lridge, Sigma)
      metrics["d", b, s, m, "LRIDGE"] <- (sum(abs(sigma.lridge) > zero) - p)/(p*(p - 1))
      metrics["n.cv.par", b, s, m, "LRIDGE"] <- selected.models$n.of.pars
      
      metrics["iters", b, s, m, "LRIDGE"] <- iters.lridge
      metrics["cv.iters", b, s, m, "LRIDGE"] <- cv.iters
      
      for(tname in names(time)){
        metrics[paste0("time.", tname), b, s, m, "LRIDGE"] <- time[tname]
        metrics[paste0("cv.time.", tname), b, s, m, "LRIDGE"] <- cv.time[tname]
      }
    }
    
    if(s %% backup.period == 0){
      sims.so.far <- 1:s
      bands.so.far <- b.seq[1:b]
      
      fname <- paste0(
        "simulation_lasso__b", 
        paste0(bands.so.far, collapse = "_"),
        "__sim_", 
        paste0(sims.so.far, collapse = "_"), 
        ".csv"
      )
      
      df <- array2DF(metrics)
      colnames(df) <- c("metric", "n_bands", "simulation", "n", "method", "value")
      write.table(df, path.join("results", "partial", fname), row.names = F)
    }
    
    print("")
  }
}

df <- array2DF(metrics)
colnames(df) <- c("metric", "n_bands", "simulation", "n", "method", "value")
write.table(df, path.join("results", "simulation_lasso.csv"), row.names = F)

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