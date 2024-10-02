# HEADER ====
rm(list = ls()) # clear environment

library(groupdata2)
source("../gicf/gicf_source.R")

set.seed(1234)

use.banded.structure <- F

if(use.banded.structure){
  bands.rocks <- 17
  bands.metals <- 31
}else{
  bands.rocks <- 60
  bands.metals <- 60
}

N.lambdas <- 30
N.kappas <- 30
N.out.folds <- 5
N.in.folds <- 10

# DATA ====

if(!file.exists("Index")){ # If data not available, download them from UCI repository
  url <- "https://archive.ics.uci.edu/static/public/151/connectionist+bench+sonar+mines+vs+rocks.zip"
  
  temp.file <- tempfile()
  download.file(url, temp.file)
  unzip(temp.file, exdir = ".")
  unlink(temp.file)
}

data.sonar.full <- read.csv("sonar.all-data", header = F)

# Fold data for nested cross validation; only outer folding (remember to set seed for reproducibility)
data.metals <- data.frame(fold(data.sonar.full[data.sonar.full$V61 == "M",], k = N.out.folds))
data.rocks <- data.frame(fold(data.sonar.full[data.sonar.full$V61 == "R",], k = N.out.folds))



# FUNCTIONS ====

# Create a banded adjacency matrix
banded.adj <- function(p, k){
  adj <- matrix(0, p, p)
  
  for(i in 1:(p - 1)){
    for(j in (i + 1):p){
      if(j - i <= k) adj[i, j] <- 1
    }
  }
  
  return(adj + t(adj))
}

## GICF  utilities====

# Finds the maximum value of lambda for which kappa_max is not numerically zero
# Needed due to numerical issues in the sample covariance matrix
empirical.max.lambda <- function(S, adj = 1 - diag(1, nrow(S)), zero = 1e-8){
  l.max <- max.lambda(S, adj = adj)
  
  a <- 0
  b <- l.max
  
  while(b - a > zero){
    d <- (a + b)/2
    k <- max.kappa(S, d, adj)
    
    if(k > zero){
      a <- d
    }else{
      b <- d
    }
  }
  
  return((a + b)/2)
}

model.selection.cv <- function(y, # data
                               N.folds, # num. of CV folds
                               N.l,
                               N.k,
                               l.max = Inf, # max value of \lambda to be explored
                               k.max = 10, # max value of \kappa to be explored
                               adj = NULL# prespecified adjacency matrix
){
  n <- nrow(y)
  p <- ncol(y)
  
  if(is.null(adj))
    adj <- 1 - diag(1, p)
  
  S <- cov(y) * (n - 1)/n
  
  folded.data <- fold(data.frame(y), k = N.folds)
  
  lambda.max <- min( empirical.max.lambda(S, adj = adj), l.max )
  lambda.seq <- seq(0, lambda.max, length.out = N.l)
  
  history <- NULL
  
  for(L in 1:N.l){
    lambda.loc <- lambda.seq[L]
    
    kappa.max.loc <- min( max.kappa(S, lambda.loc, adj = adj), k.max )
    
    if(n > p)
      kappa.seq <- seq(0, kappa.max.loc, length.out = N.k)
    else
      kappa.seq <- seq(0, kappa.max.loc, length.out = N.k + 1)[-1]
    
    #print(paste(kappa.max.loc, seq.length.loc))
    
    for(K in 1:N.k){
      kappa.loc <- kappa.seq[K]
      val <- 0
      
      for(fold in 1:N.folds){
        #print(paste(fold, kappa.loc, lambda.loc))
        train <- folded.data[folded.data$.folds != fold, ]
        validation <- folded.data[folded.data$.folds == fold, ]
        
        train <- unname(as.matrix(
          train[, names(train) != ".folds"]
        ))
        
        validation <- unname(as.matrix(
          validation[, names(validation) != ".folds"]
        ))
        
        n.train <- nrow(train)
        n.validation <- nrow(validation)
        
        S.train <- cov(train) * (n.train - 1)/n.train
        S.validation <- cov(validation) * (n.validation - 1)/n.validation
        
        D <- diag(kappa.loc, p)
        fit <- gicf.run(S = S.train + D, n = n.train, lambda = lambda.loc, adj = adj)$sigma
        
        val <- val + gicf.likelihood(fit, S.validation, n.validation, 0)
      }
      
      history <- rbind(history, c(kappa.loc, lambda.loc, val))
    }
  }
  
  colnames(history) <- c("K", "L", "val")
  return(history)
}

nested.model.selection.cv <- function(y.folded, # folded data
                               N.in.folds, # num. of inner CV folds
                               N.l,
                               N.k,
                               l.max = Inf, # max value of \lambda to be explored
                               k.max = 10, # max value of \kappa to be explored
                               adj = NULL# prespecified adjacency matrix
){
  N.out.folds <- length(unique(y.folded$.folds))
  history <- array(0, c(N.l*N.k, 3, N.out.folds))
  
  for(fold.out in 1:N.out.folds){
    trainval <- y.folded[y.folded$.folds != fold.out, -c(61, 62)]
    
    history.in <- model.selection.cv(
      trainval,
      N.in.folds,
      N.l,
      N.k,
      l.max,
      k.max,
      adj
    )
    
    history[,,fold.out] <- history.in
  }
  
  
  dimnames(history)[[2]] <- c("K", "L", "val")
  return(history)
}


## Quadratic discriminant analysis ====

# kappa (lambda) is a vector kappa = c(kappa.rock, kappa.metal)
G <- function(x, pi.rocks, mu.rocks, omega.rocks, pi.metals, mu.metals, omega.metals){
  val.rocks <- 0.5*(log(det(omega.rocks)) - 
                     t(x - mu.rocks) %*% omega.rocks %*% (x - mu.rocks)) + 
              log(pi.rocks)
  
  val.metal <- 0.5*(log(det(omega.metals)) - 
                      t(x - mu.metals) %*% omega.metals %*% (x - mu.metals)) + 
              log(pi.metals)
  
  return(ifelse(val.rocks > val.metal, "R", "M"))
}

# Quadratic discriminant analysis (compute mean error)
# kappa (lambda) as in G
QDA.err <- function(data.test, pi.rocks, mu.rocks, omega.rocks, pi.metals, mu.metals, omega.metals){
  errors <- 0
  
  for(i in 1:nrow(data.test)){
    est <- G(
      x = unlist(data.test[i, -61]),
      pi.rocks,
      mu.rocks,
      omega.rocks,
      pi.metals,
      mu.metals,
      omega.metals
    )[1,1]
    
    if(est != data.test[i, 61]) errors <- errors + 1
  }
  
  return(errors/nrow(data.test))
}

QDA.cv <-  function(y.folded, cv.scores.rocks, cv.scores.metals, use.k, use.l){
  N.folds <- length(unique(y.folded$.folds))
  
  err <- 0
  
  for(fold in 1:N.folds){
    train <- y.folded[y.folded$.folds != fold, colnames(y.folded) != ".folds"]
    test <- y.folded[y.folded$.folds == fold, colnames(y.folded) != ".folds"]
    
    pi.rocks <- sum(train$V61 == "R")/nrow(train)
    mu.rocks <- colMeans(train[train$V61 == "R", -61])
    
    pi.metals <- sum(train$V61 == "M")/nrow(train)
    mu.metals <- colMeans(train[train$V61 == "M", -61])
    
    cv.scores.rocks.loc <- cv.scores.rocks[,,fold]
    cv.scores.metals.loc <- cv.scores.metals[,,fold]
    
    if(!use.k){
      cv.scores.rocks.loc <- cv.scores.rocks.loc[
        cv.scores.rocks.loc[,"K"] == 0,
      ]
      cv.scores.metals.loc <- cv.scores.metals.loc[
        cv.scores.metals.loc[,"K"] == 0,
      ]
    }
    
    if(!use.l){
      cv.scores.rocks.loc <- cv.scores.rocks.loc[
        cv.scores.rocks.loc[,"L"] == 0,
      ]
      cv.scores.metals.loc <- cv.scores.metals.loc[
        cv.scores.metals.loc[,"L"] == 0,
      ]
    }
    
    if(is.null(dim(cv.scores.rocks.loc))){
      cv.scores.rocks.loc <- matrix(cv.scores.rocks.loc, 1, 3)
      colnames(cv.scores.rocks.loc) <- c("K", "L", "val")
    }
    
    if(is.null(dim(cv.scores.metals.loc))){
      cv.scores.metals.loc <- matrix(cv.scores.metals.loc, 1, 3)
      colnames(cv.scores.metals.loc) <- c("K", "L", "val")
    }
    
    best.rocks <- which.max(cv.scores.rocks.loc[, "val"])
    best.metals <- which.max(cv.scores.metals.loc[, "val"])
    
    pars.rocks <- cv.scores.rocks.loc[best.rocks, -3]
    pars.metals <- cv.scores.metals.loc[best.metals, -3]
    
    n.rocks <- sum(train$V61 == "R")
    S.rocks <- cov(train[train$V61 == "R", -61]) * (n.rocks - 1)/n.rocks
    sigma.rocks <- gicf.run(S = S.rocks + diag(pars.rocks["K"], 60), n = n.rocks,
                            lambda = pars.rocks["L"],
                            adj = banded.adj(60, bands.rocks))$sigma
    
    n.metals <- sum(train$V61 == "M")
    S.metals <- cov(train[train$V61 == "M", -61]) * (n.metals - 1)/n.metals
    sigma.metals <- gicf.run(S = S.metals + diag(pars.metals["K"], 60), n = n.metals,
                            lambda = pars.metals["L"],
                            adj = banded.adj(60, bands.metals))$sigma
    
    err <- err + QDA.err(test, pi.rocks, mu.rocks, solve(sigma.rocks), pi.metals, mu.metals, solve(sigma.metals))
  }
  
  return(err/N.folds)
}

# MODEL SELECTION ====


## Model selection (rocks) ====
print("Computing CV scores for ROCKS...")

cv.scores.rocks.lambda <- nested.model.selection.cv(data.rocks,
                                             N.in.folds,
                                             N.lambdas,
                                             1,
                                             k.max = 0,
                                             adj = banded.adj(60, bands.rocks))
cv.scores.rocks.kappa <- nested.model.selection.cv(data.rocks,
                                             N.in.folds,
                                             1,
                                             N.kappas,
                                             l.max = 0,
                                             adj = banded.adj(60, bands.rocks))
cv.scores.rocks <- nested.model.selection.cv(data.rocks,
                                      N.in.folds,
                                      N.lambdas,
                                      N.kappas,
                                      l.max = 20,
                                      k.max = 6e-4,
                                      adj = banded.adj(60, bands.rocks)) 

cv.scores.rocks <- abind::abind(cv.scores.rocks, 
                                cv.scores.rocks.kappa,
                                cv.scores.rocks.lambda,
                                along = 1)

## Model selection (metals) ====
print("Computing CV scores for METALS...")
cv.scores.metals.lambda <- nested.model.selection.cv(data.metals,
                                             N.in.folds,
                                             N.lambdas,
                                             1,
                                             k.max = 0,
                                             adj = banded.adj(60, bands.metals))
cv.scores.metals.kappa <- nested.model.selection.cv(data.metals,
                                            N.in.folds,
                                            1,
                                            N.kappas,
                                            l.max = 0,
                                            adj = banded.adj(60, bands.metals))
cv.scores.metals <- nested.model.selection.cv(data.metals,
                                       N.in.folds,
                                       N.lambdas,
                                       N.kappas,
                                       l.max = 20,
                                       k.max = 6e-4,
                                       adj = banded.adj(60, bands.metals)) 


cv.scores.metals <- abind::abind(cv.scores.metals, 
                                cv.scores.metals.kappa,
                                cv.scores.metals.lambda,
                                along = 1)

# OUTPUT ====
source("../gicf/gicf_source.R")

data.folded <- rbind(data.rocks, data.metals)
error.rate.mle <- QDA.cv(data.folded, 
                         cv.scores.rocks, cv.scores.metals,
                         F, F)
error.rate.lambda <- QDA.cv(data.folded, 
                         cv.scores.rocks, cv.scores.metals,
                         F, T)
error.rate.kappa <- QDA.cv(data.folded, 
                         cv.scores.rocks, cv.scores.metals,
                         T, F)
error.rate.gicf <- QDA.cv(data.folded, 
                         cv.scores.rocks, cv.scores.metals,
                         T, T)

print(error.rate.mle)
print(error.rate.lambda)
print(error.rate.kappa)
print(error.rate.gicf)