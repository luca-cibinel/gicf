source("../gicf/gicf.R")
library(groupdata2)

#==========================
#
# FUNCTIONS 


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
  l.max <- lambdamax(S, adj = adj)
  
  a <- 0
  b <- l.max
  
  while(b - a > zero){
    d <- (a + b)/2
    k <- kappamax(S, d, adj)
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
                               k.max = 5, # max value of \kappa to be explored
                               adj = NULL,# prespecified adjacency matrix
                               s.num = s.num
){
  n <- nrow(y)
  p <- ncol(y)
  
  if(is.null(adj))
    adj <- 1 - diag(1, p)
  
  S <- cov(y) * (n - 1)/n
  
  set.seed(s.num)
  folded.data <- fold(data.frame(y), k = N.folds)
  
  #lambda.max <- min(lambdamax(S, adj = adj), l.max)
  lambda.max <- min(empirical.max.lambda(S, adj = adj), l.max)
  lambda.seq <- seq(0, lambda.max, length.out = N.l)
  history <- NULL
  
  for(L in 1:N.l){
    lambda.loc <- lambda.seq[L]
    
    kappa.max.loc <- min(kappamax(S, lambda.loc, adj = adj), k.max)
    #
    
    # if(n > p)
    #    kappa.seq <- seq(0,log(kappa.max.loc+1), length.out = N.k)
    #  else
    #    kappa.seq <- seq(0, log(kappa.max.loc+1), length.out = N.k + 1)[-1]
    kappa.seq <- seq(0,log(kappa.max.loc+1), length.out = N.k)
    kappa.seq <- exp(kappa.seq) - 1
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
        
        fit <- gicf(S = S.train, n = n.train, lambda = lambda.loc, kappa = kappa.loc, adj = adj)$sigma
        
        val <- val + gcgmloglik(fit, S.validation, n.validation)
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
                                      k.max = 5, # max value of \kappa to be explored
                                      adj = NULL,# prespecified adjacency matrix
                                      s.num = s.num
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
      adj,# prespecified adjacency matrix
      s.num = s.num
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
    sigma.rocks <- gicf(S = S.rocks, n = n.rocks,
                        lambda = pars.rocks["L"],
                        kappa = pars.rocks["K"],
                        adj = banded.adj(60, bands.rocks))$sigma
    
    n.metals <- sum(train$V61 == "M")
    S.metals <- cov(train[train$V61 == "M", -61]) * (n.metals - 1)/n.metals
    sigma.metals <- gicf(S = S.metals, n = n.metals,
                         lambda = pars.metals["L"],
                         kappa = pars.metals["K"],
                         adj = banded.adj(60, bands.metals))$sigma
    
    err <- err + QDA.err(test, pi.rocks, mu.rocks, solve(sigma.rocks), pi.metals, mu.metals, solve(sigma.metals))
  }
  
  return(err/N.folds)
}

# END FUNCTIONS 
#
#============================================



N.out.folds <- 5
N.in.folds <- 5


# DATA ====

if(!file.exists("Index")){ # If data not available, download them from UCI repository
  url <- "https://archive.ics.uci.edu/static/public/151/connectionist+bench+sonar+mines+vs+rocks.zip"
  
  temp.file <- tempfile()
  download.file(url, temp.file)
  unzip(temp.file, exdir = ".")
  unlink(temp.file)
}

data.sonar.full <- read.csv("sonar.all-data", header = F)
#===========

unst.result <- matrix(NA, 4, 10)
diag.result <- matrix(NA, 4, 10)

for (i in 1:10){
  s.num <- 100+i*10-10
  print(paste0("Iteration n. ", i))
  print(paste0("seed: ", s.num))
  #
  # Fold data for nested cross validation; only outer folding (remember to set seed for reproducibility)
  set.seed(s.num)
  data.metals <- data.frame(fold(data.sonar.full[data.sonar.full$V61 == "M",], k = N.out.folds))
  data.rocks <- data.frame(fold(data.sonar.full[data.sonar.full$V61 == "R",], k = N.out.folds))
  #  
  # unstructured
  use.banded.structure <- F
  N.lambdas <- 15
  N.kappas  <- 15
  source("sub-sonar-analysis.R")
  unst.result[,i] <- errors
  #names(errors) <- c("MLE", "LAMBDA", "KAPPA", "GICF")
  #write.table(errors, paste0("sonar_analysis_band_structure_", use.banded.structure, i, ".csv"))
  
  # diagonal
  use.banded.structure <- T
  N.lambdas <- 15
  N.kappas  <- 30
  source("sub-sonar-analysis.R")
  diag.result[,i] <- errors
  #names(errors) <- c("MLE", "LAMBDA", "KAPPA", "GICF")
  #write.table(errors, paste0("sonar_analysis_band_structure_", use.banded.structure, i, ".csv"))
  
  print(unst.result)
  print(diag.result)
}
row.names(unst.result) <- c("MLE", "LAMBDA", "KAPPA", "GICF")
row.names(diag.result) <- c("MLE", "LAMBDA", "KAPPA", "GICF")
write.table(unst.result, "sonar_analysis_unstractured.csv", col.names = FALSE)
write.table(diag.result, "sonar_analysis_diagonal.csv", col.names = FALSE)

apply(unst.result, 1, "mean")
apply(unst.result, 1, "sd")

apply(unst.result, 1, "mean")-2*apply(unst.result, 1, "sd")
apply(unst.result, 1, "mean")+2*apply(unst.result, 1, "sd")


apply(diag.result, 1, "mean")
apply(diag.result, 1, "sd")


apply(diag.result, 1, "mean")-2*apply(diag.result, 1, "sd")
apply(diag.result, 1, "mean")+2*apply(diag.result, 1, "sd")



main <- "comparison"
res.tmp <- unst.result
delta.mle <- res.tmp[1,]-res.tmp[4,]
mean(delta.mle)
sd(delta.mle)

delta.l <- res.tmp[2,]-res.tmp[4,]
mean(delta.l)
sd(delta.l)

delta.k  <- res.tmp[3,]-res.tmp[4,]   
mean(delta.k)
sd(delta.k)


delta.diag <- diag.result[4,]-res.tmp[4,]
boxplot(delta.mle, delta.l, delta.k, delta.diag, main=main, names=c("U-MLE", "U-covglasso", "U-ridge", "D-rrcovglasso"), ylim=c(-0.05, 0.30), pch=20)
abline(h=0, lty=2)

