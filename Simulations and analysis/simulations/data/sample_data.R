library(mvtnorm)

# [1] "1 0.04"
# [1] "2 0.0791836734693878"
# [1] "3 0.117551020408163"
# [1] "4 0.155102040816327"
# [1] "5 0.191836734693878"
# [1] "6 0.227755102040816"
# [1] "7 0.262857142857143"
# [1] "8 0.297142857142857"
# [1] "9 0.330612244897959"
# [1] "10 0.363265306122449"
# [1] "11 0.395102040816327"
# [1] "12 0.426122448979592"
# [1] "13 0.456326530612245"
# [1] "14 0.485714285714286"
# [1] "15 0.514285714285714"
# [1] "16 0.542040816326531"
# [1] "17 0.568979591836735"
# [1] "18 0.595102040816327"
# [1] "19 0.620408163265306"
# [1] "20 0.644897959183673"
# [1] "21 0.668571428571429"
# [1] "22 0.691428571428571"
# [1] "23 0.713469387755102"
# [1] "24 0.73469387755102"
# [1] "25 0.755102040816326"
# [1] "26 0.77469387755102"
# [1] "27 0.793469387755102"
# [1] "28 0.811428571428571"
# [1] "29 0.828571428571429"
# [1] "30 0.844897959183673"
# [1] "31 0.860408163265306"
# [1] "32 0.875102040816326"
# [1] "33 0.888979591836735"
# [1] "34 0.902040816326531"
# [1] "35 0.914285714285714"
# [1] "36 0.925714285714286"
# [1] "37 0.936326530612245"
# [1] "38 0.946122448979592"
# [1] "39 0.955102040816327"
# [1] "40 0.963265306122449"
# [1] "41 0.970612244897959"
# [1] "42 0.977142857142857"
# [1] "43 0.982857142857143"
# [1] "44 0.987755102040816"
# [1] "45 0.991836734693878"
# [1] "46 0.995102040816327"
# [1] "47 0.997551020408163"
# [1] "48 0.999183673469388"
# [1] "49 1"
# [1] "50 1"

path.join <-  function(...){
  sep <- .Platform$file.sep
  pat <- paste0("\\", sep, "$")
  
  paths <- lapply(list(...), function(x){gsub(pat, "", x)})
  
  return(paste0(paths, collapse = sep))
}

density <- function(p, b){
  n.tot.edges <- 0.5*p*(p-1)
  n.edges <- sum(p - 1:b)
  
  return(n.edges/n.tot.edges)
}

rsigma <- function(p, b){
  Sigma <- matrix(0, p, p)
  
  band.idxs <- NULL
  for(i in 1:b){
    band.idxs <- rbind(band.idxs, cbind( 1:(p - i), (i + 1):p ))
  }
  
  Sigma[band.idxs] <- 2*rbinom(sum(p - 1:b), 1, 0.5) - 1
  Sigma <- Sigma + t(Sigma)
  
  eig <- range(eigen(Sigma)$values)
  diag(Sigma) <- (eig[2] - p*eig[1])/(p - 1)
  
  return(Sigma)
}

# Data simulation ====
n <- 2000
N.sim <- 50
bands <- c(3, 10)
ps <- c(25, 50, 100, 200, 300)

for(p in ps){
  model.directory <- path.join(".", paste0("p", p))
  
  for(b in bands){
    Sigma <- rsigma(p, b)
    band.directory <- path.join(model.directory, paste0("b", b))
    
    if(!file.exists(band.directory)){
      dir.create(band.directory, recursive = T)
    }
    
    write.table(Sigma, path.join(band.directory, "sigma.dat"))
    
    for(s in 1:N.sim){
      data <- as.matrix(rmvnorm(n, sigma = Sigma))
      
      write.table(data, path.join(band.directory, paste0("sim", s, ".dat")))
    }
  }
}