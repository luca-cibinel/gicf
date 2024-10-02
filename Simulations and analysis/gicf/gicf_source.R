library(Rcpp)
library(groupdata2)

sourceCpp("../gicf/gicf_source.cpp")

# GICF UTILITY ----
#

#' gicf.likelihood
#'
#' Computes the lasso penalized log-likelihood, namely
#' 
#' $$ -\frac{np}{2} - \frac{n}{2}\log{|\Sigma|} -\frac{n}{2}tr[\Sigma^{-1}S] - \lambda\frac{n}{2}\sum_{i \neq j}|\sigma_{ij}|. $$
#'
#' @param Sigma The matrix $\Sigma$
#' @param S The matrix $S$
#' @param n The number of data points
#' @param l the penalization coefficient $\lambda$
#'
#' @return The value of the $\ell_1-$penalized log.likelihood
#'
gicf.likelihood <- function(Sigma, S, n, l){
  P <- 1 - diag(nrow(S))
  
  out <- profileloglik(Sigma, S, n)
  
  return(out$loglik - 0.5*n*l*sum(abs(P * Sigma)))
}

#' max.lambda
#' 
#' Computes the maximum effective value of the penalization parameter $\lambda$ for the GICF algorithm
#' basing on the input sample covariance matrix and the initial estimate
#' 
#' @param S The sample covariance matrix of the data 
#' @param Sigma.init The initial estimate provided to GICF
#' @param adj an optional adjacency matrix
#'
#' @return A scalar representing the maximum effecting value of the penalization parameter $\lambda$ in GICF
#'
#' @examples
max.lambda <- function( S, Sigma.init = diag(diag(S)), adj = 1 - diag(1, nrow(S)) ){
  DSl <- diag(1 / diag(S))
  DSr <- diag(1 / pmin(diag(S), diag(Sigma.init)) )
  
  ML <- DSl %*% abs(S) %*% DSr
  diag(ML) <- 0
  ML <- ML * adj
  
  lambda.max <- max(ML[upper.tri(ML)])
  
  return(lambda.max)
}

#' max.kappa
#' 
#' Computes the maximum effective value of the penalization parameter $\kappa$ for the GICF algorithm
#' basing on the input sample covariance matrix and the parameter $\lambda$
#' 
#' @param S The sample covariance matrix of the data 
#' @param l The parameter $\lambda$
#' @param adj an optional adjacency matrix
#'
#' @return A scalar representing the maximum effecting value of the penalization parameter $\kappa$ in GICF
#'
#' @examples
max.kappa <- function( S, l, adj = 1 - diag(1, nrow(S)) ){
  p <- nrow(S)
  
  k.max <- -Inf
  
  Q <- diag(S)
  
  for(i in 1:(p - 1)){
    for(j in (i + 1):p){
      if(adj[i, j] == 1){
        g.ij <- (abs(S[i,j])/l - S[i,i]*Q[j])
        g.ij <- round(g.ij, 5) # for numerical instabilities
        #print(g.ij)
        
        if(g.ij >= 0){
          f.ij <- (S[i,i] + Q[j])/2
          
          k.ij <- -f.ij + sqrt(f.ij^2 + g.ij)
          
          if(k.ij > k.max)
            k.max <- k.ij
        }
      }
    }
  }
  
  return(abs(k.max))
}

# SPARSE ITERATIVE CONDITIONAL FITTING ----
# 
# Implementation of the SICF algorithm
#
# ************

#' Infers a path of Gaussian covariance graph models using the Sparse Iterative Conditional Fitting algorithm
#'
#' @param y The data
#' @param lambda A vector of values of lambda for which the inference has to be performed
#' @param max.iter The maximum number of iterations allowed (the default is 2500)
#' @param tol A tolerance below which numbers are considered 0 (the default is 1e-4)
#' @param Sigma.init An initial estimate. If NULL, the diagonal sample covariance matrix will be used (the default is NULL)
#' @param adj An adjacency matrix to brvided in case that a desired sparsity pattern is to be enforced (optional)
#'
#' @return A list in which each lambda (converted to a string) is paired to the output of the corresponding SICF call
#' @export
#'
#' @examples
gicf.run <- function(y = NULL, S = NULL, n = NULL, lambda = 0, 
                       max.iter = 2500, tol = 1e-4, Sigma.init = NULL, adj = NULL){
  if(is.null(y) && (is.null(S) || is.null(n)))
    stop("If y is not provided, BOTH S and n must be provided!")
  
  if(is.null(y)){
    p <- ncol(S)
  }else{
    n <- nrow(y) # Extract the values of n and p
    p <- ncol(y)
    S <- cov(y) * (n - 1)/n # Compute the sample covariance matrix (data are empirically centered)
  }
  
  if(is.null(adj)){
    adj <- matrix(1, p, p)
    diag(adj) <- 0
  }
  
  if(is.null(Sigma.init)) # Ensure that an initial estimate is given, if not pick diag(S)
    Sigma <- diag(diag(S))
  else
    Sigma <- Sigma.init
  
  lambda.max <- max.lambda(S, Sigma, adj)
  
  fit <- gicf_wrapper(Sigma, adj, n, S, lambda, lambda.max, tol, tol, max.iter, max.iter)
  
  if(length(lambda) == 1)
    return(fit$out[[1]])
  
  return(fit)
}