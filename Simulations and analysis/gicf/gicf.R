library(Rcpp)
library(groupdata2)

sourceCpp("../gicf/gicf.cpp")

#' Maximum effective range of regularisation parameters
#'
#' Compute the effective range of the regularisation paramerters \eqn{\kappa} and \eqn{\lambda}.
#'
#' These utility functions describe the boundary of the region
#' \deqn{\mathcal{H} = \{(\kappa, \lambda) \in \mathbb{R}_{\geq 0}^2: \lambda \leq \lambda_{MAX}(\kappa)\},}
#' with
#' \deqn{\lambda \leq \lambda_{MAX}(\kappa) \Longleftrightarrow \kappa \leq \kappa_{MAX}(\lambda),}
#' \deqn{\lambda_{MAX}(\kappa) = \max_{{i,j}\in\mathcal{G}}\frac{|s_{ij}|}{(s_{ii} + \kappa)(s_{jj} + \kappa)},}
#' \deqn{\kappa_{MAX}(\lambda) = \max_{\substack{{i,j}\in\mathcal{G}\\ g_{ij}(\lambda)\geq 0}}\left\{\sqrt{\frac{1}{4}(s_{ii} + s_{jj}) + g_{ij}(\lambda)} - \frac{1}{2}(s_{ii}+s_{jj})\right\},}
#' \deqn{g_{ij}(\lambda) = \frac{|s_{ij}|}{\lambda} - s_{ii}s_{jj}.}
#' Here \eqn{S} is the sample covariance matrix and \eqn{\mathcal{G}} is a graph whose adjacency matrix has the same sparsity pattern as \code{adj}.
#' If the parameters \eqn{(\kappa, \lambda)} lay outside of \eqn{\mathcal{H}}, and the starting
#' point of the Generalised Iterative Conditional Fitting algorithm is \eqn{\text{diag}(S + \kappa I)},
#' then the output will also be \eqn{\text{diag}(S + \kappa I)}.
#'
#'
#' @param S The sample covariance matrix.
#' @param kappa,lambda The non-negative ridge regularisation/lasso shrinkage parameters.
#' @param adj An optional matrix whose pattern of zeroes is to be enforced
#'  onto the final output of the Generalised Iterative Conditional Fitting algorithm.
#'
#' @returns \code{lambdamax} returns a scalar value representing \eqn{\lambda_{MAX}(\kappa)}. \code{kappamax} returns a scalar value representing \eqn{\kappa_{MAX}(\lambda)}
#'
#' @example inst/examples/example_hyperparameters.R
#'
#' @name Hyperparameters
#'
#' @export
lambdamax <- function( S, kappa = 0, adj = 1 - diag(1, nrow(S)) ){
  adj <- abs(sign(adj)) # Ensure that adj has the right format
  diag(adj) <- 0

  DS <- diag(1 / diag(S + kappa))

  ML <- DS %*% abs(S) %*% DS
  diag(ML) <- 0
  ML <- ML * adj

  lambda.max <- max(ML[upper.tri(ML)])

  return(lambda.max)
}

#' @rdname Hyperparameters
#' @export
kappamax <- function( S, lambda, adj = 1 - diag(1, nrow(S)) ){
  adj <- abs(sign(adj)) # Ensure that adj has the right format
  diag(adj) <- 0

  p <- nrow(S)

  k.max <- -Inf

  Q <- diag(S)

  for(i in 1:(p - 1)){
    for(j in (i + 1):p){
      if(adj[i, j] == 1){
        g.ij <- (abs(S[i,j])/lambda - S[i,i]*Q[j])
        g.ij <- round(g.ij, 5) # for numerical instabilities

        if(g.ij >= 0){
          f.ij <- (S[i,i] + Q[j])/2

          k.ij <- -f.ij + sqrt(f.ij^2 + g.ij)

          if(k.ij > k.max)
            k.max <- k.ij
        }
      }
    }
  }

  return(unname(abs(k.max)))
}

#' Gaussian Covariance Graphical Model Loglikelihood function
#'
#' Computes the penalised loglikelihood function of a Gaussian covariance graph model.
#'
#' When imposing sparsity on the covariance matrix of a multivariate Gaussian distribution, the resulting
#' model can be interpreted as a covariance graphical model, i.e., the independence structure of the components
#' of the random vector can be encoded by a graph in which the nodes are identified with the variables and
#' a missing edge between two nodes implies that the corresponding variables are marginally independent.
#'
#' In particular, this model admits both a ridge and a lasso penalty, resulting in the loglikelihood function
#' \deqn{-\text{log}|\Sigma| - \text{trace}(\Sigma^{-1}S) - \lambda\|\Sigma - \text{diag}(\Sigma)\|_1 - \kappa\|\Sigma^{-1}\|_1,}
#' where \eqn{\lambda, \kappa \geq 0}.
#'
#' @param Sigma The covariance matrix.
#' @param S The sample covariance matrix.
#' @param n The size of the observed dataset.
#' @param lambda A non-negative lasso parameter.
#' @param kappa A non-negative ridge regularisation parameter.
#'
#' @returns The value of the penalised loglikelihood function.
#' @export
#'
#' @example inst/examples/example_loglik.R
gcgmloglik <- function(Sigma, S, n, lambda = 0, kappa = 0){
  P <- 1 - diag(nrow(S))

  out <- profileloglik(Sigma, S + diag(kappa, nrow(S)), n)

  return(out$loglik - 0.5*n*lambda*sum(abs(P * Sigma)))
}

#' Penalised maximum likelihood covariance matrix estimation
#'
#' Estimation of a sparse covariance matrix via
#' the ridge-regularised covglasso estimator described in Cibinel et al. (2024).
#'
#' This function computes the ridge-regularised covglasso estimator
#' of the covariance matrix of a multivariate normal distribution, that is
#' it computes the maximum of the penalised log-likelihood
#' \deqn{-\text{log}|\Sigma| - \text{trace}(\Sigma^{-1}S) - \lambda\|\Sigma - \text{diag}(\Sigma)\|_1 - \kappa\|\Sigma^{-1}\|_1,}
#' where \eqn{\lambda, \kappa \geq 0}.
#' The optimum is computed via a coordinate descent algorithm, resulting
#' in an approach which unifies and extends the methods of Chaudhuri et. al
#' (2007), Warton (2008), Bien and Tibshirani (2011) and Wang (2014).
#'
#' @param data A numerical matrix whose rows contain
#'  the observations of multivariate normal random vector.
#'  If \code{NULL}, the sample covariance matrix S and the dataset size n
#'  must be provided.
#' @param S The sample covariance matrix. Must be provided if data is \code{NULL}.
#' @param n The dataset size. Must be provided if data is \code{NULL}.
#' @param lambda A vector of non-negative lasso parameters. For efficency purposes,
#'  should be sorted from largest to smallest.
#' @param kappa A non-negative ridge regularisation parameter.
#' @param max.iter The maximum number of iterations
#'  allowed for the coordinate descent algorithm.
#' @param tol A numerical tolerance below which quantities are treated as zero.
#' @param Sigma.init The initial guess for the coordinate descent algorithm.
#'  Defaults to the diagonal of the sample covariance matrix.
#' @param adj An optional matrix whose pattern of zeroes is enforced
#'  onto the final output of the algorithm.
#'
#' @returns If a scalar value for \code{lambda} is provided, a list containing the following elements.
#' \tabular{ll}{
#'  \code{sigma} \tab The estimate of the covariance matrix. \cr\tab\cr
#'  \code{omega} \tab The inverse of the estimated covariance matrix. \cr\tab\cr
#'  \code{loglik} \tab The (unpenalised) log-likelihood at the optimum. \cr\tab\cr
#'  \code{loglikpen} \tab The (penalised) log-likelihood at the optimum. \cr\tab\cr
#'  \code{it} \tab The number of iterations needed to reach convergence.\cr\tab\cr
#' }
#' If a vector of values of \code{lambda} is provided, the output is
#' a list in which each entry is itself a list, structured as above,
#' associated with the corresponding value of \code{lambda}.
#'
#' @section References:
#'
#' Chaudhuri, S., M. Drton, and T. S. Richardson (2007). Estimation of a covariance matrix with zeros. Biometrika 94 (1), 199–216.
#'
#' Cibinel, L., A. Roverato, and V. Vinciotti (2024). A unified approach to penalized likelihood estimation of covariance matrices in high dimensions. arXiv, arXiv:2410.02403.
#'
#' Bien, J. and R. J. Tibshirani (2011). Sparse estimation of a covariance matrix. Biometrika 98 (4), 807–820.
#'
#' Wang, H. (2014). Coordinate descent algorithm for covariance graphical lasso. Statistics and Computing 24, 521–529.
#'
#' Warton, D. I. (2008). Penalized normal likelihood and ridge regularization of correlation and covariance matrices. Journal of the American Statistical Association 103 (481), 340–349.
#'
#' @name gicf
#' @export
#'
#' @example inst/examples/example_gicf.R
gicf <- function(data = NULL, S = NULL, n = NULL, lambda = 0, kappa = 0,
                    max.iter = 2500, tol = 1e-4, Sigma.init = NULL, adj = NULL){
  if(is.null(data) && (is.null(S) || is.null(n)))
    stop("If y is not provided, BOTH S and n must be provided!")

  if(is.null(data)){
    p <- ncol(S)
  }else{
    n <- nrow(data) # Extract the values of n and p
    p <- ncol(data)
    S <- cov(data) * (n - 1)/n # Compute the sample covariance matrix (data are empirically centered)
  }

  S <- S + diag(kappa, p)

  if(is.null(adj)){
    adj <- matrix(1, p, p)
    diag(adj) <- 0
  }else{
    adj <- abs(sign(adj)) # Ensure that adj has the right format
    diag(adj) <- 0
  }

  if(is.null(Sigma.init)){ # Ensure that an initial estimate is given, if not pick diag(S)
    Sigma <- diag(diag(S))
    lambda.max <- lambdamax(S, adj)
  }else{
    Sigma <- Sigma.init
    lambda.max <- -1 # Initial condition cannot be guaranteed to be diag(S): disable lambdamax
  }
  
  if(!is.matrix(lambda)){ # If lambda is a scalar
    if(length(lambda) > 1){
      stop("Lambda can only be a scalar or a matrix.")
    }
    lambda = matrix(lambda, p, p)
    diag(lambda) <- 0
  }else{
    if(!isSymmetric(lambda)){
      stop("Provided lambda matrix is not symmetric.")
    }
  }
  
  fit <- gicf_core(Sigma, adj, n, S, lambda, tol, tol, max.iter, max.iter)
  
  return(fit)
}
