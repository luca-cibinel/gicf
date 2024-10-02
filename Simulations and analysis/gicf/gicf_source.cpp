/*
 * Based on the C++ code of covglasso and the R code of ICF (in the ggm package)
 */

#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List profileloglik(arma::mat sigma, arma::mat S, int n) {
  int p = sigma.n_rows;
  
  arma::mat inv = inv_sympd(sigma);
  double val;
  double sign;
  log_det(val, sign, sigma);
  val = -n*p*0.5*log(2*M_PI) -n*0.5*val -n*0.5*trace( inv*S ) ;
  
  return Rcpp::List::create( Named("loglik") = val,
                             Named("inv") =  inv);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List gicf_core(arma::mat start, arma::umat adj, 
                     int n, arma::mat S, double lambda, 
                     double lambda_max, 
                     double tolout, double tolin, 
                     double iterout, double iterin) {
  int p = S.n_rows;
  
  // Utilities
  arma::uvec ind = linspace<uvec>(0, p-1, p); // indexes from 1 to p
  arma::umat zr(1, 1); // a utility 1x1 matrix containing 0 (used to index a single row)
  zr(0, 0) = 0;
  
  // initial quantities
  arma::mat sigma = start;
  arma::mat sigmaprev = start;
  arma::mat inv = inv_sympd(sigma);
  
  // iteration utility quantities
  bool crit = true;
  double err;
  int it = 0;
  
  // check trivial cases
  if( start.is_diagmat() && (lambda > lambda_max - tolout) ){
    crit = false;
    sigma = arma::diagmat(S.diag(0));
  }
  
  if( (lambda < tolout) && (accu(adj) >= p*p - p) ){
    crit = false;
    sigma = S;
  }
  
  // main iteration cycle
  while ( crit ) {
    
    // iterate over the variables
    for ( int i = 0; i < p; i++ ) {
      //inner iteration utility quantities
      arma::urowvec row_i = adj.row(i);   // i-th row of adjacency matrix
      
      arma::uvec m_v = find(ind != i);    // elements [-i]
      arma::uvec v = find(ind == i);      // element [i]
      
      arma::mat beta(p - 1, 1, arma::fill::zeros); // place holder for \beta
      arma::mat beta_inv_v_beta(1,1, arma::fill::zeros); //place holder for \beta^T(\Sigma_{-i,-i})^{-1}\beta
      
      double tau;                         // place holder for \tau
      
      arma::mat inv_v = inv(m_v, m_v) - 
        inv.submat(m_v, v)*inv.submat(m_v, v).t() / as_scalar( inv(v, v) ); // (Sigma_{-i,-i})^{-1}
      
      // If there are neighbors
      if(arma::accu(row_i) > 0){
        // select neighborhood
        arma::uvec ne_v0 = find(row_i == 1);
        arma::uvec tmp0 = find(ne_v0 > i);
        
        // adapt neighbors indexes (since we consider the -i submatrices, indices above i have to shifted by -1)
        arma::uvec ne_v = ne_v0;
        ne_v(tmp0) = ne_v(tmp0) - 1;
        
        // extract neighbor related quantities
        arma::mat inv_v_ne = inv_v.cols(ne_v);
        
        // compute linear regression utility matrices
        arma::mat ZZ = inv_v_ne.t() * S(m_v, m_v) * inv_v_ne; 
        arma::mat ZY = inv_v_ne.t() * S.submat(m_v, v);
        
        // local linear regression phase
        arma::mat beta_loc = sigma.submat(ne_v0, v); // initial condition
        int q = beta_loc.n_rows; // n. of coefficients
        
        // linreg iteration utility quantities
        bool critin = true;
        double errin;
        int itin = 0;
        double tauprev = 0;
        arma::mat betaprev = beta_loc;
        
        // linreg iteration cycle
        while ( critin ) {
          // inner linreg utility quanitities
          arma::mat ZZbeta = ZZ * beta_loc;
          
          // update for tau
          tau = as_scalar( beta_loc.t()*ZZbeta - 2*beta_loc.t()*ZY + S(v,v) );
          
          // update for beta
          for ( int j = 0; j < q; j++ ) {
            double betajprev = as_scalar(beta_loc(j));
            
            double tmp = ZZ(j, j);
            double tmp2 = ZY(j) - ZZbeta(j) + tmp * beta_loc(j); // soft-thresh argument
            
            double sgn;
            if (tmp2 < 0) sgn = -tau/tmp; else if (tmp > 0) sgn = tau/tmp; else sgn = 0;
            
            beta_loc(j) = max( 0.0, abs(tmp2/tau) - lambda ) * sgn;
            
            double betadiff = beta_loc(j) - betajprev;
            if ( betadiff != 0 ) ZZbeta = ZZbeta + betadiff * ZZ.col(j);
          }
          
          errin = (as_scalar(accu( abs(beta_loc - betaprev) )) + abs(tau - tauprev))/p;
          itin++;
          critin = ( (errin > tolin) & (itin < iterin) );
          betaprev = beta_loc;
          tauprev = tau;
        }
        
        // fill non neighbor related coefficients with zeros
        beta(ne_v) = beta_loc;
        
        // Update inverse matrix
        arma::mat inv_v_beta = inv_v*beta;
        beta_inv_v_beta = beta.t()*inv_v_beta;
        
        inv(m_v, m_v) = inv_v + ( inv_v_beta*inv_v_beta.t() )/tau;
        inv(m_v, v) = -inv_v_beta/tau;
        inv(v, m_v) = (-inv_v_beta.t())/tau;
      }else{
        //If there are no neighbors, simply update inverse matrix and set \tau = S_{ii}
        
        inv(m_v, m_v) = inv_v;
        inv(m_v, v) = arma::mat(p - 1, 1, arma::fill::zeros); // if ne not empty
        inv(v, m_v) = arma::mat(1, p - 1, arma::fill::zeros); // if ne not empty
        tau = as_scalar(S(v, v));
      }
      
      // update \Sigma using \beta and \tau
      sigma(m_v, v) = beta;
      sigma(v, m_v) = beta.t();
      sigma(v, v) = tau + beta_inv_v_beta;
      
      // complete inverse matrix update
      arma::mat tmp(1,1);
      tmp(0,0) = 1/tau;
      inv(v, v) = tmp;
    }
    
    // check main cycle iteration conditions
    err = accu( abs( sigma - sigmaprev ) )/(p*p);
    sigmaprev = sigma;
    it++;
    crit = ( (err > tolout) & (it < iterout) );
  }
  
  // compute final log likelihood value
  Rcpp::List out = profileloglik(sigma, S, n);
  double llk = out["loglik"];
  double pen = 0.5*n*accu( abs(
    lambda * ( sigma - arma::diagmat(sigma.diag(0)) )
                      ) );
  
  return Rcpp::List::create( Named("omega") = inv,
                             Named("sigma") = sigma,
                             Named("loglik") = llk,
                             Named("loglikpen") = llk - pen,
                             Named("pen") = pen,
                             Named("it") = it,
                             Named("err") = err );
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List gicf_wrapper(arma::mat start, 
                      arma::umat adj, int n, arma::mat S,
                      arma::vec lambda, double lambda_max, 
                      double tolout, double tolin, double iterout, double iterin) {
  int L = lambda.n_elem;
  int p = S.n_rows;
  
  Rcpp::List out(L);
  Rcpp::List fit;
  arma::uvec par;
  arma::uvec npar(L);
  arma::vec loglik(L);
  
  for ( int l = 0; l < L; l++ ) {
    fit = gicf_core(start, adj, n, S, lambda(l), lambda_max, tolout, tolin, iterout, iterin);
    out[l] = fit;
    arma::mat sigma = fit["sigma"];
    start = sigma;
    par = find(sigma);
    npar(l) = (par.n_elem - p)/2.0 + p;
    loglik(l) = fit["loglik"];
  }
  
  return Rcpp::List::create(
    Named("out") = out,
    Named("loglik") = loglik,
    Named("npar") = npar
  );
}