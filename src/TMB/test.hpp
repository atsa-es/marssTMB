/// @file test.hpp
// univariate time series

#ifndef test_hpp
#define test_hpp

#include "marssTMB/LOM.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type test(objective_function<Type>* obj) {

  using namespace density;
  
  DATA_VECTOR(Y);
  DATA_INTEGER(n);
  DATA_INTEGER(est_drift);
  DATA_INTEGER(est_rho);
  DATA_IVECTOR(keep);
  DATA_STRUCT(free, LOM); // list of free matrices
  DATA_STRUCT(fixed, LOM); // list of fixed matrices
  DATA_STRUCT(par, LOM); // list of par matrices
  DATA_STRUCT(par_dims, LOVi); // list of par matrices

  PARAMETER(u); // drift
  PARAMETER(logit_rho); // inv-logit(b)
  PARAMETER(log_obs_sd); // obs error
  PARAMETER(log_pro_sd); // pro error
  PARAMETER_VECTOR(x); // n x 1
  //PARAMETER(x0);
  
  Type obs_sigma = exp(log_obs_sd);
  Type pro_sigma = exp(log_pro_sd);
  Type rho;
  if(est_rho==1) rho = exp(logit_rho)/(1+exp(logit_rho));
  
  // random effects / penalties
  Type nll=0; // total function to maximize
  vector<Type> pred(n); // predictions
  
  // process deviations
  for(int i = 0; i < n; i++) {
    nll -= dnorm(x(i), Type(0.0), pro_sigma, true);
  }
  // initial conditions
  pred(0) = x(0);
  // process / random walk
  for(int i = 1; i < n; i++) {
    pred(i) = pred(i-1);
    if(est_rho == 1) pred(i) = rho * pred(i);
    if(est_drift == 1) pred(i) = pred(i) + u;
    pred(i) = pred(i) + x(i); // process stochasticity
  }
  
  // observation likelihood
  for(int i = 0; i < n; i++) {
    if(keep(i) == 1) nll -= dnorm(Y(i), pred(i), obs_sigma, true); 
  }
  
  if(est_rho) {
    ADREPORT(rho); 
    REPORT(rho); 
  }
  if(est_drift) {
    ADREPORT(u); 
    REPORT(u); 
  }
  ADREPORT(obs_sigma); 
  REPORT(obs_sigma);
  ADREPORT(pro_sigma); 
  REPORT(pro_sigma);
  ADREPORT(pred); 
  REPORT(pred);
  int nnn =par_dims(0)(0);
  int ppp =par_dims(0)(1);
  matrix<Type> DD(nnn, ppp);
  DD = parmat(fixed(0), free(0), par(0), par_dims(0));
  REPORT(DD);
  // end
  return (nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif