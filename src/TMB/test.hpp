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
  DATA_MATRIX(mY);
  DATA_MATRIX(X);
  DATA_IVECTOR(tfixed);
  DATA_IVECTOR(tfree);
  DATA_IVECTOR(numpar);
  DATA_STRUCT(free, LOM); // list of free matrices
  DATA_STRUCT(fixed, LOM); // list of fixed matrices
  DATA_STRUCT(par_dims, LOVi); // list of par matrices
  DATA_VECTOR(pars); /* vector of parameters */
  
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

  int timeSteps = mY.row(0).size();
  int nY = mY.col(0).size(); /* n x T */
  int nX = X.col(0).size(); /* m x T */
  
  
  matrix<Type> x0(nX, 1);
  x0 = parmat(fixed(4), free(4), par(pars, numpar, 4), par_dims(4));
  matrix<Type> V0(nX, nX);
  V0 = parmat(fixed(5), free(5), par(pars, numpar, 5), par_dims(5));
  matrix<Type> U(nX, timeSteps);
  U = parvec(fixed(3), free(3), par(pars, numpar, 3), par_dims(3), tfree(3), tfixed(3), timeSteps);
  matrix<Type> A(nY, timeSteps);
  A = parvec(fixed(1), free(1), par(pars, numpar, 1), par_dims(1), tfree(1), tfixed(1), timeSteps);
  matrix<Type> Z(nX*nY, timeSteps);
  Z = parvec(fixed(0), free(0), par(pars, numpar, 0), par_dims(0), tfree(0), tfixed(0), timeSteps);
  matrix<Type> B(nX*nX, timeSteps);
  B = parvec(fixed(2), free(2), par(pars, numpar, 2), par_dims(2), tfree(2), tfixed(2), timeSteps);

  // Set the non time-varying parameters
  matrix<Type> Qdiag(nX, timeSteps);
  matrix<Type> Qoffdiag(nX * (nX - 1) / 2, timeSteps);
  Qdiag = parvec(fixed(9), free(9), par(pars, numpar, 9), par_dims(9), tfree(9), tfixed(9), timeSteps);
  vector<Type> sdQ=Qdiag.col(0); /* sd of Q (diag) */
  sdQ = exp(sdQ);
  vector<Type> cholCorrQ(nX * (nX - 1) / 2);
  if(nX>1){
    Qoffdiag = parvec(fixed(10), free(10), par(pars, numpar, 10), par_dims(10), tfree(10), tfixed(10), timeSteps);
    cholCorrQ = Qoffdiag.col(0);
  }
  // Set the non time-varying parameters
  matrix<Type> Rdiag(nY, timeSteps);
  matrix<Type> Roffdiag(nY * (nY - 1) / 2, timeSteps);
  Rdiag = parvec(fixed(11), free(11), par(pars, numpar, 11), par_dims(11), tfree(11), tfixed(11), timeSteps);
  vector<Type> sdR=Rdiag.col(0); /* sd of Q (diag) */
  sdR = exp(sdR);
  vector<Type> cholCorrR(nY * (nY - 1) / 2);
  if(nY>1){
    Roffdiag = parvec(fixed(12), free(12), par(pars, numpar, 12), par_dims(12), tfree(12), tfixed(12), timeSteps);
    cholCorrR = Roffdiag.col(0);
  }
  
  REPORT(sdQ);
  REPORT(cholCorrQ);
  REPORT(Rdiag);
  REPORT(cholCorrR);
  REPORT(Z);
  REPORT(A);
  REPORT(B);
  // end;
  
  return (nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif