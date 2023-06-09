/// @file marss2.hpp
// MARSS model in vectorized marss form
// The Q and R are not split into diag and
// Allowed: covariates, time-varying and linear constraints

#ifndef marss2_hpp
#define marss2_hpp

#include "marssTMB/isNA.hpp"
#include "marssTMB/LOM.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type marss2(objective_function<Type>* obj) {
  DATA_MATRIX(Y); /* n x T */
  PARAMETER_MATRIX(X); /* State m x T */
  DATA_INTEGER(V0_is_zero);
  DATA_INTEGER(tinitx);
  DATA_IVECTOR(tfixed);
  DATA_IVECTOR(tfree);
  DATA_IVECTOR(numpar);
  DATA_STRUCT(free, LOM); // list of free matrices
  DATA_STRUCT(fixed, LOM); // list of fixed matrices
  DATA_STRUCT(par_dims, LOVi); // list of par matrices
  PARAMETER_VECTOR(pars); /* vector of parameters */
  
  int timeSteps = Y.row(0).size();
  int nY = Y.col(0).size(); /* n x T */
  int nX = X.col(0).size(); /* m x T */
  
  // Create the vec mat x T matrices
  int elem = 0;
  int TZ = std::max(tfree(elem), tfixed(elem));
  matrix<Type> Z(nX*nY, TZ);
  if(numpar(elem) != elem){
    Z = parvec(fixed(elem), free(elem), par(pars, numpar, elem), par_dims(elem), tfree(elem), tfixed(elem), TZ);
  }else{
    Z = fixed(elem);
  }
  elem = 1;
  int TA = std::max(tfree(elem), tfixed(elem));
  matrix<Type> A(nY, TA);
  if(numpar(elem) != elem){
    A = parvec(fixed(elem), free(elem), par(pars, numpar, elem), par_dims(elem), tfree(elem), tfixed(elem), TA);
  }else{
    A = fixed(elem);
  }
  elem = 2;
  int TR = std::max(tfree(elem), tfixed(elem));
  matrix<Type> R(nY*nY, TR);
  if(numpar(elem) != elem){
    R = parvec(fixed(elem), free(elem), par(pars, numpar, elem), par_dims(elem), tfree(elem), tfixed(elem), TR);
  }else{
    R = fixed(elem);
  }
  elem = 3;
  int TB = std::max(tfree(elem), tfixed(elem));
  matrix<Type> B(nX*nX, TB);
  if(numpar(elem) != elem){
    B = parvec(fixed(elem), free(elem), par(pars, numpar, elem), par_dims(elem), tfree(elem), tfixed(elem), TB);
  }else{
    B = fixed(elem);
  }
  elem = 4;
  int TU = std::max(tfree(elem), tfixed(elem));
  matrix<Type> U(nX, TU);
  if(numpar(elem) != elem){
    U = parvec(fixed(elem), free(elem), par(pars, numpar, elem), par_dims(elem), tfree(elem), tfixed(elem), TU);
  }else{
    U = fixed(elem);
  }
  elem = 5;
  int TQ = std::max(tfree(elem), tfixed(elem));
  matrix<Type> Q(nX*nX, TQ);
  if(numpar(elem) != elem){
    Q = parvec(fixed(elem), free(elem), par(pars, numpar, elem), par_dims(elem), tfree(elem), tfixed(elem), TQ);
  }else{
    Q = fixed(elem);
  }
  matrix<Type> x0(nX, 1);
  x0 = parmat2(fixed(6), free(6), par(pars, numpar, 6), par_dims(6));
  matrix<Type> V0(nX, nX);
  V0 = parmat2(fixed(7), free(7), par(pars, numpar, 7), par_dims(7));
  
  // Set the non time-varying parameters
  using namespace density;
  matrix<Type> FullCovMatQ = Q.col(0);
  FullCovMatQ.resize(nX,nX);
  FullCovMatQ = FullCovMatQ.transpose() * FullCovMatQ;
  MVNORM_t<Type> XVar(FullCovMatQ);
  
  matrix<Type> FullCovMatR = R.col(0);
  FullCovMatR.resize(nY,nY);
  FullCovMatR = FullCovMatR.transpose() * FullCovMatR;
  MVNORM_t<Type> YVar(FullCovMatR);

  // Compute the likelihoods
  matrix<Type> predX(nX,1);  /* m x 1 */

  Type ans=0; /* Define likelihood */
  if(V0_is_zero){
    if(tinitx){
      X.col(0) = x0.col(0);
    }else{
      predX.col(0) = x0.col(0) + U.col(0);
      vector<Type> differ = X.col(0)-predX.col(0);
      ans += XVar(differ); /* tinitx=1 */
    }
  }else{
    if(tinitx){
      MVNORM_t<Type> initialState(V0);
      ans += initialState(X.col(0)-x0.col(0)); /* tinitx=1 */
    }else{
      // Need to to check against the KF code
      predX.col(0) = x0.col(0) + U.col(0);
      vector<Type> differ0 = X.col(0)-predX.col(0);
      MVNORM_t<Type> initialState(V0+FullCovMatQ);
      ans += initialState(differ0); /* tinitx=1 */
    }
  }
int ui = 0;
for(int i=1;i<timeSteps;i++){ 
    // Process likelihood Note marss form so covariates appear in the
    // time-varying U
    if(TU != 1){
      ui = i;
    }
    predX = X.col(i-1) + U.col(ui);
    vector<Type> differ = X.col(i)-predX.col(0);
    ans += XVar(differ); /* tinitx=1 */
}
  
  matrix<Type> predY(nY,1);  
  predY.setZero();
  matrix<Type> ZZ(nY*nX, 1);
  
  //  std::cout << std::scientific;
  //  std::cout << Xi.rows() << std::endl;
  //  std::cout << Xi.cols() << std::endl;
  ZZ = Z.col(0);
  ZZ.resize(nY, nX);
  int ai = 0;
  
  for(int i=0;i<timeSteps;i++){ //move one time step at a time
if(TA != 1){ 
  ai = i;
}
if(TZ != 1){
  ZZ = Z.col(i);
  ZZ.resize(nY, nX);
}
predY = ZZ * X.col(i) + A.col(ai);

    int nonNAcount = 0; //start at zero NA values
    vector<int> GoodVals(nY);
    for(int j=0;j<nY;j++){//loop over all time series for this time step
      if(!isNA(Y.col(i)(j))){//if value is not NA
        GoodVals(nonNAcount) = j; //add position to good values (good values only stored in beginning of vector)
        nonNAcount++; //increment the values of
      }
    }
    if(nonNAcount == 0) continue; // no data
    if(nonNAcount<nY){ 
      //if NA values present and some data
      matrix<Type> subCov(nonNAcount,nonNAcount);
      vector<Type> subData(nonNAcount);
      vector<Type> subPred(nonNAcount);
      
      for(int j=0;j<nonNAcount;j++){
        subData(j) = Y.col(i)(GoodVals(j));
        subPred(j) = predY.col(0)(GoodVals(j));
        for(int k=0;k<nonNAcount;k++){
          subCov(j,k) = FullCovMatR(GoodVals(j),GoodVals(k));
        }//end of loop through for truncated covmat
      }//end of removal of NA's from sets
      vector<Type> subDiffer = subData-subPred;
      ans += MVNORM(subCov)(subDiffer); /* tinitx=1 */
    }else{
      vector<Type> differ = Y.col(i)-predY.col(0);
      ans += YVar(differ); /* tinitx=1 */
    }//end of data likelihood for this time step
  }//end of loop over time steps
  
  // Parameters with derivatives
  ADREPORT(X);
  ADREPORT(pars);

  // Report wo derivatives
  REPORT(X);
  REPORT(FullCovMatQ);
  REPORT(FullCovMatR);
  
  return ans;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif