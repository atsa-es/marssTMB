/// @file test.hpp
// MARSS model with covariates

#ifndef test_hpp
#define test_hpp

#include "marssTMB/isNA.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type test(objective_function<Type>* obj) {
  DATA_MATRIX(Y); /* n x T */
  PARAMETER_MATRIX(X); /* State m x T */
  PARAMETER_VECTOR(Rdiag);
  PARAMETER_VECTOR(Roffdiag);
  PARAMETER_VECTOR(Qdiag);
  PARAMETER_VECTOR(Qoffdiag);

  int timeSteps = Y.row(0).size();
  int nY = Y.col(0).size(); /* n x T */
  int nX = X.col(0).size(); /* m x T */
  
  // Create the vec mat x T matrices
  matrix<Type> x0(nX, 1);
  x0.setZero();
  matrix<Type> V0(nX, nX);
  V0.setIdentity();

  vector<Type> sdQ = exp(Qdiag);
  vector<Type> sdR = exp(Rdiag);
  
  // https://kaskr.github.io/adcomp/classdensity_1_1UNSTRUCTURED__CORR__t.html
  using namespace density;
  UNSTRUCTURED_CORR_t<Type> corMatGenR(Roffdiag);// This is the lower tri
  matrix<Type> FullCorrMatR = corMatGenR.cov(); /* full corr matrix has 1 on diag */
  UNSTRUCTURED_CORR_t<Type> corMatGenQ(Qoffdiag); // This is the lower tri
  matrix<Type> FullCorrMatQ = corMatGenQ.cov(); /* full corr matrix has 1 on diag */
  
  // Compute the likelihoods
  Type ans=0; /* Define likelihood */
  
  for(int i=1;i<timeSteps;i++){ 
    vector<Type> differ = X.col(i)-X.col(i-1);
    ans += VECSCALE(corMatGenQ,sdQ)(differ);
  }
  
  matrix<Type> predY(nY,1);  

  for(int i=0;i<timeSteps;i++){ //move one time step at a time
    predY = X.col(i);

    int nonNAcount = 0; //start at zero NA values
    vector<int> GoodVals(nY);
    for(int j=0;j<nY;j++){//loop over all time series for this time step
      if(!isNA(Y.col(i)(j))){//if value is not NA
        GoodVals(nonNAcount) = j; //add position to good values (good values only stored in beginning of vector)
        nonNAcount++; //increment the values of
      }
    }
    if(nonNAcount<nY){ //if NA values present
      matrix<Type> subCorr(nonNAcount,nonNAcount);
      vector<Type> subSds(nonNAcount);
      vector<Type> subData(nonNAcount);
      vector<Type> subPred(nonNAcount);
      
      for(int j=0;j<nonNAcount;j++){
        subData(j) = Y.col(i)(GoodVals(j));
        subPred(j) = X.col(i)(GoodVals(j));
        subSds(j) = sdR(GoodVals(j));
        for(int k=0;k<nonNAcount;k++){
          subCorr(j,k) = FullCorrMatR(GoodVals(j),GoodVals(k));
        }//end of loop through for truncated cormat
      }//end of removal of NA's from sets
      vector<Type> subDiffer = subData-subPred;
      ans += VECSCALE(MVNORM(subCorr),subSds)(subDiffer);
    }else{
      vector<Type> differ = Y.col(i)-X.col(i);
      ans += VECSCALE(corMatGenR, sdR)(differ);
    }//end of data likelihood for this time step
  }//end of loop over time steps
  
  // Parameters with derivatives
  ADREPORT(X);
  ADREPORT(Rdiag);
  ADREPORT(Roffdiag);
  ADREPORT(Qdiag);
  ADREPORT(Qoffdiag);
  
  return ans;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif