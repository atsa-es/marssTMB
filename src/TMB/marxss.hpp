/// @file marxss.hpp
// MARSS model with covariates

#ifndef marxss_hpp
#define marxss_hpp

#include "marssTMB/isNA.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type marxss(objective_function<Type>* obj) {
  DATA_MATRIX(Y); /* n x T */
  DATA_MATRIX(d_Covar);
  DATA_MATRIX(c_Covar);
  DATA_INTEGER(no_c_covars);
  PARAMETER_MATRIX(X); /* State m x T */
  PARAMETER_MATRIX(x0);
  PARAMETER_MATRIX(V0); /* x[1] */
//  PARAMETER_MATRIX(Q); /* x[t] - x[t-1] */
  PARAMETER_MATRIX(C);
  PARAMETER_VECTOR(logsdQ); /* log of the sqrt of diag Q*/
  PARAMETER_VECTOR(cholCorrQ);
  PARAMETER_MATRIX(Z);
  PARAMETER_MATRIX(D);
  PARAMETER_VECTOR(logsdR); /* log of the sqrt of diag R*/
  PARAMETER_VECTOR(cholCorrR);
  
  int timeSteps = Y.row(0).size();
  int nY = Y.col(0).size(); /* n x T */
  int nX = X.col(0).size(); /* m x T */
  
  vector<Type> sdR=exp(logsdR); /* sd of R (diag) */
  vector<Type> sdQ=exp(logsdQ); /* sd of Q (diag) */
  
  // https://kaskr.github.io/adcomp/classdensity_1_1UNSTRUCTURED__CORR__t.html
  using namespace density;
  UNSTRUCTURED_CORR_t<Type> corMatGenR(cholCorrR);// This is the lower tri
  matrix<Type> FullCorrMatR = corMatGenR.cov(); /* full corr matrix has 1 on diag */
  UNSTRUCTURED_CORR_t<Type> corMatGenQ(cholCorrQ);// This is the lower tri
  matrix<Type> FullCorrMatQ = corMatGenQ.cov(); /* full corr matrix has 1 on diag */
  
  matrix<Type> predX(nX,1);  /* m x 1 */

  MVNORM_t<Type> initialState(V0);
  //MVNORM_t<Type> neg_log_density_process(Q);
  /* Define likelihood */
  Type ans=0;
  //ans -= dnorm(vector<Type>(u.row(0)),Type(0),Type(1),1).sum();
  ans += initialState(X.col(0)); /* tinitx=1 */
  for(int i=1;i<timeSteps;i++){ 
    //ans+= neg_log_density_process(u.row(i)-u.row(i-1)); // Process likelihood
    //vector<Type> differ = u.row(i)-u.row(i-1);
    // predX is m x 1 so u must be transposed
    // if statement is temporary until I can figure how create a
    // a diagonal matrix with 1 on the -1 diagonal
    // diag(1:(timeSteps+1))[1:timeSteps, 2:(timeSteps+1)]
    if(no_c_covars){
      predX = X.col(i-1) + C * c_Covar;
    }else{
      predX = X.col(i-1) + C * c_Covar.col(i);
    }
    vector<Type> differ = X.col(i)-predX;
    ans += VECSCALE(corMatGenQ,sdQ)(differ);
  }
  
  matrix<Type> predY(nY, timeSteps);  
  predY = Z * X + D * d_Covar;
  
  for(int i=0;i<timeSteps;i++){ //move one time step at a time
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
        subPred(j) = predY.col(i)(GoodVals(j));
        subSds(j) = sdR(GoodVals(j));
        for(int k=0;k<nonNAcount;k++){
          subCorr(j,k) = FullCorrMatR(GoodVals(j),GoodVals(k));
        }//end of loop through for truncated cormat
      }//end of removal of NA's from sets
      vector<Type> subDiffer = subData-subPred;
      ans += VECSCALE(MVNORM(subCorr),subSds)(subDiffer);
    }else{
      vector<Type> differ = Y.col(i)-predY.col(i);
      ans += VECSCALE(corMatGenR,sdR)(differ);
    }//end of data likelihood for this time step
  }//end of loop over time steps
  
  // Compute the full covariance matrices
  matrix<Type> FullCovMatR(nY,nY);
  matrix<Type> dSDR(nY,1);
  dSDR = sdR;
  FullCovMatR = dSDR.asDiagonal() * FullCorrMatR * dSDR.asDiagonal();
  matrix<Type> FullCovMatQ(nX,nX);
  matrix<Type> dSDQ(nX,1);
  dSDQ = sdQ;
  FullCovMatQ = dSDQ.asDiagonal() * FullCorrMatQ * dSDQ.asDiagonal();
  
  // Parameters with derivatives
  ADREPORT(X);
  ADREPORT(C);
  ADREPORT(FullCovMatQ);
  ADREPORT(Z);
  ADREPORT(D);
  ADREPORT(FullCovMatR);
  
  // Report wo derivatives
  REPORT(FullCorrMatQ);
  REPORT(FullCovMatQ);
  REPORT(FullCorrMatR);
  REPORT(FullCovMatR);
  
  return ans;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif