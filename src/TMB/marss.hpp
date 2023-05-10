/// @file marss.hpp
// MARSS model in vectorized marss form
// Allowed: covariates, time-varying and linear constraints

#ifndef marss_hpp
#define marss_hpp

#include "marssTMB/isNA.hpp"
#include "marssTMB/LOM.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type marss(objective_function<Type>* obj) {
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
  matrix<Type> x0(nX, 1);
  x0 = parmat2(fixed(4), free(4), par(pars, numpar, 4), par_dims(4));
  matrix<Type> V0(nX, nX);
  V0 = parmat2(fixed(5), free(5), par(pars, numpar, 5), par_dims(5));
  matrix<Type> U(nX, timeSteps);
  U = parvec(fixed(3), free(3), par(pars, numpar, 3), par_dims(3), tfree(3), tfixed(3), timeSteps);
  matrix<Type> A(nY, timeSteps);
  A = parvec(fixed(1), free(1), par(pars, numpar, 1), par_dims(1), tfree(1), tfixed(1), timeSteps);
  matrix<Type> Z(nX*nY, timeSteps);
  Z = parvec(fixed(0), free(0), par(pars, numpar, 0), par_dims(0), tfree(0), tfixed(0), timeSteps);
  matrix<Type> B(nX*nX, timeSteps);
  B = parvec(fixed(2), free(2), par(pars, numpar, 2), par_dims(2), tfree(2), tfixed(2), timeSteps);

  // Set the non time-varying parameters
  matrix<Type> QdiagMat(nX, timeSteps);
  matrix<Type> QoffdiagMat(nX * (nX - 1) / 2, timeSteps);
  QdiagMat = parvec(fixed(9), free(9), par(pars, numpar, 9), par_dims(9), tfree(9), tfixed(9), timeSteps);
  vector<Type> sdQ=QdiagMat.col(0); /* log sd of Q (diag) */
  sdQ = exp(sdQ);
  vector<Type> cholCorrQ(nX * (nX - 1) / 2);
  if(nX>1){
    QoffdiagMat = parvec(fixed(10), free(10), par(pars, numpar, 10), par_dims(10), tfree(10), tfixed(10), timeSteps);
    cholCorrQ = QoffdiagMat.col(0);
  }
  // Set the non time-varying parameters
  matrix<Type> RdiagMat(nY, timeSteps);
  matrix<Type> RoffdiagMat(nY * (nY - 1) / 2, timeSteps);
  RdiagMat = parvec(fixed(11), free(11), par(pars, numpar, 11), par_dims(11), tfree(11), tfixed(11), timeSteps);
  vector<Type> sdR=RdiagMat.col(0); /* log sd of R (diag) */
  sdR = exp(sdR);
  vector<Type> cholCorrR(nY * (nY - 1) / 2);
  if(nY>1){
    RoffdiagMat = parvec(fixed(12), free(12), par(pars, numpar, 12), par_dims(12), tfree(12), tfixed(12), timeSteps);
    cholCorrR = RoffdiagMat.col(0);
  }

  // https://kaskr.github.io/adcomp/classdensity_1_1UNSTRUCTURED__CORR__t.html
  using namespace density;
  UNSTRUCTURED_CORR_t<Type> corMatGenR(cholCorrR);// This is the lower tri
  matrix<Type> FullCorrMatR = corMatGenR.cov(); /* full corr matrix has 1 on diag */
  UNSTRUCTURED_CORR_t<Type> corMatGenQ(cholCorrQ);// This is the lower tri
  matrix<Type> FullCorrMatQ = corMatGenQ.cov(); /* full corr matrix has 1 on diag */
  
  // Compute the likelihoods
  matrix<Type> predX(nX,1);  /* m x 1 */

  Type ans=0; /* Define likelihood */
  if(V0_is_zero){
    if(tinitx){
      X.col(0) = x0.col(0);
    }else{
      predX = x0.col(0) + U.col(0);
      vector<Type> differ0 = X.col(0)-predX.col(0);
      ans += VECSCALE(corMatGenQ,sdQ)(differ0);
    }
  }else{
    // I think this is wrong if tinitx=0
    MVNORM_t<Type> initialState(V0);
    ans += initialState(X.col(0)-x0.col(0)); /* tinitx=1 */
  }

  for(int i=1;i<timeSteps;i++){ 
    // Process likelihood Note marss form so covariates appear in the
    // time-varying U
    predX = X.col(i-1) + U.col(i);
    vector<Type> differ = X.col(i)-predX.col(0);
    ans += VECSCALE(corMatGenQ,sdQ)(differ);
  }
  
  matrix<Type> predY(nY,1);  
//  predY.setZero();
  matrix<Type> I(nY,nY);
  I.setIdentity();
  matrix<Type> Xi(1,nX);

//  std::cout << std::scientific;
//  std::cout << Xi.rows() << std::endl;
//  std::cout << Xi.cols() << std::endl;
  
  for(int i=0;i<timeSteps;i++){ //move one time step at a time
    Xi = X.col(i).transpose();
    predY = tmbutils::kronecker(I, Xi) * Z.col(i) + A.col(i);

    int nonNAcount = 0; //start at zero NA values
    vector<int> GoodVals(nY);
    for(int j=0;j<nY;j++){//loop over all time series for this time step
      if(!isNA(Y.col(i)(j))){//if value is not NA
        GoodVals(nonNAcount) = j; //add position to good values (good values only stored in beginning of vector)
        nonNAcount++; //increment the values of
      }
    }
    if(nonNAcount == 0) continue; // no data
    if(nonNAcount<nY){ //if NA values present
      matrix<Type> subCorr(nonNAcount,nonNAcount);
      vector<Type> subSds(nonNAcount);
      vector<Type> subData(nonNAcount);
      vector<Type> subPred(nonNAcount);
      
      for(int j=0;j<nonNAcount;j++){
        subData(j) = Y.col(i)(GoodVals(j));
        subPred(j) = predY.col(0)(GoodVals(j));
        subSds(j) = sdR(GoodVals(j));
        for(int k=0;k<nonNAcount;k++){
          subCorr(j,k) = FullCorrMatR(GoodVals(j),GoodVals(k));
        }//end of loop through for truncated cormat
      }//end of removal of NA's from sets
      vector<Type> subDiffer = subData-subPred;
      ans += VECSCALE(MVNORM(subCorr),subSds)(subDiffer);
    }else{
      vector<Type> differ = Y.col(i)-predY.col(0);
      ans += VECSCALE(corMatGenR, sdR)(differ);
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