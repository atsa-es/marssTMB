/// @file marxss.hpp
// MARSS model with covariates

#ifndef marxss_hpp
#define marxss_hpp

#include "marssTMB/isNA.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type marxss(objective_function<Type>* obj) {
  DATA_MATRIX(obs); /*  timeSteps x stateDim*/
  DATA_MATRIX(Covar);
  DATA_INTEGER(has_covars);
  PARAMETER_VECTOR(logsdR); /* log of the diag of chol R*/
  PARAMETER_VECTOR(cholCorrR);
  PARAMETER_MATRIX(Q); /* x[t] - x[t-1] */
  PARAMETER_MATRIX(V0); /* x[1] */
  PARAMETER_MATRIX(D);
  PARAMETER_MATRIX(Z);
  PARAMETER_MATRIX(x0);
  PARAMETER_MATRIX(u); /* State */
  
  int timeSteps=obs.col(0).size();
  int obsDim=obs.row(0).size();
  
  vector<Type> sdR=exp(logsdR); /* sd of R (diag) */
  
  using namespace density;
  UNSTRUCTURED_CORR_t<Type> corMatGen(cholCorrR);// This is the full Cormat
  matrix<Type> FullCorrMatR=corMatGen.cov(); /* has 1 on diag */
  
  MVNORM_t<Type> initialState(V0);
  MVNORM_t<Type> neg_log_density_process(Q);
  /* Define likelihood */
  Type ans=0;
  //ans -= dnorm(vector<Type>(u.row(0)),Type(0),Type(1),1).sum();
  ans += initialState(u.row(0));
  for(int i=1;i<timeSteps;i++){ 
    ans+= neg_log_density_process(u.row(i)-u.row(i-1)); // Process likelihood 
  }
  
  matrix<Type> pred(timeSteps,obsDim);  
  pred = Z * u.transpose();
  if(has_covars) pred += D * Covar;
  
  for(int i=0;i<timeSteps;i++){ //move one time step at a time
    int nonNAcount = 0; //start at zero NA values
    vector<int> GoodVals(obs.row(i).size());
    for(int j=0;j<obs.row(i).size();j++){//loop over all time series for this time step
      if(!isNA(obs.row(i)(j))){//if value is not NA
        GoodVals(nonNAcount) = j; //add position to good values (good values only stored in beginning of vector)
        nonNAcount++; //increment the values of
      }
    }
    if(nonNAcount<obs.row(i).size()){ //if NA values present
      matrix<Type> subCorr(nonNAcount,nonNAcount);
      vector<Type> subSds(nonNAcount);
      vector<Type> subData(nonNAcount);
      vector<Type> subPred(nonNAcount);
      
      for(int j=0;j<nonNAcount;j++){
        subData(j) = obs.row(i)(GoodVals(j));
        subPred(j) = pred.transpose().row(i)(GoodVals(j));
        subSds(j) = sdR(GoodVals(j));
        for(int k=0;k<nonNAcount;k++){
          subCorr(j,k) = FullCorrMatR(GoodVals(j),GoodVals(k));
        }//end of loop through for truncated cormat
      }//end of removal of NA's from sets
      vector<Type> subDiffer = subData-subPred;
      ans += VECSCALE(MVNORM(subCorr),subSds)(subDiffer);
    }else{
      vector<Type> differ = obs.row(i)-pred.transpose().row(i);
      ans += VECSCALE(corMatGen,sdR)(differ);
    }//end of data likelihood for this time step
  }//end of loop over time steps
  
  matrix<Type> FullCovMatR(obsDim,obsDim);
  matrix<Type> dSD(obsDim,1);
  dSD = sdR;
  FullCovMatR = dSD.asDiagonal() * FullCorrMatR * dSD.asDiagonal();
  ADREPORT(Z);
  ADREPORT(D);
  ADREPORT(u);
  ADREPORT(FullCovMatR); /* with derivative */
  REPORT(FullCorrMatR); /* wo derivative */
  REPORT(FullCovMatR); /* wo derivative */
  
  return ans;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif