# Tom Cline (USGS) wrote this while a grad student at SAFS
# EEH. I editted to include x0 distribution as in MARSS.dfa
library(TMB)
library(Matrix)

dfa_model <- "// Dynamic Factor Analysis for multivariate time series
#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(obs); /*  timeSteps x stateDim*/
  DATA_MATRIX(Covar);
  PARAMETER_VECTOR(logsdObs);
  PARAMETER_VECTOR(cholCorr);
  PARAMETER_MATRIX(covState); /* x[t] - x[t-1] */
  PARAMETER_MATRIX(covinitState); /* x[1] */
  PARAMETER_MATRIX(D);
  PARAMETER_MATRIX(Z);
  PARAMETER_MATRIX(u); /* State */
  
  
  int timeSteps=obs.col(0).size();
  int obsDim=obs.row(0).size();
  
  vector<Type> sdObs=exp(logsdObs);
  
  using namespace density;
  UNSTRUCTURED_CORR_t<Type> corMatGen(cholCorr);// This is the full Cormat
  matrix<Type> FullCorrMat=corMatGen.cov();
  
 
  MVNORM_t<Type> initialState(covinitState);
  MVNORM_t<Type> neg_log_density_process(covState);
  /* Define likelihood */
  Type ans=0;
  //ans -= dnorm(vector<Type>(u.row(0)),Type(0),Type(1),1).sum();
  ans += initialState(u.row(0));
  for(int i=1;i<timeSteps;i++){ 
    ans+= neg_log_density_process(u.row(i)-u.row(i-1)); // Process likelihood 
  }
  
  matrix<Type> pred(timeSteps,obsDim);  
  pred = (Z * u.transpose()) + (D * Covar);
  
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
			subSds(j) = sdObs(GoodVals(j));
			for(int k=0;k<nonNAcount;k++){
				subCorr(j,k) = FullCorrMat(GoodVals(j),GoodVals(k));
			}//end of loop through for truncated cormat
		}//end of removal of NA's from sets
		vector<Type> subDiffer = subData-subPred;
		ans += VECSCALE(MVNORM(subCorr),subSds)(subDiffer);
	 }else{
	   	vector<Type> differ = obs.row(i)-pred.transpose().row(i);
	 	ans += VECSCALE(corMatGen,sdObs)(differ);
	 }//end of data likelihood for this time step
  }//end of loop over time steps
  
  matrix<Type> FullCovMat(obsDim,obsDim);
  matrix<Type> dSD(obsDim,1);
  dSD = sdObs;
  FullCovMat = dSD.asDiagonal() * FullCorrMat * dSD.asDiagonal();
  ADREPORT(Z);
  ADREPORT(D);
  ADREPORT(u);
  ADREPORT(FullCovMat);
  
  REPORT(FullCorrMat);
  
  return ans;
}"
  
  #A nasty bit of code that basically checks if the the model above is EXACTLY the same as the already compiled version. Saves a small amount of time for rerunning code.
  if(file.exists('dfa1tmb.cpp')){
    t1<-read.table('dfa1tmb.cpp',sep="\t",stringsAsFactors=F,comment.char='', quote = "")
    write(dfa_model, file = "dfa1tmb_temp.cpp")
    t2<-read.table('dfa1tmb_temp.cpp',sep="\t",stringsAsFactors=F,comment.char='', quote = "")
    file.remove('dfa1tmb_temp.cpp')
    #t1<-data.frame(t1[-which(t1==''),])
    #t2<-data.frame(strsplit(dfa_model,fixed=T,'\n'))
    #t2<-data.frame(t2[-which(t2==''),])
    if(identical(t1,t2)){
      dyn.load(dynlib("dfa1tmb"))
    }else{
      write(dfa_model, file = "dfa1tmb.cpp")
      compile("dfa1tmb.cpp")
      dyn.load(dynlib("dfa1tmb"))
    }
  }else{
    write(dfa_model, file = "dfa1tmb.cpp")
    compile("dfa1tmb.cpp")
    dyn.load(dynlib("dfa1tmb"))
  }
  
  # if(!exists('dfamodel')){
  # 	write(dfa_model, file = "dfa1tmb.cpp")
  # 	dfamodel<-dfa_model
  # 	compile("dfa1tmb.cpp")
  # 	dyn.load(dynlib("dfa1tmb"))
  # }else if(dfa_model != dfamodel){
  # 	write(dfa_model, file = "dfa1tmb.cpp")
  # 	dfamodel<-dfa_model
  # 	compile("dfa1tmb.cpp")
  # 	dyn.load(dynlib("dfa1tmb"))
  # }
  
  #This function generates the Z matrix based on the number of times series and the number of states that are estimated. It is called from with the run_dfa function.
  ZmatGen<-function(Data,NumStates){
    tZ<-matrix(0.5,nrow=ncol(Data),ncol=NumStates)
    if(NumStates>1){
      for(i in 1:(NumStates-1)){
        tZ[i,(NumStates-(NumStates-2)+(i-1)):NumStates]<-rep(0,NumStates-1-(i-1))
      }
      return(tZ)
    }else{
      return(tZ)
    }
  }
  
  #This function creates a Z mat factor which allows TMB to fix certain parameters. This is required as the upper triangle of the Zmatrix must fixed at 0 to allow the model to be identifiable. It is called from within the model run.
  ZmatFactorGen<-function(Data,NumStates){
    tZ<-matrix(seq(1,ncol(Data)*NumStates),nrow=ncol(Data),ncol=NumStates)
    if(NumStates>1){
      for(i in 1:(NumStates-1)){
        tZ[i,(NumStates-(NumStates-2)+(i-1)):NumStates]<-rep(NA,NumStates-1-(i-1))
      }
      tZ[!is.na(tZ)]<-seq(1,sum(!is.na(tZ)))
      return(as.factor(tZ))
    }else{
      return(as.factor(tZ))
    }
  }
  
  ### Compute AIC for my DFA model output
  dfaAIC<-function(x,AICc=F){
    opt<-x[['Optimization']]
    NumPar<-length(opt$par)
    NLL<-opt$value#opt$value
    AIC<-2*NumPar + 2*NLL
    if(AICc){
      AIC<-AIC + (2*NumPar*(NumPar+1))/(nrow(x[['Estimates']]$Z)*ncol(x[['Estimates']]$u)-NumPar-1)
    }
    return(AIC)
  }

  dfaLL<-function(x){
    opt<-x[['Optimization']]
    NumPar<-length(opt$par)
    NLL<-opt$value#opt$value
    return(-1*NLL)
  }
  
  ### This function runs the DFA using the TMB model above
  ### This function does all of the prep work for matrices based on input about the number of states, error structure, and covariates. 
  ### !!!!!!!!! obs must be a matrix with time going across the columns (same as MARSS package !!!!!!!!).
  
  #All time series and covariates need to be Zscored prior to running the model.
  
  # When you call this function enter the number of states (default is 1). this is parameter 'm' in the MARSS package.
  # Enter the error structure. Currently the only choices are Diagonal and Equal ('DE' also the current default), Diagonal and Unequal ('DUE'), or Unconstrained 'UNC'
  # Covariates:
  # First if you want to include covariates 'EstCovar' needs to be TRUE
  # Include the covariate time series as Covars. NO MISSING VALUES.
  # The Default covariate configuration is individual paramter estimates for each time series but from a single covariate time series. Example how does PDO (one covariate time series) affect the 9 sockeye runs of bristol bay (each time series gets a PDO effect)
  # If you want individual paramters from individual covariate series then you need to set indivCovar to TRUE. Example you have temperature covariates for each of the rivers of bristol bay and want to estimated individual temperature effects.
  # A combination of these two will require manual passing of the covariate paramter matrix. You can do this by supplying Dmat and Dfac. This defaults to NULL and don't mess with it unless you need to deviate from the two covariate approaches above. BEWARE!: frustrating debugging is an absolute certainty should you go this route. But if its what you need it can be done.
  
  runDFA<-function(obs,NumStates=1,ErrStruc='DE',EstCovar=FALSE,Covars=NULL,indivCovar=FALSE,Dmat=NULL,Dfac=NULL,Rfac=NULL,logsdObs=NULL,logsdObsFac=NULL,cholCorr=NULL,cholFac=NULL,EstSE=FALSE){ 
    ##
    #  TopSection is used for Debug only.
    #
    #NumStates=1;ErrStruc='UNC';EstCovar=TRUE;Covars=Z_run;indivCovar=FALSE;Dmat=NULL;Dfac=NULL;Rfac=NULL;logsdObs=NULL;logsdObsFac=NULL;cholCorr=NULL;cholFac=NULL
    #obs<-Z_L_FW1
    #Dmat<-DmatRun
    #Dfac<-DfacRun
    
    obs<-t(obs)
    Zfac<-ZmatFactorGen(Data=obs,NumStates=NumStates) #creates the Z factor to fix the upper corner if NumStates is greater than 1.
    if(EstCovar){ # If you are estimating covariates this creates Dmat and Dfac based on the data, number of covars
      if(is.null(Dmat) & is.null(Dfac)){ #This checks that you did not supply Dmat or Dfac (Manual covariate paramter entries)
        if(!indivCovar){
          Dmat<-matrix(rep(0,ncol(obs)*nrow(Covars)),ncol=nrow(Covars),nrow=ncol(obs))
          Dfac<-as.factor(seq(1,ncol(obs)*nrow(Covars)))
        }else{
          Dmat<-matrix(0,ncol=nrow(Covars),nrow=ncol(obs))
          diag(Dmat)<-rep(0,nrow(Covars))#rnorm(nrow(Covars),0,1)
          Dfac<-matrix(NA,ncol=nrow(Covars),nrow=ncol(obs))
          diag(Dfac)<-seq(1,nrow(Covars))
          Dfac<-as.factor(Dfac)
        }
      }
      data <- list(obs=obs,NumState=NumStates,Covar=Covars)
    }else{ #If you are not estimating covariates we just pass a time series of 0's, a Dmat of 0's, and and NA factors so the paramters will not be estimated
      
      Dmat<-matrix(0,ncol=1,nrow=ncol(obs))
      Dfac<-as.factor(rep(NA,ncol(obs)))
      Covars<-matrix(0,nrow=1,ncol=nrow(obs))
      data <- list(obs=obs,NumState=NumStates,Covar=Covars)
    }
    
    #This set of if-else statements creates the proper parameter set for the error structure selected
    if(is.null(logsdObs) & is.null(logsdObsFac) & is.null(cholCorr) & is.null(cholFac)){
      if(ErrStruc=='DE'){
        cholCorr<-rep(0,ncol(obs)*(ncol(obs)-1)/2)
        logsdObs<-log(rep(0.5,ncol(obs)))
        logsdObsFac<-rep(1,ncol(obs))
        logsdObsFac<-factor(logsdObsFac)
        cholFac <-rep(NA,ncol(obs)*(ncol(obs)-1)/2)
        cholFac<-factor(cholFac)
      }else if(ErrStruc =='DUE'){
        cholCorr<-rep(0,ncol(obs)*(ncol(obs)-1)/2)
        logsdObs<-log(rep(0.5,ncol(obs)))
        logsdObsFac<-seq(1,ncol(obs))
        logsdObsFac<-factor(logsdObsFac)
        cholFac <-rep(NA,ncol(obs)*(ncol(obs)-1)/2)
        cholFac<-factor(cholFac)
      }else if(ErrStruc == 'UNC'){
        cholCorr<-rep(0,ncol(obs)*(ncol(obs)-1)/2)
        logsdObs<-log(rep(0.5,ncol(obs)))
        logsdObsFac<-seq(1,ncol(obs))
        logsdObsFac<-factor(logsdObsFac)
        cholFac <-seq(1,(ncol(obs)*(ncol(obs)-1)/2))
        cholFac<-factor(cholFac)
      }}
    
    #Creates the input parameter list
    parameters <- list(
      logsdObs = logsdObs,
      cholCorr = cholCorr,
      covState = diag(1,NumStates),
      covinitState = diag(5,NumStates),
      D=Dmat,
      Z=ZmatGen(Data=obs,NumStates=NumStates),
      u=matrix(0,nrow=nrow(obs),ncol=NumStates)
     )
    covinitStateFac<-factor(matrix(NA,nrow=NumStates,ncol=NumStates))
    covStateFac<-factor(matrix(NA,nrow=NumStates,ncol=NumStates))
    #Creates the model object and runs the optimization
    obj1 <- MakeADFun(data, parameters, random="u", DLL="dfa1tmb", silent=T, map=list(Z=Zfac, D=Dfac, cholCorr=cholFac, logsdObs=logsdObsFac, covState=covStateFac, covinitState=covinitStateFac))
    opt1 <- nlminb(obj1$par,obj1$fn,obj1$gr,control=list(iter.max=2000,eval.max=2000))
    #newtonOption(obj1,smartsearch=TRUE)
    obj1$control=list(trace=1,REPORT=1,reltol=1e-12,maxit=2000)
    obj1$fn()
    obj1$gr()
    #obj1$method='BFGS'
    obj1$par=opt1$par
    #obj1$lower
    system.time(opt1 <- do.call("optim",obj1))
    pl1 <- obj1$env$parList()#This contains all of your parameter estimates RAW as they come out of the optimizer
    if(EstSE){
      sdr<-sdreport(obj1)
    }
    
    # #Trend directions are arbitrary adjust them so that most load positively
    # for(i in 1:NumStates){
    #   if(median(pl1$Z[,i])<0){
    #     pl1$Z[,i]<--pl1$Z[,i]
    #     pl1$u[,i]<--pl1$u[,i]
    #   }
    # }
    # ScaleFac<-as.vector(apply(pl1$u,2,FUN=sd))
    # pl1$u<-t(t(pl1$u)/ScaleFac)
    # pl1$Z<-t(t(pl1$Z)*ScaleFac)
    
    # # Do the Varimax rotation from models with more that one trend so I dont have to do it later.
    # if(NumStates>1){
    #   H.inv = varimax(pl1$Z)$rotmat
    #   Z.rot = pl1$Z %*% H.inv #maximum variance explained
    #   trends.rot = solve(H.inv) %*% t(pl1$u)
    #   
    #   pl1$Z<-Z.rot
    #   pl1$u<-t(trends.rot)
    # }
    pl1$u<-t(pl1$u)
    
    if(EstSE){pl1$R <- matrix(sdr$value[which(names(sdr$value)=='FullCovMat')],nrow=length(logsdObs),ncol=length(logsdObs))
    }else{pl1$R<- diag(exp(pl1$logsdObs)) %*% obj1$report()$FullCorrMat  %*% diag(exp(pl1$logsdObs))}
    
    #Fits for each time series
    FitSeries<- pl1$Z %*% pl1$u + pl1$D %*% Covars
    
    #Create plots of model fit
    #pdf('DFAfit.pdf')
    #for(i in 1:ncol(obs)){
    #plot(obs[,i])
    #points(FitSeries[i,],type='l')
    #}
    #dev.off()
    
    #Standard Errors for parameters
    if(EstSE){
      SES <- list(D=sdr$sd[which(names(sdr$value)=='D')],
                  Z=sdr$sd[which(names(sdr$value)=='Z')]*ScaleFac,
                  u=sdr$sd[which(names(sdr$value)=='u')]/ScaleFac,
                  R=matrix(sdr$sd[which(names(sdr$value)=='FullCovMat')],nrow=length(logsdObs),ncol=length(logsdObs)))}
    
    #Compute AIC.
    AIC<-2*length(opt1$par) + 2*opt1$value;AIC
    
    #print(AIC)
    if(EstSE){return(list(Optimization = opt1, Estimates = pl1, Fits = FitSeries,AIC=AIC,StdErr=SES,ParCorrs=sdr$cov.fixed))
    }else{return(list(Optimization = opt1, Estimates = pl1, Fits = FitSeries,AIC=AIC))}
    
  }
  #
  # runDFA(ZL_FW1,ErrStruc='UNC')
  #
  # newmethodR<-pl1$R
  
  ## load the data (there are 3 datasets contained here)
  data(lakeWAplankton, package = "MARSS")
  ## we want lakeWAplanktonTrans, which has been transformed
  ## so the 0s are replaced with NAs and the data z-scored
  all_dat <- lakeWAplanktonTrans
  ## use only the 10 years from 1980-1989
  yr_frst <- 1980
  yr_last <- 1989
  plank_dat <- all_dat[all_dat[, "Year"] >= yr_frst & all_dat[, "Year"] <= yr_last, ]
  ## create vector of phytoplankton group names
  phytoplankton <- c("Cryptomonas", "Diatoms", "Greens", "Unicells", "Other.algae")
  ## get only the phytoplankton
  dat_1980 <- plank_dat[, phytoplankton]
  ## transpose data so time goes across columns
  dat_1980 <- t(dat_1980)
  ## get number of time series
  N_ts <- dim(dat_1980)[1]
  ## get length of time series
  TT <- dim(dat_1980)[2]
  dat <- MARSS::zscore(dat)
  
  library(MARSS)
  m1.em<-MARSS(dat,model=list(R='unconstrained',m=1,tinitx=1),form='dfa',z.score=FALSE, control=list(minit=1500, maxit=1500))
  m1.bfgs<-MARSS(dat,model=list(R='unconstrained',m=1,x0="zero"),form='dfa',z.score=FALSE, method="BFGS")
  mym1<-runDFA(obs=dat,NumStates=1,ErrStruc='UNC');
  cbind(coefficients(m1.em)$R,as.vector(mym1$Estimates$R[lower.tri(mym1$Estimates$R,diag=T)]))
  cbind(coefficients(m1.em)$Z,as.vector(mym1$Estimates$Z)) 
  
  c(logLik(m1.em), logLik(m1.bfgs), dfaLL(mym1))
  
  mod_list = list(m = 3, R = "diagonal and unequal")
  dfa_temp <- MARSS(dat, model = mod_list, form = "dfa", z.score = FALSE, 
                    control = con_list, covariates = temp)
  
  small<-ZL_FW1[sample.int(4,n=8,replace=F),]
  small[sample.int(size=100,n=nrow(small)*ncol(small),replace=F)]<-NA

  small<-dat.class[sample.int(5,n=nrow(dat.class),replace=F),]
  #library(MARSS)
  m1<-MARSS(small,model=list(R='unconstrained',m=1),form='dfa',z.score=T)#,control=list(MCInit=T));
  mym1<-runDFA(obs=small,NumStates=1,ErrStruc='UNC');
  cbind(coefficients(m1)$R,as.vector(mym1$Estimates$R[lower.tri(mym1$Estimates$R,diag=T)]))#mym1$Estimates$Variances

  #Check MARSS vs TMB model fits
  par(mfrow=c(2,2))
  R2S<-matrix(NA,nrow=2,ncol=4)
  for(i in 1:4){
  	R1<-na.omit(cbind(small[i,],coefficients(m1)$Z[i]*m1$states[1,],mym1$Fits[i,]))
  	R2_1_MARSS<-summary(lm(R1[,1]~R1[,2]))$r.squared;R2S[1,i]<-R2_1_MARSS
  	R2_1_TMB<-summary(lm(R1[,1]~R1[,3]))$r.squared;R2S[2,i]<-R2_1_TMB
  	plot(small[i,],pch=16,ylab=paste('Series ',i,sep=''))
  	points(coefficients(m1)$Z[i]*m1$states[1,],type='l',col='blue',lwd=2)
  	points(mym1$Fits[i,],type='l',col='red',lwd=2)
  	legend('bottom',legend=c(as.character(round(R2_1_MARSS,2)),as.character(round(R2_1_TMB,2))),col=c('blue','red'),bty='n',lty=rep(1,2),lwd=rep(2,2))
  }
  mtext(paste('MARSS = ',as.character(round(mean(R2S[1,]),2))),3,outer=T,col='blue',line=-2.5)
  mtext(paste('TMB = ',as.character(round(mean(R2S[2,]),2))),3,outer=T,col='red',line=-3.5)