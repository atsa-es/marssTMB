#' Fit a DFA model with TMB.
#'
#' @param obs Vector of observations n x T.
#' @param NumStates m
#' @param ErrStruc DE=diagonal and equal, UNC=unconstrained, DUE=diagonal and unequal
#' @param EstCovar TRUE/FALSE
#' @param Covars
#' @param indivCovar
#' @param Dmat
#' @param Dfac
#' @param Rfac
#' @param logsdObs
#' @param logsdObsFac
#' @param cholCorr
#' @param cholFac
#' @param EstSE
#' @param silent Show TMB output when fitting

#' @return A list with Optimization, Estimates, Fits, and AIC
#' @example inst/examples/dfa_example.R
#' @author Tim Cline wrote most of this while a graduate student in the Fish 507 Time Series Analysis course. Eli Holmes later modified it to replicate the MARSS(x, form="dfa") model.
#' @export
dfaTMB <- function(obs, 
                   NumStates = 1, 
                   ErrStruc = "DE", EstCovar = FALSE,
                   Covars = NULL, indivCovar = FALSE, 
                   Dmat = NULL, Dfac = NULL,
                   Rfac = NULL, logsdObs = NULL, 
                   logsdObsFac = NULL,
                   cholCorr = NULL, cholFac = NULL, 
                   EstSE = FALSE,
                   silent = TRUE) {
  obs <- t(obs)
  # creates the Z factor to fix the upper corner to 0
  Zfac <- ZmatFactorGen(Data = obs, NumStates = NumStates) 
  if (EstCovar) { 
    # If you are estimating covariates this creates 
    # Dmat and Dfac based on the data, number of covars
    if (is.null(Dmat) & is.null(Dfac)) { 
      # This checks that you did not supply Dmat or Dfac 
      # (Manual covariate parameter entries)
      if (!indivCovar) {
        Dmat <- matrix(rep(0, ncol(obs) * nrow(Covars)), ncol = nrow(Covars), nrow = ncol(obs))
        Dfac <- as.factor(seq(1, ncol(obs) * nrow(Covars)))
      } else {
        Dmat <- matrix(0, ncol = nrow(Covars), nrow = ncol(obs))
        diag(Dmat) <- rep(0, nrow(Covars)) # rnorm(nrow(Covars),0,1)
        Dfac <- matrix(NA, ncol = nrow(Covars), nrow = ncol(obs))
        diag(Dfac) <- seq(1, nrow(Covars))
        Dfac <- as.factor(Dfac)
      }
    }
    data <- list(
      model = "dfa",
      obs = obs, 
      NumState = NumStates, 
      Covar = Covars)
  } else { 
    # If you are not estimating covariates we just pass a time series of 0's, 
    # a Dmat of 0's, and and NA factors so the paramters will not be estimated

    Dmat <- matrix(0, ncol = 1, nrow = ncol(obs))
    Dfac <- as.factor(rep(NA, ncol(obs)))
    Covars <- matrix(0, nrow = 1, ncol = nrow(obs))
    data <- list(
      model = "dfa", obs = obs,
      NumState = NumStates,
      Covar = Covars
    )
  }

  # This set of if-else statements creates the proper parameter set for the error structure selected
  if (is.null(logsdObs) & is.null(logsdObsFac) & is.null(cholCorr) & is.null(cholFac)) {
    if (ErrStruc == "DE") {
      cholCorr <- rep(0, ncol(obs) * (ncol(obs) - 1) / 2)
      logsdObs <- log(rep(0.5, ncol(obs)))
      logsdObsFac <- rep(1, ncol(obs))
      logsdObsFac <- factor(logsdObsFac)
      cholFac <- rep(NA, ncol(obs) * (ncol(obs) - 1) / 2)
      cholFac <- factor(cholFac)
    } else if (ErrStruc == "DUE") {
      cholCorr <- rep(0, ncol(obs) * (ncol(obs) - 1) / 2)
      logsdObs <- log(rep(0.5, ncol(obs)))
      logsdObsFac <- seq(1, ncol(obs))
      logsdObsFac <- factor(logsdObsFac)
      cholFac <- rep(NA, ncol(obs) * (ncol(obs) - 1) / 2)
      cholFac <- factor(cholFac)
    } else if (ErrStruc == "UNC") {
      cholCorr <- rep(0, ncol(obs) * (ncol(obs) - 1) / 2)
      logsdObs <- log(rep(0.5, ncol(obs)))
      logsdObsFac <- seq(1, ncol(obs))
      logsdObsFac <- factor(logsdObsFac)
      cholFac <- seq(1, (ncol(obs) * (ncol(obs) - 1) / 2))
      cholFac <- factor(cholFac)
    }
  }

  # Creates the input parameter list
  parameters <- list(
    logsdObs = logsdObs,
    cholCorr = cholCorr,
    covState = diag(1, NumStates),
    covinitState = diag(5, NumStates),
    D = Dmat,
    Z = ZmatGen(Data = obs, NumStates = NumStates),
    u = matrix(0, nrow = nrow(obs), ncol = NumStates)
  )
  covinitStateFac <- factor(matrix(NA, nrow = NumStates, ncol = NumStates))
  covStateFac <- factor(matrix(NA, nrow = NumStates, ncol = NumStates))
  maplist <- list(Z = Zfac, D = Dfac, cholCorr = cholFac, logsdObs = logsdObsFac, covState = covStateFac, covinitState = covinitStateFac)

  # Creates the model object and runs the optimization
  obj1 <- TMB::MakeADFun(data,
    parameters,
    random = "u",
    DLL = "marssTMB_TMBExports",
    silent = silent,
    map = maplist
  )
  opt1 <- stats::nlminb(obj1$par, obj1$fn, obj1$gr, control = list(iter.max = 2000, eval.max = 2000))
  # newtonOption(obj1, smartsearch=TRUE)
  obj1$control <- list(trace = 1, REPORT = 1, reltol = 1e-12, maxit = 2000)
  obj1$fn()
  obj1$gr()
  # obj1$method='BFGS'
  obj1$par <- opt1$par
  system.time(opt1 <- do.call("optim", obj1))
  pl1 <- obj1$env$parList() # This contains all of your parameter estimates RAW as they come out of the optimizer
  if (EstSE) {
    sdr <- sdreport(obj1)
  }

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
  pl1$u <- t(pl1$u)

  if (EstSE) {
    pl1$R <- matrix(sdr$value[which(names(sdr$value) == "FullCovMat")], nrow = length(logsdObs), ncol = length(logsdObs))
  } else {
    pl1$R <- diag(exp(pl1$logsdObs)) %*% obj1$report()$FullCorrMat %*% diag(exp(pl1$logsdObs))
  }

  # Fits for each time series
  FitSeries <- pl1$Z %*% pl1$u + pl1$D %*% Covars

  # Standard Errors for parameters
  if (EstSE) {
    SES <- list(
      D = sdr$sd[which(names(sdr$value) == "D")],
      Z = sdr$sd[which(names(sdr$value) == "Z")] * ScaleFac,
      u = sdr$sd[which(names(sdr$value) == "u")] / ScaleFac,
      R = matrix(sdr$sd[which(names(sdr$value) == "FullCovMat")], nrow = length(logsdObs), ncol = length(logsdObs))
    )
  }

  # Compute AIC.
  AIC <- 2 * length(opt1$par) + 2 * opt1$value
  AIC

  if (EstSE) {
    res <- list(Optimization = opt1, Estimates = pl1, Fits = FitSeries, AIC = AIC, StdErr = SES, ParCorrs = sdr$cov.fixed)
  } else {
    res <- list(Optimization = opt1, Estimates = pl1, Fits = FitSeries, AIC = AIC)
  }
  class(res) <- c("marssTMB", class(res))
  return(res)
}

# This function generates the Z matrix based on the number of times series and the number of states that are estimated. It is called from with the run_dfa function.
ZmatGen <- function(Data, NumStates) {
  tZ <- matrix(0.5, nrow = ncol(Data), ncol = NumStates)
  if (NumStates > 1) {
    for (i in 1:(NumStates - 1)) {
      tZ[i, (NumStates - (NumStates - 2) + (i - 1)):NumStates] <- rep(0, NumStates - 1 - (i - 1))
    }
    return(tZ)
  } else {
    return(tZ)
  }
}

# This function creates a Z mat factor which allows TMB to fix certain parameters. 
# This is required as the upper triangle of the Z matrix must fixed at 0 to allow the model to be identifiable. 
ZmatFactorGen <- function(Data, NumStates) {
  tZ <- matrix(seq(1, ncol(Data) * NumStates), nrow = ncol(Data), ncol = NumStates)
  if (NumStates > 1) {
    for (i in 1:(NumStates - 1)) {
      tZ[i, (NumStates - (NumStates - 2) + (i - 1)):NumStates] <- rep(NA, NumStates - 1 - (i - 1))
    }
    tZ[!is.na(tZ)] <- seq(1, sum(!is.na(tZ)))
    return(as.factor(tZ))
  } else {
    return(as.factor(tZ))
  }
}

# Compute AIC for my DFA model output
dfaAIC <- function(x, AICc = F) {
  opt <- x[["Optimization"]]
  NumPar <- length(opt$par)
  NLL <- opt$value # opt$value
  AIC <- 2 * NumPar + 2 * NLL
  if (AICc) {
    AIC <- AIC + (2 * NumPar * (NumPar + 1)) / (nrow(x[["Estimates"]]$Z) * ncol(x[["Estimates"]]$u) - NumPar - 1)
  }
  return(AIC)
}

# Compute log likelihood
dfaLL <- function(x) {
  opt <- x[["Optimization"]]
  NumPar <- length(opt$par)
  NLL <- opt$value # opt$value
  return(-1 * NLL)
}
