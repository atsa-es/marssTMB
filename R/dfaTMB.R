#' Fit a DFA model with TMB.
#'
#' @param y Vector of observations n x T.
#' @param model list with 
#'    * R  "diagonal and equal", "unconstrained", "diagonal and unequal"
#'    * m number of states (x)
#' @param inits list of initial conditions
#' @param EstCovar TRUE/FALSE
#' @param Covars
#' @param indivCovar
#' @param Dmat
#' @param Dfac
#' @param EstSE
#' @param silent Show TMB output when fitting
#' @param fun.opt function to use for optimization: `stats::nlminb()` or `stats::optim()`
#' @param method to pass to optim call; ignored for `fun="nlminb"`
#' @param form The equation form used in the marssTMB() call. The default is "dfa". 

#' @return A list with Optimization, Estimates, Fits, and AIC
#' @example inst/examples/dfa_example.R
#' @author Tim Cline wrote most of this while a graduate student in the Fish 507 Time Series Analysis course. Eli Holmes later modified it to replicate the MARSS(x, form="dfa") model.
#' @export
dfaTMB <- function(y, 
                   model = list(m = 1, R="diagonal and equal"),
                   inits = list(R=0.05),
                   EstCovar = FALSE,
                   Covars = NULL, indivCovar = FALSE, 
                   Dmat = NULL, Dfac = NULL,
                   EstSE = FALSE,
                   silent = TRUE,
                   fun.opt = c("nlminb", "optim", "nlminb+optim"),
                   method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
                   control = NULL,
                   form = c("dfa", "marxss")){
  ty <- t(y)
  m <- model$m
  n <- ncol(ty)
  TT <- nrow(ty)
  fun.opt <- match.arg(fun.opt)
  method <- match.arg(method)
  # set up control defaults
  if(fun.opt == "nlminb"){
    if(is.null(control$iter.max)) control$iter.max = 2000
    if(is.null(control$eval.max)) control$eval.max = 2000
  }
  if(fun.opt == "optim" | fun.opt == "nlminb+optim"){
    if(is.null(control$reltol)) control$reltol = 1e-12
    if(is.null(control$maxit)) control$maxit = 2000
  }
  # creates the Z factor to fix the upper corner to 0
  Zfac <- ZmatFactorGen(n, m) 
  if (EstCovar) { 
    # If you are estimating covariates this creates 
    # Dmat and Dfac based on the data, number of covars
    if (is.null(Dmat) & is.null(Dfac)) { 
      # This checks that you did not supply Dmat or Dfac 
      # (Manual covariate parameter entries)
      if (!indivCovar) {
        Dmat <- matrix(rep(0, n * nrow(Covars)), ncol = nrow(Covars), nrow = n)
        Dfac <- as.factor(seq(1, n * nrow(Covars)))
      } else {
        Dmat <- matrix(0, ncol = nrow(Covars), nrow = n)
        #diag(Dmat) <- rep(0, nrow(Covars)) # rnorm(nrow(Covars),0,1)
        Dfac <- matrix(NA, ncol = nrow(Covars), nrow = n)
        # EEH: this assumes you have n covariates
        diag(Dfac) <- seq(1, nrow(Covars))
        Dfac <- as.factor(Dfac)
      }
    }
    data <- list(
      model = "dfa",
      obs = ty, 
      Covar = Covars)
  } else { 
    # If you are not estimating covariates we just pass a time series of 0's, 
    # a Dmat of 0's, and and NA factors so the paramters will not be estimated

    Dmat <- matrix(0, ncol = 1, nrow = n)
    Dfac <- as.factor(rep(NA, n))
    Covars <- matrix(0, nrow = 1, ncol = TT)
    data <- list(
      model = "dfa", 
      obs = ty,
      Covar = Covars
    )
  }

  # Creates the proper parameter set for the error structure selected
    cholCorr <- rep(0, n * (n - 1) / 2)
    logsdObs <- log(rep(sqrt(inits$R), n))
    if (model$R == "diagonal and equal") {
      logsdObsFac <- rep(1, n)
      cholFac <- rep(NA, n * (n - 1) / 2)
    } else if (model$R == "diagonal and unequal") {
      logsdObsFac <- seq(1, n)
      cholFac <- rep(NA, n * (n - 1) / 2)
    } else if (model$R == "unconstrained") {
      logsdObsFac <- seq(1, n)
      cholFac <- seq(1, (n * (n - 1) / 2))
    }
    logsdObsFac <- factor(logsdObsFac)
    cholFac <- factor(cholFac)

  # Creates the input parameter list
  parameters <- list(
    logsdObs = logsdObs,
    cholCorr = cholCorr,
    covState = diag(1, m),
    covinitState = diag(5, m),
    D = Dmat,
    Z = ZmatGen(n, m),
    u = matrix(0, nrow = TT, ncol = m)
  )
  covinitStateFac <- factor(matrix(NA, nrow = m, ncol = m))
  covStateFac <- factor(matrix(NA, nrow = m, ncol = m))
  maplist <- list(Z = Zfac, D = Dfac, cholCorr = cholFac, logsdObs = logsdObsFac, covState = covStateFac, covinitState = covinitStateFac)

  # Creates the model object and runs the optimization
  obj1 <- TMB::MakeADFun(data,
    parameters,
    random = "u",
    DLL = "marssTMB_TMBExports",
    silent = silent,
    map = maplist
  )

  # Optimization
  if(fun.opt == "nlminb"){
    opt1 <- stats::nlminb(obj1$par, obj1$fn, obj1$gr, control = control)
  }
  if(fun.opt == "optim"){
    obj1$control <- control
    opt1 <- do.call("optim", obj1)
    opt1$objective <- opt1$value
  }
  if(fun.opt == "nlminb+optim"){
    opt1 <- stats::nlminb(obj1$par, obj1$fn, obj1$gr, control = list(iter.max = 2000, eval.max = 2000))
    obj1$par <- opt1$par
    obj1$control <- control
    opt1 <- do.call("optim", obj1)
    opt1$objective <- opt1$value
  }
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
  # AIC <- 2 * length(opt1$par) + 2 * opt1$value
  AIC <- 2 * length(opt1$par) + 2 * opt1$objective
  AIC
  
  

  if (EstSE) {
    res <- list(Optimization = opt1, Estimates = pl1, Fits = FitSeries, AIC = AIC, StdErr = SES, ParCorrs = sdr$cov.fixed, logLik = -1*opt1$objective)
  } else {
    res <- list(Optimization = opt1, Estimates = pl1, Fits = FitSeries, AIC = AIC, logLik = -1*opt1$objective)
  }
  
  # Set up the return list to match marssMLE objects
  res$method <- "TMB"
  class(res) <- c("marssTMB", class(res))
  return(res)
}

# This function generates the Z matrix based on the number of times series and the number of states that are estimated. It is called from with the run_dfa function.
ZmatGen <- function(n, m) {
  tZ <- matrix(0.5, nrow = n, ncol = m)
  if (m > 1) {
    for (i in 1:(m - 1)) {
      tZ[i, (m - (m - 2) + (i - 1)):m] <- rep(0, m - 1 - (i - 1))
    }
    return(tZ)
  } else {
    return(tZ)
  }
}

# This function creates a Z mat factor which allows TMB to fix certain parameters. 
# This is required as the upper triangle of the Z matrix must fixed at 0 to allow the model to be identifiable. 
ZmatFactorGen <- function(n, m) {
  tZ <- matrix(seq(1, n * m), nrow = n, ncol = m)
  if (m > 1) {
    for (i in 1:(m - 1)) {
      tZ[i, (m - (m - 2) + (i - 1)):m] <- rep(NA, m - 1 - (i - 1))
    }
    tZ[!is.na(tZ)] <- seq(1, sum(!is.na(tZ)))
    return(as.factor(tZ))
  } else {
    return(as.factor(tZ))
  }
}

# Compute AIC for my DFA model output
dfaAIC <- function(x, AICc = FALSE) {
  opt <- x[["Optimization"]]
  NumPar <- length(opt$par)
  NLL <- opt$objective # opt$value (optim)
  AIC <- 2 * NumPar + 2 * NLL
  if (AICc) {
    AIC <- AIC + (2 * NumPar * (NumPar + 1)) / (nrow(x[["Estimates"]]$Z) * ncol(x[["Estimates"]]$u) - NumPar - 1)
  }
  return(AIC)
}

# Compute log likelihood
dfaLL <- function(x) {
  opt <- x[["Optimization"]]
  NLL <- opt$objective # opt$value
  return(-1 * NLL)
}
