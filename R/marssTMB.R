#' Parameter estimation using TMB
#'
#' Minimal error checking is done in this function.  
#' Normal calling function is `MARSS.tmb()` which in
#' turn calls this function. Restrictions
#' 
#' * No time-varying parameters
#' * Currently only DFA models are coded up
#' * x0 and V0 fixed (stochastic prior)
#' * Q is fixed (not estimated)
#' 
#' @param MLEobj A properly formatted MARSS model as output by MARSS.tmb()

#' @return A list with Optimization, Estimates, Fits, and AIC
#' @example inst/examples/dfa_example.R
#' @author This function is based off of dfaTMB.R written by Tim Cline while a graduate student in the Fish 507 Time Series Analysis course. Eli Holmes later modified it to replicate the MARSS(x, form="dfa") model.
#' @export
MARSStmb <- function(MLEobj) {
  MODELobj <- MLEobj[["model"]]
  ty <- t(MODELobj[["y"]])
  model.dims <- attr(MODELobj, "model.dims")
  n <- model.dims[["data"]][1]
  TT <- model.dims[["data"]][2]
  m <- model.dims[["x"]][1]
  Covars <- coef(fit, type="matrix")$c
  if(ncol(Covars)==1) Covars <- matrix(Covars, nrow=nrow(Covars), ncol=TT)
  
  # Set up the initial matrices
  eleminits <- list()
  for(elem in c("Z", "D", "R", "Q", "V0")){
    eleminits[[elem]] <- coef(fit, type="matrix", what="start")[[elem]]
  }
  # Set up the maps
  elemmaps <- list()
  for(elem in c("Z", "D")){
    elemmaps[[elem]] <- create.elem.maps(MLEobj, elem=elem)[["map"]]
  }
  for(elem in c("R", "Q")){
    elemmaps[[elem]] <- list(
      diag = create.varcov.maps(MLEobj, elem=elem)[["map.diag"]],
      offdiag = create.varcov.maps(MLEobj, elem=elem)[["map.offdiag"]]
    )
  }

  # Creates the input data list
  data <- list(
    model = "dfa", 
    obs = ty,
    Covar = coef(fit, type="matrix")[["c"]]
  )
  
  # Creates the input parameter list
  R <- eleminits[["R"]]
  parameters <- list(
    logsdObs = log(diag(R)),
    cholCorr = chol(R)[upper.tri(R)],
    covState = eleminits[["Q"]],
    covinitState = eleminits[["V0"]],
    D = eleminits[["D"]],
    Z = eleminits[["Z"]],
    u = matrix(0, nrow = TT, ncol = m)
  )
  
  # Create the map (mask) that indicates what parameters to not estimate
  # Note x0 and V0 are fixed (stochastic prior)
  maplist <- list(
    logsdObs = elemmaps[["R"]][["diag"]], 
    cholCorr = elemmaps[["R"]][["offdiag"]], 
    covState = factor(matrix(NA, nrow = m, ncol = m)), 
    covinitState = factor(matrix(NA, nrow = m, ncol = m)),
    D = elemmaps[["D"]], 
    Z = elemmaps[["Z"]]
    )
  
  # Creates the model object and runs the optimization
  obj1 <- TMB::MakeADFun(
    data, parameters, random = "u",
    DLL = "marssTMB_TMBExports",
    silent = silent, map = maplist
  )
  
  
}
