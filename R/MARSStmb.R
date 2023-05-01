#' Parameter estimation using TMB
#' 
#' Status 4/30/23 fits model. Now work on [to_marssMLE()] to format a [MARSS::marssMLE object].
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

#' @return A [MARSS::marssMLE] object
#' @example inst/examples/dfa_example.R
#' @author Eli Holmes. This function is inspired by dfaTMB.R written by Tim Cline while a graduate student in the Fish 507 Time Series Analysis course. 
#' @seealso [MARSS::MARSSoptim()], [MARSS::MARSSkem()]
#' @export
MARSStmb <- function(MLEobj) {
  MODELobj <- MLEobj[["model"]]
  ty <- t(MODELobj[["data"]])
  model.dims <- attr(MODELobj, "model.dims")
  n <- model.dims[["data"]][1]
  TT <- model.dims[["data"]][2]
  m <- model.dims[["x"]][1]
  Covars <- MODELobj[["fixed"]][["c"]]
  control <- MLEobj[["control"]]
  fun.opt <- control[["fun.opt"]]
  tmb.silent <- control[["tmb.silent"]]
  optim.method <- control[["optim.method"]]
  control <- control[!(names(control) %in% c("fun.opt", "optim.method", "tmb.silent"))]
  # Expand out to full covariate matrix
  if(ncol(Covars)==1) Covars <- matrix(Covars, nrow=nrow(Covars), ncol=TT)
  
  # Set up the initial matrices
  eleminits <- list()
  for(elem in c("Z", "D", "R", "Q", "V0", "x0")){
    eleminits[[elem]] <- coef(MLEobj, type="matrix", what="start")[[elem]]
  }
  # Set up the maps
  elemmaps <- list()
  for(elem in c("Z", "D", "x0")){
    elemmaps[[elem]] <- create.elem.maps(MLEobj, elem=elem)[["map"]]
  }
  # maps for var-cov matrices have diagonal separate from off-diagonal
  for(elem in c("R", "Q", "V0")){
    elemmaps[[elem]] <- list(
      diag = create.varcov.maps(MLEobj, elem=elem)[["map.diag"]],
      offdiag = create.varcov.maps(MLEobj, elem=elem)[["map.offdiag"]]
    )
  }

  # Creates the input data list
  # For now, we will assume V0 is a fixed (diagonal) matrix
  data <- list(
    model = "marxss", 
    obs = ty,
    Covar = Covars,
    V0 = eleminits[["V0"]]
  )
  
  # Note x0 and V0 are fixed (stochastic prior) for DFA
  # But x0 might be estimated in the future so is here as well

    # Creates the list of initial (start) values of parameter list
  R <- eleminits[["R"]]
  parameters <- list(
    logsdObs = log(diag(R)), # log of diagonal of R
    cholCorr = chol(R)[upper.tri(R)], # off-diagonal of chol of R
    covState = eleminits[["Q"]],
    covinitState = eleminits[["V0"]],
    D = eleminits[["D"]],
    Z = eleminits[["Z"]],
    x0 = eleminits[["x0"]],
    u = matrix(0, nrow = TT, ncol = m) # states
  )
  
  # Create the map (mask) that indicates what parameters to not estimate
  # the states are treated as a parameter but are all estimated
  # so do not appear here
  maplist <- list(
    logsdObs = elemmaps[["R"]][["diag"]], 
    cholCorr = elemmaps[["R"]][["offdiag"]], 
    covState = factor(matrix(NA, nrow = m, ncol = m)), 
    covinitState = factor(matrix(NA, nrow = m, ncol = m)),
    D = elemmaps[["D"]], 
    Z = elemmaps[["Z"]],
    x0 = elemmaps[["x0"]]
    )
  
  # Creates the model object and runs the optimization
  obj1 <- TMB::MakeADFun(
    data, parameters, random = "u",
    DLL = "marssTMB_TMBExports",
    silent = MLEobj[["control"]][["tmb.silent"]], map = maplist
  )
  
  # Optimization
  if(MLEobj[["control"]][["fun.opt"]] == "nlminb"){
    opt1 <- stats::nlminb(obj1$par, obj1$fn, obj1$gr, control = control)
    # add output also found in optim output
    opt1$value <- opt1$objective
    opt1$counts <- opt1$evaluations
    # add this on for printing
    MLEobj[["control"]][["maxit"]] <-  MLEobj[["control"]][["iter.max"]]
  }
  if(MLEobj[["control"]][["fun.opt"]] == "optim"){
    opt1 <- stats::optim(obj1$par, obj1$fn, gr=obj1$gr, control = control)
    #obj1$control <- MLEobj[["control"]]
    #opt1 <- do.call("optim", obj1)
    opt1$objective <- opt1$value
    opt1$iterations <- opt1$counts[2]
  }
  
  # Add names to the par output
  model.elem <- names(MODELobj[["fixed"]])
  for(elem in model.elem){
    names(opt1$par)[names(opt1$par)==elem] <- paste0(elem, ".", levels(obj1$env$map[[elem]]))
  }
  
  return(list(obj.function = obj1, opt.output = opt1, MLEobj = MLEobj))

}
