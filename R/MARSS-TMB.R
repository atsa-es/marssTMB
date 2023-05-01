#' Main marssTMB function
#' 
#' Fit a MARSS model with TMB + optimization. TMB does fast computation of the
#' likelihood gradients and this is used
#' with [stats::nlminb()] or [stats::optim()] to optimize the likelihood. The 
#' former is much faster than the later and is the default.
#' 
#' Note. Goal is to have the user call with the MARSS package as 
#' `MARSS(y, method="tmb")`, which then calls `MARSStmb()`. This is the mimics
#' the behavior of `method="BFGS"` and `method="kem"` in [MARSS::MARSS()] which
#' looks for a fitting function called `MARSSxyz`, where `xyz` is the method.
#' Further arguments for the optimization method can be passed into `control`.
#' 
#' @param y Vector of observations n x T.
#' @param model list with 
#'    * R  "diagonal and equal", "unconstrained", "diagonal and unequal"
#'    * m the number of states (factors). default is 1
#' @param inits list of initial conditions
#' @param method to pass to optim call; ignored for `fun="nlminb"`
#' @param form The equation form used in the marssTMB() call. The default is "dfa". 
#' @param fit Whether to fit the model.
#' @param silent Show TMB output when fitting
#' @param control list for the optimization function. [stats::nlminb()] or [stats::optim()], `control$fun.opt` allows you to choose optim or nlminb as the optimization function. `control$optim.method` allows you to choose method for `optim()`.
#' 
#' @return The output list from [MARSStmb()]
#' @example inst/examples/MARSS_TMB_example.R
#' @author Eli Holmes
#' @export
MARSS_tmb <- function(y,
                     model = NULL,
                     inits = NULL,
                     miss.value = as.numeric(NA),
                     method = "TMB",
                     form = "dfa",
                     fit = TRUE,
                     silent = FALSE,
                     control = NULL,
                     ...) {
  pkg <- "marssTMB"
  method <- match.arg(method)
  # This section is temporary
  # Needed now because currently MARSS does not allow fun.opt in control
  if(is.null(control$fun.opt)){
    fun.opt <- "optim"
    optim.method <- ifelse(is.null(control$optim.method), "BFGS", control$optim.method)
  }else{
    fun.opt <- control$fun.opt
  }
  control <- control[!(names(control) %in% c("fun.opt", "optim.method"))]
  
  # Set up a MARSS() call list and call to get a properly set up MARSS model
  MARSS.call <- list(y = y, inits = inits, model = model, control = control, method = "BFGS", form = form, silent = TRUE, fit = FALSE, ...)
  # x is now a marssMLE object but no par element since not fit
  x <- do.call(MARSS::MARSS, MARSS.call)
  
  # Set up the control list a TMB fit
  # Temporary until method="tmb" is added to MARSS
  x[["method"]] <- "TMB"
  allowed.fun.opt = c("optim", "nlminb")
  if (!(fun.opt %in% allowed.fun.opt)) {
    stop(paste("control$fun.opt must be one of:", paste(allowed.fun.opt, collapse=", ")))
  }
  if(fun.opt == "optim"){
    allowed.optim.methods <- c("BFGS")
    # Some error checks depend on an allowable method
    if (!(optim.method %in% allowed.optim.methods)) {
      stop(paste("control$optim.method must be one of:", paste(allowed.optim.methods, collapse=", ")))
    }
  }
  control[["fun.opt"]] <- fun.opt
  if(fun.opt == "optim") control[["optim.method"]] <- optim.method
  control[["tmb.silent"]] <- TRUE
  # set up control defaults
  if(fun.opt == "nlminb"){
    allowed.in.opt.fun <- c("eval.max", "iter.max", "trace", "abs.tol", "rel.tol", "x.tol", "xf.tol", "step.min", "step.max", "sing.tol", "scale.init", "diff.g")
    if(is.null(control$iter.max)) control$iter.max = control$maxit
    if(is.null(control$eval.max)) control$eval.max = control$maxit
  }
  if(fun.opt == "optim"){
    allowed.in.opt.fun <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol", "alpha", "beta", "gamma", "REPORT", "warn.1d.NelderMead", "type", "lmm", "factr", "pgtol", "temp", "tmax")
    if(is.null(control$reltol)) control$reltol = 1e-12
    if(is.null(control$maxit)) control$maxit = control$maxit
  }
  allowed.in.control <- c("fun.opt", "optim.method", "tmb.silent", "silent")
  control <- control[names(control) %in% c("fun.opt", "optim.method", "tmb.silent", allowed.in.opt.fun)]
  x[["control"]] <- control
  
  
  # Error check for the DFA model
  # This section is temporary. Currently only DFA model is allowed.
  model.descrip <- MARSS:::describe.marssMODEL(x$model)
  if(form == "dfa"){
  is.unconstrained <- function(elem) substr(model.descrip[[elem]], 1, 5) == "uncon"
  is.diagonal <- function(elem) substr(model.descrip[[elem]], 1, 8) == "diagonal"
  is.fixed <- function(elem) substr(model.descrip[[elem]], 1, 5) == "fixed"
  is.identity <- function(elem){
    val <- model.descrip[[elem]]
    substr(val, 1, 8) == "identity" | val == "fixed and all one (1 x 1)"
  }
  is.zero <- function(elem) substr(model.descrip[[elem]], 1, 14) == "fixed and zero"
  
  # Check that R is allowed
  elem <- "R"
  ok <- is.diagonal(elem) | is.fixed(elem) | is.unconstrained(elem)
  if(!ok){
    stop(paste0(pkg, ": R must be diagonal, fixed or unconstrained"))
  }

  # Check that Q and B are identity
  for(elem in c("B", "Q")){
  ok <- is.identity(elem)
  if(!ok) stop(paste0(pkg, ": ", elem, " must be identity"))
  }
  
  # Check that u, a, and C are zero
  for(elem in c("U", "A", "C")){
  ok <- is.zero(elem)
  if(!ok) stop(paste0(pkg, ": ", elem, " must be zero"))
  }

  # Check that x0 and V0 are fixed
  for(elem in c("x0", "V0")){
    ok <- is.fixed(elem)
    if(!ok) stop(paste0(pkg, ": ", elem, " must be fixed"))
  }
  }
  
  obj <- MARSStmb(x)
  
  return(obj)
}
