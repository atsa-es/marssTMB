#' Main marssTMB function
#' 
#' Fit a MARSS model with TMB + optimization. TMB does fast computation of the
#' likelihood gradients and this is used
#' with stats::nlminb()` or `stats::optim()` to optimize the likelihood. The 
#' former is much faster than the later and is the default.
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
#' @param control list for the optimization function. `stats::nlminb()` or `stats::optim()`, control$fun.opt allows you to choose optim or nlminb as the optimization function. control$optim.method allows you to choose method for `optim()`.
#' 
#' @return A list with Optimization, Estimates, Fits, and AIC
#' @example inst/examples/marssTMB_example.R
#' @author Eli Holmes
#' @export
MARSS.tmb <- function(y,
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
  allowed.fun.opt = c("optim", "nlminb")
  fun.opt <- ifelse(is.null(control$fun.opt), "optim", control$fun.opt)
  if (!fun.opt %in% allowed.fun.opt) {
    stop(paste("control$fun.opt must be one of:", fun.opt))
  }
  if(fun.opt == "optim"){
    optim.method <- ifelse(is.null(control$optim.method), "BFGS", control$optim.method)
    allowed.optim.methods <- c("BFGS")
    # Some error checks depend on an allowable method
    if (!optim.method %in% allowed.optim.methods) {
      stop(paste("control$optim.method must be one of:", allowed.optim.methods))
    }
  }
  
  # Check the model forms
  
  
  MARSS.call <- list(y = y, inits = inits, model = model, control = control, method = "kem", form = form, silent = silent, fit = FALSE, ...)
  
  x <- do.call(MARSS::MARSS, MARSS.call)
  x$method <- "TMB"
  
  # Error check for the DFA model
  if(form == "dfa"){
  is.unconstrained <- function(elem) substr(MARSS:::describe.marssMODEL(x$model)[[elem]], 1, 5) == "uncon"
  is.diagonal <- function(elem) substr(MARSS:::describe.marssMODEL(x$model)[[elem]], 1, 8) == "diagonal"
  is.fixed <- function(elem) substr(MARSS:::describe.marssMODEL(x$model)[[elem]], 1, 5) == "fixed"
  is.identity <- function(elem){
    val <- MARSS:::describe.marssMODEL(x$model)[[elem]]
    substr(val, 1, 8) == "identity" | val == "fixed and all one (1 x 1)"
  }
  is.zero <- function(elem) substr(MARSS:::describe.marssMODEL(x$model)[[elem]], 1, 14) == "fixed and zero"
  
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
  
  return(x)
}
