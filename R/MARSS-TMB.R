#' Main marssTMB function
#' 
#' Fit a MARSS model with TMB + optimization. TMB does fast computation of the
#' likelihood gradients and this is used
#' with [stats::nlminb()] or [stats::optim()] to optimize the likelihood. The 
#' former is much faster than the later and is the default.
#' 
#' Note in the development version of MARSS() on GitHub, you can use `MARSS()` as normal.
#' 
#' @param y Vector of observations n x T.
#' @param model list with 
#'    * R  "diagonal and equal", "unconstrained", "diagonal and unequal"
#'    * tinitx = 1 (must be)
#'    * (if form="dfa") m the number of states (factors). default is 1
#' @param inits list of initial conditions
#' @param miss.value A parameter for backcompatibility. Not used.
#' @param method "TMB" (default), "nlminb.TMB", or "BFGS.TMB". See details.
#' @param form The equation form used in the marssTMB() call. The default is "dfa". 
#' @param fit Whether to fit the model.
#' @param silent Show output when fitting
#' @param control list for the optimization function. See [stats::nlminb()] or [stats::optim()] for the options.
#' @param ... Extra parameters. Not used.
#' 
#' @return The output list from [MARSStmb()]
#' @example inst/examples/MARSS_TMB_example.R
#' @author Eli Holmes
#' @export
MARSS_tmb <- function(y,
                     model = NULL,
                     inits = NULL,
                     miss.value = as.numeric(NA),
                     method = c("TMB", "BFGS.TMB", "nlminb.TMB"),
                     form = c("marxss", "dfa"),
                     fit = TRUE,
                     silent = FALSE,
                     control = NULL,
                     ...) {
  pkg <- "marssTMB"
  form <- match.arg(form)
  method <- match.arg(method)
  # This section is temporary
  fun.opt <- ifelse(method %in% c("TMB", "nlminb.TMB"), "nlminb", "optim")
  if(fun.opt == "optim") optim.method <- strsplit(method, "[.]")[[1]][1]

  # Set up a MARSS() call list and call to get a properly set up MARSS model
  MARSS.call <- list(y = y, inits = inits, model = model, control = control, method = "BFGS", form = form, silent = TRUE, fit = FALSE, ...)
  # x is now a marssMLE object but no par element since not fit
  x <- do.call(MARSS::MARSS, MARSS.call)
  control <- x[["control"]] # set to MARSSoptim defaults
  x[["method"]] <- method

  control[["tmb.silent"]] <- TRUE
  # set up control defaults
  if(fun.opt == "nlminb"){
    allowed.in.opt.fun <- c("eval.max", "iter.max", "trace", "abs.tol", "rel.tol", "x.tol", "xf.tol", "step.min", "step.max", "sing.tol", "scale.init", "diff.g")
    if(is.null(control[["iter.max"]])){
      control[["iter.max"]] = control[["maxit"]]
    }else{ 
      control[["maxit"]] <- control[["iter.max"]]
    }
    if(is.null(control[["eval.max"]])) control[["eval.max"]] = control[["maxit"]]
    for(val in allowed.in.opt.fun) if(is.null(control[[val]])) control[[val]] <- NULL
  }
  if(fun.opt == "optim"){
    allowed.in.opt.fun <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol", "alpha", "beta", "gamma", "REPORT", "warn.1d.NelderMead", "type", "lmm", "factr", "pgtol", "temp", "tmax")
    if(is.null(control$reltol)) control$reltol = 1e-12
    for(val in allowed.in.opt.fun) if(is.null(control[[val]])) control[[val]] <- NULL
  }
  # Temporary until TMB added to the optim methods in {MARSS}
  control[["minit"]] <- 0
  # These are not allowed in the optimization function but keep in control
  allowed.in.control <- c("tmb.silent", "silent", "maxit", "minit")
  control <- control[names(control) %in% c(allowed.in.control, allowed.in.opt.fun)]
  x[["control"]] <- control
  
  # Error check for the DFA model
  model.descrip <- MARSS:::describe.marssMODEL(x[["model"]])

  is.unconstrained <- function(elem) substr(model.descrip[[elem]], 1, 5) == "uncon"
  is.diagonal <- function(elem){
    substr(model.descrip[[elem]], 1, 8) == "diagonal" | model.descrip[[elem]] == "scalar (1 x 1)" }
  is.fixed <- function(elem) substr(model.descrip[[elem]], 1, 5) == "fixed"
  is.identity <- function(elem){
    val <- model.descrip[[elem]]
    substr(val, 1, 8) == "identity" | val == "fixed and all one (1 x 1)"
  }
  is.zero <- function(elem) substr(model.descrip[[elem]], 1, 14) == "fixed and zero"

  # Check that R and Q are allowed
  for(elem in c("R", "Q")){
    ok <- is.diagonal(elem) | is.fixed(elem) | is.unconstrained(elem) | is.identity(elem)
    if(!ok) stop(paste0(pkg, ": ", elem, " must be diagonal, fixed or unconstrained"))
  }

  # Check that B is identity
  for(elem in c("B")){
  ok <- is.identity(elem)
  if(!ok) stop(paste0(pkg, ": ", elem, " must be identity"))
  }
  
  # Check that V0 is fixed
  for(elem in c("V0")){
    ok <- is.fixed(elem) | is.identity(elem)
    if(!ok) stop(paste0(pkg, ": ", elem, " must be fixed"))
  }

  if(fit) return(MARSStmb(x))
  return(x)
}
