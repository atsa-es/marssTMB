#' Internal function: MARSS parameter estimation using TMB
#' 
#' Minimal error checking is done in this function.
#' Normal calling function is [MARSS::MARSS()] with `method="TMB"`. 
#' 
#' The main thing this does is 
#' * collapse the 3D fixed and free matrices into 2D
#' * separate out the diag and offdiag parameter elements of R and Q
#'
#' Restrictions
#' * No time-varying parameters
#' * V0 fixed (not estimated)
#' * No B estimation
#'
#' @param MLEobj A properly formatted MARSS model as output by [MARSS()]

#' @return A list with the objective and optimization objects.
#' * `obj` is the raw output from the [TMB::MakeADFun()] call.
#' * `op` is the raw output from the optimization call (optim or nlminb)

#' @author Eli Holmes. This function is inspired by dfaTMB.R written by Tim Cline while a graduate student in the Fish 507 Time Series Analysis course.
#' @seealso [MARSS::MARSSoptim()], [MARSS::MARSSkem()]
#' @export
estimate_marss <- function(MLEobj) {
  if(!inherits(MLEobj, "marssMLE"))
    stop("marssTMB::estimate_marss_parameters requires a marssMLE object from the MARSS package.")
  pkg <- "estimate_marss"
  MODELobj <- MLEobj[["model"]]
  y <- MODELobj[["data"]]
  model.dims <- attr(MODELobj, "model.dims")
  n <- model.dims[["data"]][1]
  TT <- model.dims[["data"]][2]
  m <- model.dims[["x"]][1]

  control <- MLEobj[["control"]]
  tmb.silent <- control[["tmb.silent"]]
  fun.opt <- ifelse(MLEobj[["method"]] %in% c("TMB", "nlminb.TMB"), "nlminb", "optim")
  if(fun.opt == "optim") optim.method <- strsplit(MLEobj[["method"]], "[.]")[[1]][1]
  if(fun.opt == "nlminb")
    opt.control <- control[!(names(control) %in% c("fun.opt", "optim.method", "tmb.silent", "silent", "maxit", "minit"))]
  if(fun.opt == "optim")
    opt.control <- control[!(names(control) %in% c("fun.opt", "optim.method", "tmb.silent", "silent", "minit"))]

  # Set up the initial matrices
  eleminits <- list()
  model.elem <- attr(MODELobj, "par.names")
  for (elem in model.elem[!(model.elem %in% c("c", "d"))]) {
    eleminits[[elem]] <- coef(MLEobj, type = "matrix", what = "start")[[elem]]
  }
  # Check that no 0s on diagonal of V0 unless V0 is all zero
  # Note V0 is fixed by definition
  V0_is_zero <- all(eleminits[["V0"]]==0)
  if(!V0_is_zero && any(diag(eleminits[["V0"]])==0))
    stop(paste0(pkg, ": V0 can only have 0s on the diagonal if it is all zero"))
  # No zeros on diagonal of Q or R by definition
  for(elem in c("R", "Q")){
  if(any(diag(eleminits[[elem]])==0))
    stop(paste0(pkg, ": No zeros allowed on the diagonal of ", elem, "."))
  }
  
  free <- MLEobj$marss$free
  fixed <- MLEobj$marss$fixed
  pars <- MLEobj$start
  par_dims <- attr(MLEobj$marss, "model.dims")
  # Helper code to find the diagonals and off-diagonal pars in the vector
  # Need to fix if offdiags of start are not 0!
  for(elem in c("Q", "R")){
    dname <- paste0(elem, "diag")
    oname <- paste0(elem, "offdiag")
    mn <- sqrt(dim(free[[elem]])[1])
    d <- seq(1, mn*mn, mn+1) # diagonal rows
    ut <- matrix(1:(mn*mn), mn, mn); ut <- ut[upper.tri(ut)] # offdiagonal rows
    od <-  (1:(mn*mn)); od <- od[!(od %in% d)] # offdiagonal rows
    free[[dname]] <- free[[elem]][d,,,drop=FALSE]
    fixed[[dname]] <- fixed[[elem]][d,,,drop=FALSE]
    free[[oname]] <- free[[elem]]
    free[[oname]][ut,,] <- 0
    free[[oname]] <- free[[elem]][od,,,drop=FALSE]
    fixed[[oname]] <- fixed[[elem]][od,,,drop=FALSE]
    pars[[dname]] <- pars[[elem]][colSums(free[[dname]])>0,,drop=FALSE]
    pars[[oname]] <- pars[[elem]][colSums(free[[oname]])>0,,drop=FALSE]
    pars[[dname]] <- log(sqrt(pars[[dname]]))
    par_dims[[dname]] <- par_dims[[elem]]
    par_dims[[oname]] <- par_dims[[elem]]
  }
  free <- free[!(names(free) %in% c("Q", "R"))]
  fixed <- fixed[!(names(fixed) %in% c("Q", "R"))]
  pars <- pars[!(names(pars) %in% c("Q", "R"))]

  tfixed <- unlist(lapply(fixed, function(x){dim(x)[3]}))
  tfree <- unlist(lapply(free, function(x){dim(x)[3]}))
  free = lapply(free, function(x){new <- matrix(x, dim(x)[1]*max(dim(x)[2],1), dim(x)[3]); new[is.na(new)] <- 0; new})
  fixed = lapply(fixed, function(x){matrix(x, dim(x)[1], dim(x)[3])})
  pars0 = lapply(pars, function(x){new <- matrix(x, max(dim(x)[1],1), dim(x)[2]); new[is.na(new)]<- 0; new})
  par_dims = par_dims[names(fixed)]
  par_dims = lapply(par_dims, as.integer)
  npar <- unlist(lapply(pars, nrow))
  # a vector of the parameters
  pars <- unlist(lapply(pars, function(x){x[,1]}))
  
  # Creates the input data list
  # For now, we will assume V0 is a fixed (diagonal) matrix
  data <- list(
    model = "marxss2",
    Y = y,
    V0_is_zero = as.numeric(V0_is_zero),
    tinitx = MODELobj[["tinitx"]],
    free = free,
    fixed = fixed,
    par_dims = par_dims,
    npar = npar,
    tfixed = tfixed,
    tfree = tfree
  )

  # Note x0 and V0 are fixed (stochastic prior) for DFA
  # But x0 might be estimated in the future so is here as well

  # Creates the list of initial (start) values of parameter list
  parameters <- list(
    X = matrix(0, ncol = TT, nrow = m), # states
    pars = pars
  )

  # Create the map (mask) that indicates what parameters to not estimate
  # the states are treated as a parameter but are all estimated
  # so do not appear here
  maplist <- list()
  if(MODELobj[["tinitx"]]==1 & V0_is_zero){
    mat <- matrix(1:(m*TT), m, TT)
    mat[,1] <- NA
    maplist$X <- mat |> unlist() |> as.factor()
  }

  MLEobj[["control"]][["tmb.silent"]] <- TRUE
  
  # Creates the model object and runs the optimization
  obj1 <- TMB::MakeADFun(
    data, parameters,
    random = "X",
    DLL = "marssTMB_TMBExports",
    silent = MLEobj[["control"]][["tmb.silent"]], 
    map = maplist
  )

  # Optimization
  # remove any NULLs
  opt.control <- opt.control[unlist(lapply(opt.control, function(x) !is.null(x)))]
  if (fun.opt == "nlminb") {
    opt1 <- stats::nlminb(obj1$par, obj1$fn, obj1$gr, control = opt.control)
    # add output also found in optim output
    opt1$value <- opt1$objective
    opt1$counts <- opt1$evaluations
  }
  if (fun.opt == "optim") {
    opt1 <- stats::optim(obj1$par, obj1$fn, gr = obj1$gr, control = opt.control, method = optim.method)
    opt1$objective <- opt1$value
    opt1$iterations <- opt1$counts[2]
  }

  return(list(obj = obj1, opt = opt1))
}
