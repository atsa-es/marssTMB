#' Internal function: Parameter estimation using TMB
#'
#' This function was the first working version and used a 
#' non-vectorized marxss form (aka the form you would write out on paper).
#' It is replaced by [estimate_marss()] but is kept here because the 
#' code might be easier to follow.
#' 
#' Minimal error checking is done in this function.
#' Normal calling function is [MARSS_tmb()] which in
#' turn calls this function. Note when MARSS is updated, the
#' normal calling function will be [MARSS::MARSS()]. Restrictions
#'
#' * No time-varying parameters
#' * Currently only DFA models are coded up
#' * x0 and V0 fixed (stochastic prior)
#' * Q is fixed (not estimated)
#'
#' @details
#' This function returns [MARSS::marssMLE] object. For the `iter.record` element
#' of the object, the following are returned as a list:
#'
#' * `obj.function` is the raw output from the [TMB::MakeADFun()] call.
#' * `opt.output` is the raw output from the optimization call (optim or nlminb)
#'
#' @param MLEobj A properly formatted MARSS model as output by [MARSS_tmb()]

#' @return A [MARSS::marssMLE] object
#' @example inst/examples/estimate_marxss.R
#' @author Eli Holmes. This function is inspired by dfaTMB.R written by Tim Cline while a graduate student in the Fish 507 Time Series Analysis course.
#' @seealso [MARSS::MARSSoptim()], [MARSS::MARSSkem()]
#' @export
estimate_marxss <- function(MLEobj) {
  pkg <- "estimate_marxss"
  MODELobj <- MLEobj[["model"]]
  y <- MODELobj[["data"]]
  model.dims <- attr(MODELobj, "model.dims")
  n <- model.dims[["data"]][1]
  TT <- model.dims[["data"]][2]
  m <- model.dims[["x"]][1]
  d_Covars <- MODELobj[["fixed"]][["d"]]
  c_Covars <- MODELobj[["fixed"]][["c"]]
  
  control <- MLEobj[["control"]]
  tmb.silent <- control[["tmb.silent"]]
  fun.opt <- ifelse(MLEobj[["method"]] %in% c("TMB", "nlminb.TMB"), "nlminb", "optim")
  if(fun.opt == "optim") optim.method <- strsplit(MLEobj[["method"]], "[.]")[[1]][1]
  if(fun.opt == "nlminb")
    opt.control <- control[!(names(control) %in% c("fun.opt", "optim.method", "tmb.silent", "silent", "maxit", "minit"))]
  if(fun.opt == "optim")
    opt.control <- control[!(names(control) %in% c("fun.opt", "optim.method", "tmb.silent", "silent", "minit"))]

  # Expand out to full covariate matrix
  d_Covars <- matrix(d_Covars[,1,], nrow = nrow(d_Covars))
  c_Covars <- matrix(c_Covars[,1,], nrow = nrow(c_Covars))
  
  # Set up the initial matrices
  eleminits <- list()
  model.elem <- attr(MODELobj, "par.names")
  for (elem in model.elem[!(model.elem %in% c("c", "d"))]) {
    eleminits[[elem]] <- coef(MLEobj, type = "matrix", what = "start")[[elem]]
  }
  # Check that no 0s on diagonal of V0 unless V0 is all zero
  # Note V0 is fixed by definition
  V0_is_zero <- all(eleminits[["V0"]]==0)
  if(!V0_is_zero && any(diag(eleminits[[elem]])==0))
    stop(paste0(pkg, ": V0 can only have 0s on the diagonal if it is all zero"))
  # No zeros on diagonal of Q or Rdefinition
  for(elem in c("R", "Q")){
  if(any(diag(eleminits[[elem]])==0))
    stop(paste0(pkg, ": No zeros allowed on the diagonal of ", elem, "."))
  }
  
  # Set up the maps
  elemmaps <- list()
  for (elem in model.elem) {
    elemmaps[[elem]] <- create.elem.maps(MLEobj, elem = elem)[["map"]]
  }
  # maps for var-cov matrices have diagonal separate from off-diagonal
  for (elem in c("R", "Q", "V0")) {
    elemmaps[[elem]] <- list(
      diag = create.varcov.maps(MLEobj, elem = elem)[["map.diag"]],
      offdiag = create.varcov.maps(MLEobj, elem = elem)[["map.offdiag"]]
    )
  }

  # Creates the input data list
  # For now, we will assume V0 is a fixed (diagonal) matrix
  data <- list(
    model = "marxss",
    Y = y,
    d_Covar = d_Covars,
    c_Covar = c_Covars,
    has_c_covars = as.numeric(ncol(c_Covars) != 1),
    has_d_covars = as.numeric(ncol(d_Covars) != 1),
    V0_is_zero = as.numeric(V0_is_zero),
    tinitx = MODELobj[["tinitx"]]
  )

  # Note x0 and V0 are fixed (stochastic prior) for DFA
  # But x0 might be estimated in the future so is here as well

  # Creates the list of initial (start) values of parameter list
  R <- eleminits[["R"]]
  sdR <- sqrt(diag(R))
  corrR <- diag(1/sdR, n) %*% R %*% diag(1/sdR, n) # correlation matrix
  Q <- eleminits[["Q"]]
  sdQ <- sqrt(diag(Q))
  corrQ <- diag(1/sdQ, m) %*% Q %*% diag(1/sdQ, m) # correlation matrix
  parameters <- list(
    X = matrix(0, ncol = TT, nrow = m), # states
    x0 = eleminits[["x0"]],
    V0 = eleminits[["V0"]],
    logsdQ = log(sdQ), # log of sqrt of diagonal of Q
    cholCorrQ = chol(corrQ)[upper.tri(Q)], # off-diagonal of chol of corr Q
    U = eleminits[["U"]],
    C = eleminits[["C"]],
    Z = eleminits[["Z"]],
    A = eleminits[["A"]],
    D = eleminits[["D"]],
    logsdR = log(sdR), # log of sqrt of diagonal of R
    cholCorrR = chol(corrR)[upper.tri(R)] # off-diagonal of chol of corr R
  )

  # Create the map (mask) that indicates what parameters to not estimate
  # the states are treated as a parameter but are all estimated
  # so do not appear here
  maplist <- list(
    x0 = elemmaps[["x0"]],
    V0 = factor(matrix(NA, nrow = m, ncol = m)),
    logsdQ = elemmaps[["Q"]][["diag"]],
    cholCorrQ = elemmaps[["Q"]][["offdiag"]],
    U = elemmaps[["U"]],
    C = elemmaps[["C"]],
    Z = elemmaps[["Z"]],
    A = elemmaps[["A"]],
    D = elemmaps[["D"]],
    logsdR = elemmaps[["R"]][["diag"]],
    cholCorrR = elemmaps[["R"]][["offdiag"]]
  )
  if(MODELobj[["tinitx"]]==1 & V0_is_zero){
    mat <- matrix(1:(m*TT), m, TT)
    mat[,1] <- NA
    maplist$X <- mat |> unlist() |> as.factor()
  }

  # Creates the model object and runs the optimization
  obj1 <- TMB::MakeADFun(
    data, parameters,
    random = "X",
    DLL = "marssTMB_TMBExports",
    silent = MLEobj[["control"]][["tmb.silent"]], map = maplist
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

  # Add names to the par output; par vec needs to be in this order
  # In MARSS(), Dd and Cc are within the A and U matrices in marss form
  # The par element of a MLEobj is in this marss form
  # {MARSS} has a helper function to convert from marxss (with Dd and Cc) to marss 
  parlist <- list()
  # These include D, d, C, c
  model.elem <- attr(MODELobj, "par.names")
  for (elem in model.elem[!(model.elem %in% c("R", "Q", "V0"))]) {
    if (!MARSS:::is.fixed(MODELobj$free[[elem]])) {
      parlist[[elem]] <- matrix(opt1$par[names(opt1$par) == elem], ncol = 1)
      parnames <- paste0(elem, ".", levels(obj1$env$map[[elem]]))
      names(opt1$par)[names(opt1$par) == elem] <- parnames
      rownames(parlist[[elem]]) <- levels(obj1$env$map[[elem]])
    } else {
      parlist[[elem]] <- matrix(0, nrow = 0, ncol = 1)
    }
  }
  for (elem in c("R", "Q", "V0")) {
    if (!MARSS:::is.fixed(MODELobj$free[[elem]])) { # get a new par if needed
      val <- paste0("FullCovMat", elem)
      the.par <- obj1$report()[[val]]
      d <- MARSS:::sub3D(MLEobj$model$free[[elem]], t = 1)
      # A bit of a hack but I want to allow any varcov contraints (d mat)
      # Also ensures that the par names are in the right order;
      # They might not be since TMB code split out the diag separate from offdiag
      parlist[[elem]] <- solve(crossprod(d)) %*% t(d) %*% MARSS:::vec(the.par)
    } else {
      parlist[[elem]] <- matrix(0, nrow = 0, ncol = 1)
    }
  }
  # par is in marxss form with D, d, C, c
  MLEobj[["par"]] <- parlist[model.elem]
  # We need this in marss form; only.par means the marss element is already ok
  MLEobj <- MARSS:::marxss_to_marss(MLEobj, only.par = TRUE)
  
  MLEobj$iter.record <-
    list(obj.function = obj1, opt.output = opt1)
  MLEobj$start <- MLEobj$start
  MLEobj$convergence <- opt1$convergence

  if (fun.opt == "optim" && opt1$convergence > 1) {
    if (!control$silent) cat(paste0(pkg, "() stopped with errors. No parameter estimates returned.\n"))
    MLEobj$par <- MLEobj$kf <- MLEobj$logLik <- NULL
    MLEobj$errors <- opt1$message
    return(MLEobj)
  }
  
  MLEobj$states <- obj1$env$parList()$X
  rownames(MLEobj$states) <- attr(MLEobj$model, "X.names")
  MLEobj$numIter <- opt1$iterations
  MLEobj$logLik <- -1*opt1$objective
  ## Add AIC and AICc to the object
  MLEobj <- MARSS::MARSSaic(MLEobj)

  funinfo <- paste0("Function ", fun.opt, " used for optimization and TMB for likelihood calculation.\n")
  if ((!control$silent || control$silent == 2) && opt1$convergence == 0) cat(paste0("Success! Converged in ", opt1$iterations, " iterations.\n", funinfo))
  if ((!control$silent || control$silent == 2) && opt1$convergence == 1) cat(paste0("Warning! Max iterations of ", control$maxit, " reached before convergence.\n", funinfo))

  return(MLEobj)
}
