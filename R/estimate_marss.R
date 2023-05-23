#' Internal function: MARSS parameter estimation using TMB
#'
#' This model is in the general "marss" vectorized form. The diagonals and offdiagonals of R and Q are split apart like Tim Cline did. In [estimate_marss2()], I use instead the approach in [MARSS::MARSSoptim()] where I just use the [chol()] of these matrices. Technically there are 2 equal solutions since the diagonals appear as the square so -a and a are the same. But I have not observed that this affects the behavior of optim().
#'
#' Minimal error checking is done in this function.
#' Normal calling function is [MARSS::MARSS()] with `method="TMB"`.
#'
#' The main thing this does is
#' * collapse the 3D fixed and free matrices into 2D
#' * separate out the diag and offdiag parameter elements of R and Q
#'
#' Restrictions
#' * V0 fixed (not estimated)
#' * Q and R cannot be time-varying (at the moment)
#'
#' @param MLEobj A properly formatted MARSS model as output by [MARSS()]
#' @param method Normally passed in as MLEobj$method, but allows user to pass in a new method if they want to use MLEobj with another method. Allowed values are "TMB", "nlminb.TMB", "BFGS.TMB".
#' @param opt.control Normally this is passed in as MLEobj$control, but if the MLEobj was set up using a different method, then you will need to set the opt.control options. See details.
#' @param ... not used
#' @details
#' `opt.control` is what is passed to the control argument in [nlminb()] or [optim()]. If you use `MARSS(x, method="TMB")`, this will be set to appropriate defaults which you can see via `MLEobj$control`. But if you call `estimate_marss()` with a MLEobj from a call such as `MARSS(x, method="kem")` (so not a TMB method), you will need to set `opt.control` if you want values different from the base defaults for those functions. Note as a shortcut for `nlminb()`, you can set both `eval.max`, `iter.max` to the same value with `opt.control=list(maxit=1000)`. Note, if you pass in `opt.control`, this will replace all values currently in `MLEobj$control` that are associated with the optimizer function.
#'
#' The defaults set in [MARSS::MARSS()] are
#'
#' * `nlminb`: `eval.max = 5000`, `iter.max = 5000` and `trace = 0`.
#' * `optim`: `maxit = 5000` and `trace = 0`
#'
#' All other controls for the optimization function are left at NULL.

#' @return A list with the objective and optimization objects.
#' * `obj` is the raw output from the [TMB::MakeADFun()] call.
#' * `op` is the raw output from the optimization call (optim or nlminb). Note that the function is minimizing the negative log-likelihood so the sign will be opposite of the log-likelihood reported by MARSS()
#' * `opt.control` is the controls sent to the optimization function.
#' * `method` method used for optimization

#' @author Eli Holmes.
#' @example inst/examples/estimate_marss.R
#' @seealso [MARSS::MARSSoptim()], [MARSS::MARSSkem()]
#' @export
estimate_marss <- function(MLEobj, method = c("TMB", "nlminb_TMB", "BFGS_TMB"), opt.control = NULL, ...) {
  if (!inherits(MLEobj, "marssMLE")) {
    stop("marssTMB::estimate_marss_parameters requires a marssMLE object from the MARSS package.")
  }
  if (missing(method)) {
    if (MLEobj[["method"]] %in% c("TMB", "nlminb_TMB", "BFGS_TMB")) {
      method <- MLEobj[["method"]]
    } else {
      cat("MARSSfit.TMB: The method in the marssMLE object is not TMB. Setting to nlminb_TMB for fitting.\n")
      method <- "nlminb_TMB"
    }
  } else {
    method <- match.arg(method)
  }
  pkg <- "estimate_marss"
  MODELobj <- MLEobj[["marss"]]
  y <- MODELobj[["data"]]
  model.dims <- attr(MODELobj, "model.dims")
  n <- model.dims[["data"]][1]
  TT <- model.dims[["data"]][2]
  m <- model.dims[["x"]][1]

  control <- MLEobj[["control"]]
  tmb.silent <- ifelse(is.null(control[["tmb.silent"]]), TRUE, control[["tmb.silent"]])
  fun.opt <- ifelse(method %in% c("TMB", "nlminb_TMB"), "nlminb", "optim")
  if (fun.opt == "optim") {
    optim.method <- strsplit(method, "[_]")[[1]][1]
  }
  if (fun.opt == "nlminb") {
    if (!is.null(opt.control)) { # user passed in a value
      if (!is.null(opt.control$maxit)) {
        opt.control$eval.max <- opt.control$iter.max <- opt.control$maxit
        opt.control$maxit <- NULL
      }
      opt.control <- opt.control
    } else {
      opt.control <- control
    }
    allowed.in.opt.control <- c("eval.max", "iter.max", "trace", "abs.tol", "rel.tol", "x.tol", "xf.tol", "step.min", "step.max", "sing.tol", "scale.init", "diff.g")
    opt.control <- opt.control[names(opt.control) %in% allowed.in.opt.control]
  }
  if (fun.opt == "optim") {
    opt.control <- ifelse(!is.null(opt.control), opt.control, control)
    allowed.in.opt.control <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol", "alpha", "beta", "gamma", "REPORT", "warn.1d.NelderMead", "type", "lmm", "factr", "pgtol", "temp", "tmax")
    opt.control <- opt.control[names(opt.control) %in% allowed.in.opt.control]
  }

  # Set up the initial matrices
  eleminits <- list()
  model.elem <- attr(MODELobj, "par.names")
  for (elem in model.elem) {
    eleminits[[elem]] <- coef(MLEobj, type = "matrix", what = "start")[[elem]]
  }
  # Check that no 0s on diagonal of V0 unless V0 is all zero
  # Note V0 is fixed by definition
  V0_is_zero <- all(eleminits[["V0"]] == 0)
  if (!V0_is_zero && any(diag(eleminits[["V0"]]) == 0)) {
    stop(paste0(pkg, ": V0 can only have 0s on the diagonal if it is all zero"))
  }
  # No zeros on diagonal of Q or R by definition
  for (elem in c("R", "Q")) {
    if (is.na(dim(eleminits[[elem]])[3])) {
      dim(eleminits[[elem]]) <- c(dim(eleminits[[elem]]), 1)
    }
    bad <- apply(eleminits[[elem]], 3, function(x) {
      any(diag(x) == 0)
    }) |>
      unlist() |>
      any()
    if (bad) stop(paste0(pkg, ": No zeros allowed on the diagonal of ", elem, "."))
  }

  free <- MODELobj$free
  fixed <- MODELobj$fixed
  pars <- MLEobj$start
  par_dims <- attr(MODELobj, "model.dims")
  # Helper code to find the diagonals and off-diagonal pars in the vector
  # Need to fix if offdiags of start are not 0
  # Why is this so horrible?? Because the free matrix is a general linear
  # constraint matrix and I don't know what it looks like. It have 1 or 0 for
  # RowSums and ColSums but the ordering depends on how the user constructed
  # their Q and R constraints. This painful code makes no assumptions and uses
  # the free matrix to back calculate which par are assoc with diagonal or
  # off diagonal elements.
  for (elem in c("Q", "R")) {
    dname <- paste0(elem, "diag")
    oname <- paste0(elem, "offdiag")
    mn <- sqrt(dim(free[[elem]])[1])
    d <- seq(1, mn * mn, mn + 1) # diagonal rows
    ut <- matrix(1:(mn * mn), mn, mn)
    ut <- ut[upper.tri(ut)]
    od <- (1:(mn * mn))
    od <- od[!(od %in% d)] # offdiagonal rows
    free[[dname]] <- free[[elem]][d, , , drop = FALSE]
    fixed[[dname]] <- fixed[[elem]][d, , , drop = FALSE]
    free[[oname]] <- free[[elem]]
    free[[oname]] <- free[[elem]][ut, , , drop = FALSE] # only upper tri
    fixed[[oname]] <- fixed[[elem]][ut, , , drop = FALSE]
    # Drop cols and rows associated with values not on d or o
    pars[[dname]] <- pars[[elem]][colSums(free[[dname]]) > 0, , drop = FALSE]
    pars[[oname]] <- pars[[elem]][colSums(free[[oname]]) > 0, , drop = FALSE]
    free[[dname]] <- free[[dname]][, colSums(free[[dname]]) > 0, , drop = FALSE]
    free[[oname]] <- free[[oname]][, colSums(free[[oname]]) > 0, , drop = FALSE]
    # Works since no linear constraints or mixing of fixed/free on
    # diagonal allowed
    pars[[dname]] <- log(sqrt(pars[[dname]]))
    fixed[[dname]][fixed[[dname]] != 0] <- log(sqrt(fixed[[dname]][fixed[[dname]] != 0]))
    # Chols of corr mat; these are either in fixed or in par
    if (mn > 1) { # there are off-diagonals
      no <- mn * (mn - 1) / 2 # number of upper tri off diagonals
      chols <- apply(eleminits[[elem]], 3, function(x) {
        tmp <- chol(diag(1 / sqrt(diag(x))) %*% x %*% diag(1 / sqrt(diag(x))))
        tmp[upper.tri(tmp)]
      })
      # guard against R dropping the dims
      if (is.null(dim(chols))) chols <- matrix(chols, nrow = 1)
      parind <- apply(free[[oname]], 3, function(x) {
        rowSums(x) > 0
      })
      fixed[[oname]] <- array(chols * !parind, dim = c(no, 1, dim(chols)[2]))
      tmp <- 0 * pars[[oname]]
      if (dim(free[[oname]])[2] != 0) {
        for (t in 1:dim(chols)[2]) tmp <- tmp + t(free[[oname]][, , t]) %*% (chols * parind)[, t]
        pars[[oname]] <- tmp
      }
    }
    par_dims[[dname]] <- c(mn, 1, 1)
    par_dims[[oname]] <- c(mn * (mn - 1) / 2, 1, 1)
  }
  free <- free[!(names(free) %in% c("Q", "R"))]
  fixed <- fixed[!(names(fixed) %in% c("Q", "R"))]
  pars <- pars[!(names(pars) %in% c("Q", "R"))]

  tfixed <- unlist(lapply(fixed, function(x) {
    dim(x)[3]
  }))
  tfree <- unlist(lapply(free, function(x) {
    dim(x)[3]
  }))
  # Make time-varying matrices in 2D vec form d1*d2*npar x T for free and d1*d2 x T for fixed
  free <- lapply(free, function(x) {
    new <- matrix(x, max(dim(x)[1], 1) * max(dim(x)[2], 1), dim(x)[3])
    new[is.na(new)] <- 0
    new
  })
  fixed <- lapply(fixed, function(x) {
    matrix(x, dim(x)[1], dim(x)[3])
  })
  pars0 <- lapply(pars, function(x) {
    new <- matrix(x, max(dim(x)[1], 1), dim(x)[2])
    new[is.na(new)] <- 0
    new
  })
  par_dims <- par_dims[names(fixed)]
  par_dims <- lapply(par_dims, as.integer)
  numpar <- unlist(lapply(pars, nrow))
  # a vector of the parameters
  pars <- unlist(lapply(pars, function(x) {
    colnames(x) <- NULL
    x[, 1]
  }))

  # Creates the input data list
  # For now, we will assume V0 is a fixed (diagonal) matrix
  data <- list(
    model = "marss",
    Y = y,
    V0_is_zero = as.numeric(V0_is_zero),
    tinitx = MODELobj[["tinitx"]],
    tfixed = tfixed,
    tfree = tfree,
    numpar = numpar,
    free = free,
    fixed = fixed,
    par_dims = par_dims
  )

  # Creates the list of initial (start) values of parameter list
  parameters <- list(
    X = matrix(0, ncol = TT, nrow = m), # states
    pars = pars
  )

  # The map (mask) that indicates what parameters to not estimate, however the
  # par list only includes estimated parameters. Thus we only need to worry
  # about the case where X(0) is fixed via x0
  maplist <- list()
  if (MODELobj[["tinitx"]] == 1 & V0_is_zero) {
    mat <- matrix(1:(m * TT), m, TT)
    mat[, 1] <- NA
    maplist$X <- mat |>
      unlist() |>
      as.factor()
  }

  if (is.null(MLEobj[["control"]][["tmb.silent"]])) MLEobj[["control"]][["tmb.silent"]] <- TRUE
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
    # Only show output if trace = 1
    if (opt.control$trace > 1) opt.control$trace <- 0
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
  names(opt1$par) <- names(pars)

  return(list(obj = obj1, opt = opt1, fun.opt = fun.opt, opt.control = opt.control, method = method))
}
