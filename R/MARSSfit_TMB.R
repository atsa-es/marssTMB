#' MARSS parameter estimation using TMB
#'
#' This takes a [MARSS::marssMLE] object (fitted or not) and estimates 
#' parameters. Note this is the TMB method for the [MARSS::MARSSfit] generic. 
#' The typical use would be to call as `MARSS(data, method="TMB")`.
#'
#' Restrictions
#' * V0 fixed (not estimated)
#' * Q and R cannot be time-varying (at the moment)
#'
#' @param x A properly formatted [MARSS::marssMLE] object ready for fitting.
#' @param fun A debugging option to switch between estimate_marss() and estimate_marss2()
#' @param ... not used
#'
#' @details
#' `opt.control` is what is passed to the control argument in [nlminb()] or [optim()]. If you use `fit <- MARSS(data, method="TMB")`, this will be set to appropriate defaults which you can see via `fit$control`. But if you call `estimate_marss()` with a marssMLE object from a call such as `MARSS(data, method="kem")` (so not a TMB method), you will need to set `opt.control` if you want values different from the base defaults for those functions. Note as a shortcut for `nlminb()`, you can set both `eval.max`, `iter.max` to the same value with `opt.control=list(maxit=1000)`. Note, if you pass in `opt.control`, this will replace all values currently in `fit$control` that are associated with the optimizer function.
#'
#' The defaults set in [MARSS::MARSS()] are
#'
#' * `nlminb`: `eval.max = 5000`, `iter.max = 5000` and `trace = 0`.
#' * `optim`: `maxit = 5000` and `trace = 0`
#'
#' All other controls for the optimization function are left at NULL. You can
#' set other controls in the call `MARSS(..., control=list(...))`.

#' @return The [MARSS::marssMLE] object which was passed in, with additional components:
#'
#' * `method`: From the call or argument method if user passed that in.
#' * `kf`: Kalman filter output.
#' * `iter.record`: If \code{x$control$trace = TRUE}, then this is the \code{$message} value from [stats::optim()] or [stats::nlminb()] plus the output from the [TMB::MakeADFun()] call and the output from the optimization function.
#' * `numIter`: Number of iterations needed for convergence.
#' * `convergence`: Did estimation converge successfully? 
#'   - `convergence=0`: Converged in less than \code{x$control$maxit} iterations and no evidence of degenerate solution. 
#'   - `convergence=3`: No convergence diagnostics were computed because all parameters were fixed thus no fitting required.
#'   - `convergence=-1`: No convergence diagnostics were computed because the MLE object was not fit (called with fit=FALSE). This is not a convergence error just information. There is not par element so no functions can be run with the object.
#'   - `convergence=1`: Maximum number of iterations \code{x$control$maxit} was reached before convergence.
#'   - For other convergence errors, see[stats::optim()] or [stats::nlminb()].
#' * `logLik`: Log-likelihood.
#' * `states`: State estimates from the Kalman smoother.
#' * `states.se`: Confidence intervals based on state standard errors, see caption of Fig 6.3 (p. 337) in Shumway & Stoffer (2006).
#' * `errors`: Any error messages.
#'
#' @author Eli Holmes. 
#' @example inst/examples/MARSS_TMB_example.R
#' @seealso [MARSS::MARSSoptim()], [MARSS::MARSSkem()]
#' @export
MARSSfit.TMB <- function(x, fun=1, ...) {
  
  if(fun==1) out <- marssTMB::estimate_marss(x)
  if(fun==2) out <- marssTMB::estimate_marss2(x)
  obj1 <- out$obj
  opt1 <- out$opt
  
  parvec <- opt1$par
  nonvar_elems <- c("B","U", "Z", "A", "G", "H", "L", "x0", "V0")
  keep <- sapply(names(parvec), function(x){strsplit(x, "[.]")[[1]][1]}) %in% nonvar_elems
  parvec <- parvec[keep]
  
  parlist <- list()
  for (elem in c("R", "Q")) {
    if (dim(x$model$free[[elem]])[2]!=0) { # get a new par if needed
      val <- paste0("FullCovMat", elem)
      the.par <- obj1$report()[[val]]
      d <- sub3D(x$model$free[[elem]], t = 1)
      # A bit of a hack but I want to allow any varcov contraints (d mat)
      # Also ensures that the par names are in the right order;
      # They might not be since TMB code split out the diag separate from offdiag
      tmp <- solve(crossprod(d)) %*% t(d) %*% vec(the.par)
      tmpnames <- paste0(elem, ".", rownames(tmp))
      attr(tmp, "dim") <- NULL
      attr(tmp, "names") <- tmpnames
      parvec <- c(parvec, tmp)
    } 
  }
  ord <- names(coef(x, what="start", type="vector", form="marss"))
  parvec <- parvec[ord]
  x <- MARSS::MARSSvectorizeparam(x, parvec = parvec)
  
  x$iter.record <- list(message = opt1$message, opt.output = opt1)
  if(x$control$trace > 0)
    x$iter.record <-c(x$iter.record, list(obj.function = obj1))
  x$convergence <- opt1$convergence
  x$numIter <- opt1$iterations
  x$logLik <- -1*x$iter.record$opt.output$objective

  return(x)
}