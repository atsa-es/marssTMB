#' Internal function: helper function to make marssMLE obj
#' 
#' Takes output from [MARSStmb()] and makes into [MARSS::marssMLE] object.
#' 
#' @param x list as output from [MARSStmb()]: has obj.function and opt.output
#' @return A [MARSS::marssMLE] object
#' @export
to_marssMLE <- function(x){
  obj1 <- x$obj.function
  opt.output <- x$opt.output
  MLEobj <- x$MLEobj
  control <- MLEobj$control
  fun.opt <-paste0(c(control$fun.opt, control$optim.method), collapse=" method=")
  
MLEobj.return <- MLEobj
MLEobj.return$iter.record <- opt.output$message
MLEobj.return$start <- MLEobj$start 
MLEobj.return$convergence <- opt.output$convergence
if (opt.output$convergence %in% c(1, 0)) {
  funinfo <- paste0("Function ", fun.opt, "used for optimization and TMB for likelihood calculation.\n")
  if ((!control$silent || control$silent == 2) && opt.output$convergence == 0) cat(paste0("Success! Converged in ", opt.output$iterations, " iterations.\n", funinfo))
  if ((!control$silent || control$silent == 2) && opt.output$convergence == 1) cat(paste0("Warning! Max iterations of ", control$maxit, " reached before convergence.\n", funinfo))
  
  browser()
  tmp.MLEobj <- MARSSvectorizeparam(tmp.MLEobj, opt.output$par)
  # par has the fixed and estimated values using t chol of Q and R
  
  # back transform Q, R and V0 if needed from chol form to usual form
  for (elem in c("Q", "R", "V0")) { # this works because by def fixed and free blocks of var-cov mats are independent
    if (!is.fixed(MODELobj$free[[elem]])) # get a new par if needed
    {
      d <- sub3D(tmp.MLEobj$marss$free[[elem]], t = 1) # this will be the one with the upper tri zero-ed out but ok since symmetric
      par.dim <- par.dims[[elem]][1:2]
      L <- unvec(tmp.MLEobj$marss$free[[elem]][, , 1] %*% tmp.MLEobj$par[[elem]], dim = par.dim) # this by def will have 0 row/col at the fixed values
      the.par <- tcrossprod(L) # L%*%t(L)
      tmp.MLEobj$par[[elem]] <- solve(crossprod(d)) %*% t(d) %*% vec(the.par)
    }
  } # end for
  
  pars <- MARSSvectorizeparam(tmp.MLEobj) # now the pars values have been adjusted back to normal scaling
  # now put the estimated values back into the original MLEobj; fixed and free matrices as in original
  MLEobj.return <- MARSSvectorizeparam(MLEobj.return, pars)
  kf.out <- try(MARSSkf(MLEobj.return), silent = TRUE)
  
  if (inherits(kf.out, "try-error")) {
    MLEobj.return$numIter <- optim.output$counts[1]
    MLEobj.return$logLik <- -1 * optim.output$value
    MLEobj.return$errors <- c(paste0("\nWARNING: optim() successfully fit the model but ", kf.function, " returned an error with the fitted model. Try MARSSinfo('optimerror54') for insight.", sep = ""), "\nError: ", kf.out[1])
    MLEobj.return$convergence <- 54
    MLEobj.return <- MARSSaic(MLEobj.return)
    kf.out <- NULL
  }
} else {
  if (optim.output$convergence == 10) optim.output$message <- c("degeneracy of the Nelder-Mead simplex\n", paste("Function ", kf.function, " used for likelihood calculation.\n", sep = ""), optim.output$message)
  optim.output$counts <- NULL
  if (!control$silent) cat("MARSSoptim() stopped with errors. No parameter estimates returned.\n")
  if (control$silent == 2) cat("MARSSoptim() stopped with errors. No parameter estimates returned. See $errors in output for details.\n")
  
  MLEobj.return$par <- NULL
  MLEobj.return$errors <- optim.output$message
  kf.out <- NULL
}

if (!is.null(kf.out)) {
  if (control$trace > 0) MLEobj.return$kf <- kf.out
  MLEobj.return$states <- kf.out$xtT
  MLEobj.return$numIter <- optim.output$counts[1]
  MLEobj.return$logLik <- kf.out$logLik
}
MLEobj.return$method <- MLEobj$method

## Add AIC and AICc to the object
if (!is.null(kf.out)) MLEobj.return <- MARSSaic(MLEobj.return)

return(MLEobj.return)
}
