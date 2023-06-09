% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MARSSfit_TMB.R
\name{MARSSfit.TMB}
\alias{MARSSfit.TMB}
\title{MARSS parameter estimation using TMB}
\usage{
\method{MARSSfit}{TMB}(x, fun = 1, ...)
}
\arguments{
\item{x}{A properly formatted \link[MARSS:marssMLE-class]{MARSS::marssMLE} object ready for fitting.}

\item{fun}{A debugging option to switch between estimate_marss() and estimate_marss2()}

\item{...}{not used}
}
\value{
The \link[MARSS:marssMLE-class]{MARSS::marssMLE} object which was passed in, with additional components:
\itemize{
\item \code{method}: From the call or argument method if user passed that in.
\item \code{kf}: Kalman filter output.
\item \code{iter.record}: If \code{x$control$trace = TRUE}, then this is the \code{$message} value from \code{\link[stats:optim]{stats::optim()}} or \code{\link[stats:nlminb]{stats::nlminb()}} plus the output from the \code{\link[TMB:MakeADFun]{TMB::MakeADFun()}} call and the output from the optimization function.
\item \code{numIter}: Number of iterations needed for convergence.
\item \code{convergence}: Did estimation converge successfully?
\itemize{
\item \code{convergence=0}: Converged in less than \code{x$control$maxit} iterations and no evidence of degenerate solution.
\item \code{convergence=3}: No convergence diagnostics were computed because all parameters were fixed thus no fitting required.
\item \code{convergence=-1}: No convergence diagnostics were computed because the MLE object was not fit (called with fit=FALSE). This is not a convergence error just information. There is not par element so no functions can be run with the object.
\item \code{convergence=1}: Maximum number of iterations \code{x$control$maxit} was reached before convergence.
\item For other convergence errors, see\code{\link[stats:optim]{stats::optim()}} or \code{\link[stats:nlminb]{stats::nlminb()}}.
}
\item \code{logLik}: Log-likelihood.
\item \code{states}: State estimates from the Kalman smoother.
\item \code{states.se}: Confidence intervals based on state standard errors, see caption of Fig 6.3 (p. 337) in Shumway & Stoffer (2006).
\item \code{errors}: Any error messages.
}
}
\description{
This takes a \link[MARSS:marssMLE-class]{MARSS::marssMLE} object (fitted or not) and estimates
parameters. Note this is the TMB method for the \link[MARSS:MARSSfit]{MARSS::MARSSfit} generic.
The typical use would be to call as \code{MARSS(data, method="TMB")}.
}
\details{
Restrictions
\itemize{
\item V0 fixed (not estimated)
\item Q and R cannot be time-varying (at the moment)
}

\code{opt.control} is what is passed to the control argument in \code{\link[=nlminb]{nlminb()}} or \code{\link[=optim]{optim()}}. If you use \code{fit <- MARSS(data, method="TMB")}, this will be set to appropriate defaults which you can see via \code{fit$control}. But if you call \code{estimate_marss()} with a marssMLE object from a call such as \code{MARSS(data, method="kem")} (so not a TMB method), you will need to set \code{opt.control} if you want values different from the base defaults for those functions. Note as a shortcut for \code{nlminb()}, you can set both \code{eval.max}, \code{iter.max} to the same value with \code{opt.control=list(maxit=1000)}. Note, if you pass in \code{opt.control}, this will replace all values currently in \code{fit$control} that are associated with the optimizer function.

The defaults set in \code{\link[MARSS:MARSS]{MARSS::MARSS()}} are
\itemize{
\item \code{nlminb}: \code{eval.max = 5000}, \code{iter.max = 5000} and \code{trace = 0}.
\item \code{optim}: \code{maxit = 5000} and \code{trace = 0}
}

All other controls for the optimization function are left at NULL. You can
set other controls in the call \code{MARSS(..., control=list(...))}.
}
\seealso{
\code{\link[MARSS:MARSSoptim]{MARSS::MARSSoptim()}}, \code{\link[MARSS:MARSSkem]{MARSS::MARSSkem()}}
}
\author{
Eli Holmes.
}
