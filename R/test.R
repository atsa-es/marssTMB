#' Fit a univariate MARSS model with TMB.
#' 
#' @param y the data. Can have NAs.
#' @param estimate_drift estimate the u parameter
#' @param estimate_rho estimate b parameter
#'
#' @return A data frame with estimates and se's
#' @example inst/examples/uni_example.R
#' @author Eric Ward and edited by Eli Holmes.
#' @export
test <- function(){
  library(MARSS)
  x = cumsum(rnorm(30))
  y = x + rnorm(length(x), 0, 0.01)
  estimate_drift = TRUE # U in MARSS
  estimate_rho = FALSE # AR(1) parameter, b in MARSS
  

parameters <- list(
  log_obs_sd = 0,
  log_pro_sd = 0,
  x = rep(0, length(y)),
  u = 0,
  #x0 = 0,
  logit_rho = 0
)

# Map off parameters not being estimated
tmb_map <- list(x = as.factor(c(NA,1:(length(y)-1))))
if(estimate_drift == FALSE) tmb_map <- c(tmb_map, list(u = factor(NA)))
if(estimate_rho == FALSE) tmb_map <- c(tmb_map, list(logit_rho = factor(NA)))

dat <- t(harborSealWA); dat <- dat[2:4, ]
fit <- MARSS(dat, silent=TRUE)
# Create TMB data
free = lapply(fit$marss$free, function(x){new <- matrix(x, dim(x)[1]*max(dim(x)[2],1), dim(x)[3]); new[is.na(new)] <- 0; new})
fixed = lapply(fit$marss$fixed, function(x){matrix(x, dim(x)[1], dim(x)[3])})
pars = lapply(fit$par, function(x){new <- matrix(x, max(dim(x)[1],1), dim(x)[2]); new[is.na(new)]<- 0; new})

par_dims <- attr(fit$marss, "model.dims")[names(free)]
par_dims = lapply(par_dims, as.integer)
data_list <- list(Y = y, n = length(y),
                  est_drift = as.numeric(estimate_drift),
                  est_rho = as.numeric(estimate_rho),
                  keep = ifelse(!is.na(y),1,0),
                  free = free,
                  fixed = fixed,
                  par = pars,
                  par_dims = par_dims,
                  model = "test")

# Create object for fitting
obj <- TMB::MakeADFun(
  data = data_list,
  map = tmb_map,
  random = "x",
  parameters = parameters,
  DLL = "marssTMB_TMBExports",
  silent = TRUE
)


return(obj$report())
}