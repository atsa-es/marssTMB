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

library(MARSS)
dat <- t(harborSealWA); dat <- dat[2:4, ]
MLEobj <- MARSS(dat, silent=TRUE, model=list(tinitx=1), form="dfa")
MODELobj <- MLEobj$model

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
  par_dims[[dname]] <- c(mn,1,1)
  par_dims[[oname]] <- c(mn * (mn - 1) / 2, 1, 1)
}
free <- free[!(names(free) %in% c("Q", "R"))]
fixed <- fixed[!(names(fixed) %in% c("Q", "R"))]
pars <- pars[!(names(pars) %in% c("Q", "R"))]

tfixed <- unlist(lapply(fixed, function(x){dim(x)[3]}))
tfree <- unlist(lapply(free, function(x){dim(x)[3]}))
free = lapply(free, function(x){new <- matrix(x, max(dim(x)[1],1)*max(dim(x)[2],1), dim(x)[3]); new[is.na(new)] <- 0; new})
fixed = lapply(fixed, function(x){matrix(x, dim(x)[1], dim(x)[3])})
pars0 = lapply(pars, function(x){new <- matrix(x, max(dim(x)[1],1), dim(x)[2]); new[is.na(new)]<- 0; new})
par_dims = par_dims[names(fixed)]
par_dims = lapply(par_dims, as.integer)
numpar <- unlist(lapply(pars, nrow))
# a vector of the parameters
pars <- unlist(lapply(pars, function(x){x[,1]}))

data_list <- list(Y = y, n = length(y),
                  est_drift = as.numeric(estimate_drift),
                  est_rho = as.numeric(estimate_rho),
                  keep = ifelse(!is.na(y),1,0),
                  mY =  MLEobj$model[["data"]],
                  X = MLEobj$states,
                  tfixed = tfixed,
                  tfree = tfree, 
                  numpar = numpar, 
                  free = free,
                  fixed = fixed,
                  par_dims = par_dims,
                  pars = pars,
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