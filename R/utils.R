#' Create maps for diagonal and unconstrained varcov matrices
#' 
#' The matrix is split into the diagonal and offdiag components.
#' factor NAs map/mask the values that are fixed and should not be estimated
#' @param x a marssMLE object
#' @param elem the variance-covariance matrix: Q, R or V0
#' @return a list with the factors for the diagonal and offdiagonals. NAs are elements that are not estimated (fixed).
#' @author Eli Holmes and inspired by Tim Cline's code
#' @export
create.varcov.maps <- function(x, elem){
  fixed <- x$model$fixed[[elem]]
  free <- x$model$free[[elem]]
  var.dim <- attr(x$model, "model.dims")[[elem]][1:2]
  par.as.list <- MARSS:::fixed.free.to.formula(fixed[, , 1, drop = FALSE], free[, , 1, drop = FALSE], var.dim)
  mat <- lapply(par.as.list, function(x){ifelse(!is.character(x), NA, x)}) |>
    matrix(var.dim)
  diag.factor <- diag(mat) |> unlist() |> factor()
  offdiag.factor <- mat[upper.tri(mat)] |> unlist() |> factor()
  return(list(diag.factor = diag.factor, offdiag.factor = offdiag.factor, matrix=mat))
}