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
  return(list(map.diag = diag.factor, map.offdiag = offdiag.factor, map.matrix=mat, raw.matrix=par.as.list))
}

#' Create maps for diagonal and unconstrained varcov matrices
#' 
#' The matrix is split into the diagonal and offdiag components.
#' factor NAs map/mask the values that are fixed and should not be estimated
#' @param x a marssMLE object
#' @param elem the variance-covariance matrix: Q, R or V0
#' @return a list with the factors for the diagonal and offdiagonals. NAs are elements that are not estimated (fixed).
#' @author Eli Holmes and inspired by Tim Cline's code
#' @export
create.elem.maps <- function(x, elem="Z"){
    fixed <- x$model$fixed[[elem]]
    free <- x$model$free[[elem]]
    par.names <- colnames(free)
    var.dim <- attr(x$model, "model.dims")[[elem]][1:2]
    par.as.list <- MARSS:::fixed.free.to.formula(fixed[, , 1, drop = FALSE], free[, , 1, drop = FALSE], var.dim)
    mat <- lapply(par.as.list, function(x){ifelse(!is.character(x), NA, x)}) |>
      matrix(var.dim)
    return(list(map = mat |> unlist() |> factor(levels=par.names), map.matrix=mat, raw.matrix=par.as.list))
}

#' Take one 3rd element from a 3D matrix
#' 
#' Make sure that R does not turn the matrix into a vector
#' 
#' @param x 3D array
#' @param t what 3rd element to take
#' @return a 2D matrix
#' @author Eli Holmes
#' @keywords internal
sub3D <- function(x, t = 1) {
  x.dims <- dim(x)
  if (x.dims[1] != 1 & x.dims[2] != 1) {
    x <- x[, , t]
    return(x)
  } else {
    x <- x[, , t, drop = FALSE]
    dimns <- attr(x, "dimnames")[1:2]
    attr(x, "dim") <- attr(x, "dim")[1:2]
    attr(x, "dimnames") <- dimns
    return(x)
  }
}

#' Vectorize a 2D or 3D matrix
#' 
#' Make sure that R does not turn the matrix into a vector
#' 
#' @param x 2D matrix
#' @return a column matrix
#' @author Eli Holmes
#' @keywords internal
vec <- function(x) {
    attr(x, "dim") <- c(length(x), 1)
    return(x)
}
