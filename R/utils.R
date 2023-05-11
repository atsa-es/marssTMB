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
