#' Return R, Q and V0 with free and par for chol of matrix
#'
#' The free, fixed and par (and start) return the vec of the var-cov matrix M:
#' fixed + free%*%p = vec(M)
#' This function transforms the free, fixed, par and start so that they return the chol of the varcov matrix: vec(chol(M)).
#'
#' Why does this function solve for the new par rather than just putting the
#' chol values in? Because the user might have scrambled the order of the
#' parameter list. I won't know which par correspond to which elements of the
#' var-cov matrix. The free matrix is what does the sorting/permutation of the
#' estimated parameters into their correct places. Also the matrix might be
#' block-diagonal, partially fixed, etc.
#'
#' Restrictions
#' * Q and R cannot be time-varying (at the moment)
#'
#' @param MLEobj A properly formatted MARSS model as output by [MARSS()]
#'
#' @return A new [MARSS::marssMLE] object with new fixed, free, par and start

#' @author Eli Holmes.
#' @example inst/examples/var_to_cholvar.R
#' @export
var_decompose <- function(MLEobj) {
  for (elem in c("Q", "R", "V0")) {
    d <- MLEobj[["marss"]][["free"]][[elem]]
    f <- MLEobj[["marss"]][["fixed"]][[elem]]
    if (dim(d)[3] != 1) stop("var_to_cholvar assumes that matrix is time-invariant")
    d <- sub3D(d, t = 1)
    f <- sub3D(f, t = 1)
    # Set lower tri of free to 0
    r <- sqrt(dim(d)[1])
    lt <- matrix(1:(r * r), r, r)
    lt <- lt[lower.tri(lt)] # rows of free assoc with the lower tri
    d[lt, ] <- 0

    for (what in c("par", "start")) {
      if (is.null(MLEobj[[what]])) next
      the.par <- coef(MLEobj, type = "matrix", what = what)[[elem]]
      diag.par <- diag(the.par)
      is.zero <- diag(the.par) == 0 # where the 0s on diagonal are
      # so the chol doesn't fail if there are zeros on the diagonal
      if (any(is.zero)) diag(the.par)[is.zero] <- 1
      the.par <- chol(the.par) # convert to transpose of chol
      if (any(is.zero)) diag(the.par)[is.zero] <- 0 # set back to 0

      # from f+Dm=M so m = solve(crossprod(d))%*%t(d)%*%(vec(the.par)-f)
      # no f needed since if d!=0,then f==0. if f!-0, then d==0.
      if (dim(d)[2] != 0) { # something is estimated
        MLEobj[[what]][[elem]] <- solve(crossprod(d)) %*% t(d) %*% vec(the.par)
      } else {
        MLEobj[[what]][[elem]] <- matrix(0, 0, 1)
      }
    }
    # update so that lower tri is zero
    if (dim(d)[2] != 0) {
      MLEobj[["marss"]][["free"]][[elem]][, , 1] <- d
      MLEobj[["model"]][["free"]][[elem]] <- MLEobj[["marss"]][["free"]][[elem]]
    }
    # update to fixed chol values
    MLEobj[["marss"]][["fixed"]][[elem]][, , 1] <- (diag(1, r * r) - tcrossprod(d)) %*% vec(the.par)
    MLEobj[["model"]][["fixed"]][[elem]] <- MLEobj[["marss"]][["fixed"]][[elem]]
  }
  return(MLEobj)
}
