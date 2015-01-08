#' @title Inverse logit for matrices
#' 
#' @description
#' Apply the inverse logit function to a matrix, element-wise. It 
#' generalizes the \code{inv.logit} function from the \code{gtools} 
#' library to matrices
#'
#' @param x matrix
#' @param min Lower end of logit interval
#' @param max Upper end of logit interval
#' @examples
#' (mat = matrix(rnorm(10 * 5), nrow = 10, ncol = 5))
#' inv.logit.mat(mat)
#' @export
inv.logit.mat <- function(x, min = 0, max = 1) {
  .Call(inv_logit_mat, x, min, max)
}

#' @title Bernoulli Log Likelihood
#' 
#' @description
#' Calculate the Bernoulli log likelihood of matrix
#'
#' @param x matrix with all binary entries
#' @param theta estimated natural parameters with 
#'  same dimensions as x
#' @param q instead of x, you can input matrix q which is 
#'  -1 if \code{x = 0}, 1 if \code{x = 1}, and 0 if \code{is.na(x)}
#' @examples
#' (mat = matrix(rnorm(10 * 5), nrow = 10, ncol = 5))
#' inv.logit.mat(mat)
#' @export
log_like_Bernoulli <- function(x, theta, q) {
  if (missing(x) & !missing(q)) {
    q = 2 * as.matrix(x) - 1
    q[is.na(q)] <- 0
  }
  .Call(compute_loglik, q, theta)
}
