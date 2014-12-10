#' @title Logistic Singular Value Decomposition
#' 
#' @description 
#' Dimension reduction for binary data by extending SVD to 
#' minimize binomial deviance.
#' 
#' @param dat matrix with all binary entries
#' @param k rank of the SVD
#' @param quiet logical; whether the calculation should give feedback
#' @param use_irlba logical; if \code{TRUE}, the function uses the irlba package 
#'   to more quickly calculate the SVD
#' @param max_iters number of maximum iterations
#' @param conv_criteria convergence criteria. The difference between average deviance
#'   in successive iterations
#' @param random_start logical; whether to randomly inititalize the parameters. If \code{FALSE},
#'   function will use an SVD as starting value
#' @param start_A starting value for the orthoganal matrix. If missing, initializes 
#'   with first \code{k} left singular vectors of \code{dat}
#' @param start_B starting value for the orthoganal matrix. If missing, initializes 
#'   with first \code{k} right singular vectors of \code{dat}
#' @param start_mu starting value for mu, if \code{main_effects = TRUE}
#' @param main_effects logical; whether to include main effects in the model
#' 
#' @return An S3 object of class \code{lsvd} which is a list with the
#' following components:
#' \item{A}{a \code{k}-dimentional orthogonal matrix with the left singular vectors}
#' \item{B}{a \code{k}-dimentional orthonormal matrix with the right singular vectors}
#' \item{mu}{the main effects}
#' \item{iters}{number of iterations required for convergence}
#' \item{loss_trace}{the loss trace of the algorithm. Should be non-increasing}
#' 
#' @references 
#' de Leeuw, Jan, 2006. Principal component analysis of binary data 
#' by iterated singular value decomposition. Computational Statistics & Data Analysis 
#' 50 (1), 21--39.
#' 
#' Collins, M., Dasgupta, S., & Schapire, R. E., 2001. A generalization of principal 
#' components analysis to the exponential family. In NIPS, 617--624.
#' 
#' @examples
#' # construct a low rank matrix in the logit scale
#' rows = 100
#' cols = 10
#' set.seed(1)
#' mat_logit = outer(rnorm(rows), rnorm(cols))
#' 
#' # generate a binary matrix
#' mat = (matrix(runif(rows * cols), rows, cols) <= inv.logit.mat(mat_logit)) * 1.0
#' 
#' # run logistic PCA on it
#' lsvd = logisticSVD(mat, k = 1, main_effects = FALSE)
#' 
#' # Logistic SVD likely does a better job finding latent features
#' # than standard SVD
#' plot(svd(mat_logit)$u[, 1], lsvd$A[, 1])
#' plot(svd(mat_logit)$u[, 1], svd(mat)$u[, 1])
#' @export
logisticSVD <- function(dat, k = 2, quiet = TRUE, max_iters = 1000, conv_crit=1e-5,
                        randstart = FALSE, start_A, start_B, start_mu, 
                        use_irlba = TRUE, main_effects = TRUE) {
  # TODO: Add ALS option?
  use_irlba = use_irlba && requireNamespace("irlba", quietly = TRUE)
  q = 2 * as.matrix(dat) - 1
  q[is.na(q)] <- 0 # forces x to be equal to theta when data is missing
  n = nrow(q)
  d = ncol(q)
  
  # Initialize #
  ##################
  if (!randstart) {
    if (main_effects) {
      mu = colMeans(4 * q)
    } else {
      mu = rep(0, d)
    }
    if (missing(start_A) | missing(start_B)) {
      if (!quiet) {cat("Initializing SVD... ")}
      if (use_irlba) {
        udv = irlba::irlba(scale(4 * q, center = main_effects, scale = FALSE), nu = max(k, 2), nv = max(k, 2))
      } else {
        udv = svd(scale(4 * q, center = main_effects, scale = FALSE))
      }
      if (!quiet) {cat("Done!\n")}
      A = matrix(udv$u[, 1:k], n, k) %*% diag(udv$d[1:k], nrow = k, ncol = k)
      B = matrix(udv$v[, 1:k], d, k)
    }
  } else {
    if (main_effects) {
      mu = rnorm(d)
    } else {
      mu = rep(0, d)
    }
    A = matrix(runif(n * k, -1, 1), n, k)
    B = matrix(runif(d * k, -1, 1), d, k)
  }
  if (!missing(start_B))
    B = start_B
  if (!missing(start_A))
    A = start_A
  if (!missing(start_mu))
    mu = start_mu
  
  # row.names(A) = row.names(dat); row.names(B) = colnames(dat)
  loss_trace = numeric(max_iters + 1)
  loglike = sum(log(inv.logit.mat(q * (outer(rep(1, n), mu) + tcrossprod(A, B))))[q != 0])
  loss_trace[1] = -loglike / sum(q!=0)
  ptm <- proc.time()
  if (!quiet) {
    cat(0, "  ", loss_trace[1], "")
    cat("0 hours elapsed\n")
  }
  
  for (m in 1:max_iters) {
    last_mu = mu
    last_A = A
    last_B = B
    
    theta = outer(rep(1, n), mu) + tcrossprod(A, B)
    X = as.matrix(theta + 4 * q * (1 - inv.logit.mat(q * theta)))
    if (main_effects) {
      mu = as.numeric(colMeans(X))
    }
    
    if (use_irlba) {
      udv = irlba::irlba(scale(X, center = main_effects, scale = FALSE), nu = max(k, 2), nv = max(k, 2))
    } else {
      udv = svd(scale(X, center = main_effects, scale = FALSE))
    }
    # TODO: would sweep be faster here?
    A = matrix(udv$u[, 1:k], n, k) %*% diag(udv$d[1:k], nrow = k, ncol = k)
    B = matrix(udv$v[, 1:k], d, k)
    
    loglike = sum(log(inv.logit.mat(q * (outer(rep(1, n), mu) + tcrossprod(A, B))))[q != 0])
    loss_trace[m+1] = -loglike / sum(q != 0)
    
    if (!quiet) {
      time_elapsed = as.numeric(proc.time() - ptm)[3]
      tot_time = max_iters / m * time_elapsed
      time_remain = tot_time - time_elapsed
      cat(m, "  ", loss_trace[m+1], "")
      cat(round(time_elapsed / 3600, 1), "hours elapsed. Max", 
          round(time_remain / 3600, 1), "hours remain.\n")
    }
    if (m > 4) {
      if ((loss_trace[m] - loss_trace[m+1]) < conv_crit)
        break
    }
  }
  if (loss_trace[m] < loss_trace[m+1]) {
    mu = last_mu
    A = last_A
    B = last_B
    m = m - 1
    
    warning("Algorithm stopped because deviance increased.\nThis should not happen!")
  }
  
  object = list(mu = mu,
                A = A,
                B = B,
                iters = m,
                loss_trace = loss_trace[1:(m+1)])
}

predict.lsvd <- function(object, newdat, quiet = TRUE, max_iters = 1000, conv_crit = 1e-5,
                         randstart = FALSE, procrustes = FALSE, normalize = FALSE, start_A, ...) {
  # TODO: initialization can be improved. Assume B orthonormal
  q = 2* as.matrix(newdat) - 1
  q[is.na(q)] <- 0 # forces x to be equal to theta when data is missing
  n = nrow(q)
  d = ncol(q)
  k = ncol(object$B)
  
  mu = object$mu
  B = object$B
  if (!missing(start_A)) {
    A = start_A
  } else {
    if (!randstart) {
      udv = svd(scale(q, center = mu, scale = FALSE))
      A = matrix(udv$u[, 1:k], n, k) %*% diag(udv$d[1:k], nrow = k, ncol = k)
    } else {
      A = matrix(runif(n * k, -1, 1), n, k) 
    }
  }
  
  loss_trace = numeric(max_iters)
  
  mu_mat = outer(rep(1,n),mu)
  if (!procrustes) {
    BBtB = B %*% solve(t(B) %*% B)
  }
  for (m in 1:max_iters) {
    last_A = A
    
    theta = mu_mat + tcrossprod(A, B)
    Xstar = as.matrix(theta + 4*q*(1 - inv.logit.mat(q * theta))) - mu_mat
    if (procrustes) {
      M = svd(Xstar %*% B)
      A = M$u %*% t(M$v)
    } else {
      A = Xstar %*% BBtB
    }
    
    loglike = sum(log(inv.logit.mat(q * (outer(rep(1, n), mu) + tcrossprod(A, B))))[q != 0])
    loss_trace[m] = (-loglike) / sum(!is.na(newdat))
    
    if (!quiet) 
      cat(m," ",loss_trace[m], "\n")
    
    if (m > 4) {
      if ((loss_trace[m - 1] - loss_trace[m]) < conv_crit)
        break
    }
  }
  if (loss_trace[m - 1] < loss_trace[m]) {
    A = last_A
    m = m - 1
    
    loglike = sum(log(inv.logit.mat(q * (outer(rep(1, n), mu) + tcrossprod(A, B))))[q != 0])
  }
  if (normalize) {
    A = sweep(A, 2, sqrt(colSums(B^2)), "*")
    B = sweep(B, 2, sqrt(colSums(B^2)), "/")
  }
  
  object = list(A = A,
                iters = m,
                loss_trace = loss_trace[1:m])
  return(object)
}