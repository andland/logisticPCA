#' @title Logistic Singular Value Decomposition
#' 
#' @description 
#' Dimension reduction for binary data by extending SVD to 
#' minimize binomial deviance.
#' 
#' @param x matrix with all binary entries
#' @param k rank of the SVD
#' @param quiet logical; whether the calculation should give feedback
#' @param use_irlba logical; if \code{TRUE}, the function uses the irlba package 
#'   to more quickly calculate the SVD. When the number of columns is small, 
#'   the approximation may be less accurate
#' @param max_iters number of maximum iterations
#' @param conv_criteria convergence criteria. The difference between average deviance
#'   in successive iterations
#' @param random_start logical; whether to randomly inititalize the parameters. If \code{FALSE},
#'   algorithm will use an SVD as starting value
#' @param start_A starting value for the left singular vectors
#' @param start_B starting value for the right singular vectors
#' @param start_mu starting value for mu. Only used if \code{main_effects = TRUE}
#' @param main_effects logical; whether to include main effects in the model
#' 
#' @return An S3 object of class \code{lsvd} which is a list with the
#' following components:
#' \item{mu}{the main effects}
#' \item{A}{a \code{k}-dimentional orthogonal matrix with the scaled left singular vectors}
#' \item{B}{a \code{k}-dimentional orthonormal matrix with the right singular vectors}
#' \item{iters}{number of iterations required for convergence}
#' \item{loss_trace}{the trace of the average negative log likelihood of the algorithm. 
#'    Should be non-increasing}
#' \item{prop_deviance_expl}{the proportion of deviance explained by this model.
#'    If \code{main_effects = TRUE}, the null model is just the main effects, otherwise 
#'    the null model estimates 0 for all natural parameters.}
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
#' # run logistic SVD on it
#' lsvd = logisticSVD(mat, k = 1, main_effects = FALSE, use_irlba = FALSE)
#' 
#' # Logistic SVD likely does a better job finding latent features
#' # than standard SVD
#' plot(svd(mat_logit)$u[, 1], lsvd$A[, 1])
#' plot(svd(mat_logit)$u[, 1], svd(mat)$u[, 1])
#' @export
logisticSVD <- function(x, k = 2, quiet = TRUE, max_iters = 1000, conv_criteria = 1e-5,
                        random_start = FALSE, start_A, start_B, start_mu, 
                        use_irlba = TRUE, main_effects = TRUE) {
  # TODO: Add ALS option?
  use_irlba = use_irlba && requireNamespace("irlba", quietly = TRUE)
  q = 2 * as.matrix(x) - 1
  q[is.na(q)] <- 0 # forces Z to be equal to theta when data is missing
  n = nrow(q)
  d = ncol(q)
  
  # Initialize #
  ##################
  if (!random_start) {
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
  if (!missing(start_mu) && main_effects)
    mu = start_mu
  
  # row.names(A) = row.names(x); row.names(B) = colnames(x)
  loss_trace = numeric(max_iters + 1)
  theta = outer(rep(1, n), mu) + tcrossprod(A, B)
  loglike <- log_like_Bernoulli(q = q, theta = theta)
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
    
    Z = as.matrix(theta + 4 * q * (1 - inv.logit.mat(q * theta)))
    if (main_effects) {
      mu = as.numeric(colMeans(Z))
    }
    
    if (use_irlba) {
      udv = irlba::irlba(scale(Z, center = main_effects, scale = FALSE), nu = max(k, 2), nv = max(k, 2))
    } else {
      udv = svd(scale(Z, center = main_effects, scale = FALSE))
    }
    
    # this is faster than A = sweep(udv$u, 2, udv$d, "*")
    A = matrix(udv$u[, 1:k], n, k) %*% diag(udv$d[1:k], nrow = k, ncol = k)
    B = matrix(udv$v[, 1:k], d, k)
    
    theta = outer(rep(1, n), mu) + tcrossprod(A, B)
    loglike <- log_like_Bernoulli(q = q, theta = theta)
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
      if ((loss_trace[m] - loss_trace[m+1]) < conv_criteria)
        break
    }
    if (m == max_iters) {
      warning("Algorithm ran ", max_iters, " iterations without converging.
              You may want to run it longer.")
    }
  }
  if (loss_trace[m] < loss_trace[m+1]) {
    mu = last_mu
    A = last_A
    B = last_B
    m = m - 1
    
    warning("Algorithm stopped because deviance increased.\nThis should not happen!
            Try rerunning with use_irlba = FALSE")
  }
  
  
  # calculate the null log likelihood for % deviance explained
  if (main_effects) {
    null_proportions = colMeans(x, na.rm = TRUE)
  } else {
    null_proportions = rep(0.5, d)
  }
  null_loglikes <- null_proportions * log(null_proportions) + 
    (1 - null_proportions) * log(1 - null_proportions)
  null_loglike = sum((null_loglikes * colSums(q!=0))[!(null_proportions %in% c(0, 1))])
  
  object = list(mu = mu,
                A = A,
                B = B,
                iters = m,
                loss_trace = loss_trace[1:(m+1)],
                prop_deviance_expl = 1 - loglike / null_loglike)
  class(object) <- "lsvd"
  object
}

#' @title Predict Logistic SVD left singular values or reconstruction on new data
#' 
#' @param object logistic SVD object
#' @param newdata matrix with all binary entries. If missing, will use the 
#'  data that \code{object} was fit on
#' @param quiet logical; whether the calculation should give feedback
#' @param max_iters number of maximum iterations
#' @param conv_criteria convergence criteria. The difference between average deviance
#'   in successive iterations
#' @param random_start logical; whether to randomly inititalize the parameters. If \code{FALSE},
#'   algorithm implicitly starts \code{A} with 0 matrix
#' @param start_A starting value for the left singular vectors
#' @param type the type of fitting required. \code{type = "PCs"} gives the left singular vectors, 
#'  \code{type = "link"} gives matrix on the logit scale and \code{type = "response"} 
#'  gives matrix on the probability scale
#' @param ... Additional arguments
#' 
#' @details
#' Minimizes binomial deviance for new data by finding the optimal left singular vector
#' matrix (\code{A}), given \code{B} and \code{mu}. Assumes the columns of the right 
#' singular vector matrix (\code{B}) are orthonormal.
#' 
#' @examples
#' # construct a low rank matrices in the logit scale
#' rows = 100
#' cols = 10
#' set.seed(1)
#' loadings = rnorm(cols)
#' mat_logit = outer(rnorm(rows), loadings)
#' mat_logit_new = outer(rnorm(rows), loadings)
#' 
#' # convert to a binary matrix
#' mat = (matrix(runif(rows * cols), rows, cols) <= inv.logit.mat(mat_logit)) * 1.0
#' mat_new = (matrix(runif(rows * cols), rows, cols) <= inv.logit.mat(mat_logit_new)) * 1.0
#' 
#' # run logistic PCA on it
#' lsvd = logisticSVD(mat, k = 1, main_effects = FALSE, use_irlba = FALSE)
#' 
#' A_new = predict(lsvd, mat_new)
#' @export
predict.lsvd <- function(object, newdata, quiet = TRUE, max_iters = 1000, conv_criteria = 1e-5,
                         random_start = FALSE, start_A, type = c("PCs", "link", "response"), ...) {
  # TODO: glm option?
  type = match.arg(type)
  
  if (missing(newdata)) {
    A = object$A
  } else {
    x = as.matrix(newdata)
    q = 2* x - 1
    q[is.na(q)] <- 0 # forces Z to be equal to theta when data is missing
    n = nrow(q)
    d = ncol(q)
    k = ncol(object$B)
    
    mu = object$mu
    B = object$B
    mu_mat = outer(rep(1,n),mu)
    if (!missing(start_A)) {
      A = start_A
    } else {
      if (!random_start) {
        # assumes A is initially matrix of 0's and B is orthonormal
        A = 4 * (x - inv.logit.mat(mu_mat)) %*% B
      } else {
        A = matrix(runif(n * k, -1, 1), n, k) 
      }
    }
    
    loss_trace = numeric(max_iters)
    
    for (m in 1:max_iters) {
      last_A = A
      
      theta = mu_mat + tcrossprod(A, B)
      Z = as.matrix(theta + 4*q*(1 - inv.logit.mat(q * theta))) - mu_mat
      
      # assumes columns of B are orthonormal
      A = Z %*% B
      
      loglike = sum(log(inv.logit.mat(q * (mu_mat + tcrossprod(A, B))))[q != 0])
      loss_trace[m] = (-loglike) / sum(q != 0)
      
      if (!quiet) 
        cat(m," ",loss_trace[m], "\n")
      
      if (m > 4) {
        if ((loss_trace[m - 1] - loss_trace[m]) < conv_criteria)
          break
      }
    }
    if (loss_trace[m - 1] < loss_trace[m]) {
      A = last_A
      m = m - 1
      
      warning("Algorithm stopped because deviance increased.\nThis should not happen!")
    }
  }
  
  if (type == "PCs") {
    A
  } else {
    object$A = A
    fitted(object, type, ...)
  }
}

#' @title Fitted values using logistic SVD
#' 
#' @description 
#' Fit a lower dimentional representation of the binary matrix using logistic SVD
#' 
#' @param object logistic SVD object
#' @param type the type of fitting required. \code{type = "link"} gives output on the logit scale and
#'  \code{type = "response"} gives output on the probability scale
#' @param ... Additional arguments
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
#' # run logistic SVD on it
#' lsvd = logisticSVD(mat, k = 1, main_effects = FALSE, use_irlba = FALSE)
#' 
#' # construct fitted probability matrix
#' fit = fitted(lsvd, type = "response")
#' @export
fitted.lsvd <- function(object, type = c("link", "response"), ...) {
  type = match.arg(type)
  n = nrow(object$A)
  
  theta = outer(rep(1, n), object$mu) + tcrossprod(object$A, object$B)
  
  if (type == "link") {
    return(theta)
  } else if (type == "response") {
    return(inv.logit.mat(theta))
  }
}

#' @title Plot logistic SVD
#' 
#' @description 
#' Plots the results of a logistic SVD
#' 
#' @param object logistic SVD object
#' @param type the type of plot \code{type = "trace"} plots the algorithms progress by
#' iteration, \code{type = "loadings"} plots the first 2 principal component
#' loadings, \code{type = "scores"} plots the loadings first 2 principal component scores
#' @param ... Additional arguments
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
#' # run logistic SVD on it
#' lsvd = logisticSVD(mat, k = 2, main_effects = FALSE, use_irlba = FALSE)
#' 
#' \dontrun{
#' plot(lsvd)
#' }
#' @export
plot.lsvd <- function(object, type = c("trace", "loadings", "scores"), ...) {
  library("ggplot2")
  type = match.arg(type)
  
  if (type == "trace") {
    df = data.frame(Iteration = 0:object$iters,
                    NegativeLogLikelihood = object$loss_trace)
    p <- ggplot2::ggplot(df, aes(Iteration, NegativeLogLikelihood)) + 
      geom_line()
  } else if (type == "loadings") {
    df = data.frame(object$B)
    colnames(df) <- paste0("PC", 1:ncol(df))
    if (ncol(df) == 1) {
      p <- ggplot2::qplot(PC1, 0, data = df, ylab = NULL)
    } else {
      p <- ggplot2::ggplot(df, aes(PC1, PC2)) + geom_point()
    }
  } else if (type == "scores") {
    df = data.frame(object$A)
    colnames(df) <- paste0("PC", 1:ncol(df))
    if (ncol(df) == 1) {
      p <- ggplot2::qplot(PC1, 0, data = df, ylab = NULL)
    } else {
      p <- ggplot2::ggplot(df, aes(PC1, PC2)) + geom_point()
    }
  }
  
  return(p)
}

#' @title Print logistic SVD object
#' 
#' @param x logistic SVD object
#' @param ... Additional arguments
#' 
#' @export
print.lsvd <- function(x, ...) {
  cat(nrow(x$A), "rows and ")
  cat(nrow(x$B), "columns\n")
  cat("Rank", ncol(x$B), "solution\n")
  cat("\n")
  cat(round(x$prop_deviance_expl * 100, 1), "% of deviance explained\n", sep = "")
  cat(x$iters, "iterations to converge\n")
  
  invisible(x)
}

#' @title CV for logistic SVD
#' 
#' @description 
#' Run cross validation on dimension for logistic SVD
#' 
#' @param x matrix with all binary entries
#' @param ks the different dimensions \code{k} to try
#' @param folds if \code{folds} is a scalar, then it is the number of folds. If 
#'  it is a vector, it should be the same length as the number of rows in \code{x}
#' @param quiet logical; whether the function should display progress
#' @param ... Additional arguments passed to logisticSVD
#' 
#' @return A matrix of the CV log likelihood with \code{k} in rows
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
#' \dontrun{
#' loglikes = cv.lsvd(mat, ks = 1:9)
#' plot(loglikes)
#' }
#' @export
cv.lsvd <- function(x, ks, folds = 5, quiet = TRUE, ...) {
  q = 2 * as.matrix(x) - 1
  q[is.na(q)] <- 0
  
  if (length(folds) > 1) {
    # does this work if factor?
    if (length(unique(folds)) <= 1) {
      stop("If inputing CV split, must be more than one level")
    }
    if (length(folds) != nrow(x)) {
      stop("if folds is a vector, it should be of same length as nrow(x)")
    }
    cv = folds
  } else {
    cv = sample(1:folds, nrow(q), replace = TRUE)
  }
  
  log_likes = matrix(0, length(ks), 1,
                     dimnames = list(k = ks, M = "LSVD"))
  for (k in ks) {
    if (!quiet) {
      cat("k =", k, "\n")
    }
    for (c in unique(cv)) {
      lsvd = logisticSVD(x[c != cv, ], k = k, ...)
      pred_theta = predict(lsvd, newdat = x[c == cv, ], type = "link")
      log_likes[k == ks] = log_likes[k == ks] + 
        log_like_Bernoulli(q = q[c == cv, ], theta = pred_theta)
      #       log_likes[k == ks] = log_likes[k == ks] + 
      #         sum(log(inv.logit.mat(q[c == cv, ] * pred_theta)))
    }
  }
  class(log_likes) <- c("matrix", "cv.lpca")
  which_max = which.max(log_likes)
  if (!quiet) {
    cat("Best: k =", ks[which_max], "\n")
  }
  
  return(log_likes)
}