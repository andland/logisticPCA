#' @title Logistic Principal Component Analysis
#' 
#' @description 
#' Dimension reduction for binary data by extending Pearson's
#' PCA formulation to minimize Binomial deviance
#' 
#' @param x matrix with all binary entries
#' @param k number of principal components to return
#' @param M value to approximate the saturated model. If \code{M = 0}, M is solved for
#' @param quiet logical; whether the calculation should give feedback
#' @param use_irlba logical; if \code{TRUE}, the function uses the irlba package 
#'   to more quickly calculate the eigen-decomposition
#' @param max_iters number of maximum iterations
#' @param conv_criteria convergence criteria. The difference between average deviance
#'   in successive iterations
#' @param random_start logical; whether to randomly inititalize the parameters. If \code{FALSE},
#'   function will use an eigen-decomposition as starting value
#' @param start_U starting value for the orthogonal matrix
#' @param start_mu starting value for mu. Only used if \code{main_effects = TRUE}
#' @param main_effects logical; whether to include main effects in the model
#' 
#' @return An S3 object of class \code{lpca} which is a list with the
#' following components:
#' \item{mu}{the main effects}
#' \item{U}{a \code{k}-dimentional orthonormal matrix with the loadings}
#' \item{PCs}{the princial components}
#' \item{M}{the parameter inputed or solved for}
#' \item{iters}{number of iterations required for convergence}
#' \item{loss_trace}{the trace of the average negative log likelihood of the algorithm. 
#'    Should be non-increasing}
#' \item{prop_deviance_expl}{the proportion of deviance explained by this model.
#'    If \code{main_effects = TRUE}, the null model is just the main effects, otherwise 
#'    the null model estimates 0 for all natural parameters.}
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
#' lpca = logisticPCA(mat, k = 1, M = 4, main_effects = FALSE)
#' 
#' # Logistic PCA likely does a better job finding latent features
#' # than standard PCA
#' plot(svd(mat_logit)$u[, 1], lpca$PCs[, 1])
#' plot(svd(mat_logit)$u[, 1], svd(mat)$u[, 1])
#' @export
logisticPCA <- function(x, k = 2, M = 4, quiet = TRUE, use_irlba = FALSE,
                        max_iters = 1000, conv_criteria = 1e-5, random_start = FALSE,
                        start_U, start_mu, main_effects = TRUE) {
  # better name for k
  use_irlba = use_irlba && requireNamespace("irlba", quietly = TRUE)
  q = as.matrix(2 * x - 1)
  q[is.na(q)] <- 0 # forces Z to be equal to theta when data is missing
  n = nrow(q)
  d = ncol(q)
  if (M == 0) {
    M = 4
    solve_M = TRUE
  } else {
    solve_M = FALSE
  }
  # eta = q * abs(M)
  
  if (main_effects) {
    if (!missing(start_mu)) {
      mu = start_mu
    } else {
      mu = colMeans(M * q)
    }
  } else {
    mu = rep(0, d)
  }
  
  # Initialize #
  ##################
  if (!missing(start_U)) {
    U = sweep(start_U, 2, sqrt(colSums(start_U^2)), "/")
  } else if (random_start) {
    U = matrix(rnorm(d * k), d, k)
    U = qr.Q(qr(U))
  } else {
    if (use_irlba) {
      udv = irlba::irlba(scale(q, center = main_effects, scale = FALSE), nu = k, nv = k)
    } else {
      udv = svd(scale(q, center = main_effects, scale = FALSE))
    }
    U = matrix(udv$v[, 1:k], d, k)
  }
  
  # etaTeta = crossprod(eta)
  qTq = crossprod(q)
  
  loss_trace = numeric(max_iters + 1)
  theta = outer(rep(1, n), mu) + scale(M * q, center = mu, scale = FALSE) %*% U %*% t(U)
  loglike <- log_like_Bernoulli(q = q, theta = theta)
  loss_trace[1] = (-loglike) / sum(q!=0)
  ptm <- proc.time()
  
  if (!quiet) {
    cat(0, "  ", loss_trace[1], "")
    cat("0 hours elapsed\n")
  }
  
  for (m in 1:max_iters) {
    last_U = U
    last_M = M
    last_mu = mu
    
    if (solve_M) {
      Phat = inv.logit.mat(theta)
      # Need to check this is correct when missing data
      M_slope = sum(((x - Phat) * (q %*% U %*% t(U)))[q != 0]) 
      M_curve = -sum((Phat * (1 - Phat) * (q %*% U %*% t(U))^2)[q != 0])
      M = M - M_slope / M_curve
      
      theta = outer(rep(1, n), mu) + scale(M * q, center = mu, scale = FALSE) %*% tcrossprod(U)
    }
    
    Z = as.matrix(theta + 4 * q * (1 - inv.logit.mat(q * theta)))
    if (main_effects) {
      mu = as.numeric(colMeans(Z - (M * q) %*% U %*% t(U)))
    }
    
    mat_temp = t(scale(M * q, center = mu, scale = FALSE)) %*% Z
    mat_temp = mat_temp + t(mat_temp) - M^2 * qTq + n * outer(mu, mu)
    
    # irlba sometimes gives poor estimates of e-vectors
    # so I switch to standard eigen if it does
    repeat {
      if (use_irlba) {
        udv = irlba::irlba(mat_temp, nu=k, nv=k, adjust=3)
        U = matrix(udv$u[, 1:k], d, k)
      } else {
        eig = eigen(mat_temp, symmetric=TRUE)
        U = matrix(eig$vectors[, 1:k], d, k)
      }
      
      theta = outer(rep(1, n), mu) + scale(M * q, center = mu, scale = FALSE) %*% tcrossprod(U)
      this_loglike <- log_like_Bernoulli(q = q, theta = theta)
      
      if (!use_irlba | this_loglike>=loglike) {
        loglike = this_loglike
        break
      } else {
        use_irlba = FALSE
        if (!quiet) {cat("Quitting irlba!\n")}
      }
    }
    
    loss_trace[m + 1] = (-loglike) / sum(q!=0)
    
    if (!quiet) {
      time_elapsed = as.numeric(proc.time() - ptm)[3]
      tot_time = max_iters / m * time_elapsed
      time_remain = tot_time - time_elapsed
      cat(m, "  ", loss_trace[m + 1], "")
      cat(round(time_elapsed / 3600, 1), "hours elapsed. Max", round(time_remain / 3600, 1), "hours remain.\n")
    }
    if (m > 4) {
      # when solving for M, the monoticity does not apply
      if (solve_M) {
        if (abs(loss_trace[m] - loss_trace[m+1]) < conv_criteria) {
          break
        }
      } else {
        if ((loss_trace[m] - loss_trace[m+1]) < conv_criteria) {
          break
        }
      }
    }
  }
  
  # test if loss function increases
  if ((loss_trace[m + 1] - loss_trace[m]) > (1e-10)) {
    U = last_U
    mu = last_mu
    M = last_M
    m = m - 1
    
    warning("Algorithm stopped because deviance increased.\nThis should not happen!")
  }
  
  # calculate the null log likelihood for % deviance explained
  # assumes no missing data
  if (main_effects) {
    null_proportions = colMeans(x)
  } else {
    null_proportions = rep(0.5, d)
  }
  null_loglikes <- null_proportions * log(null_proportions) + 
    (1 - null_proportions) * log(1 - null_proportions)
  null_loglike = sum(null_loglikes[!(null_proportions %in% c(0, 1))]) * n
  
  object <- list(mu = mu,
                 U = U,
                 PCs = scale(M * q, center = mu, scale = FALSE) %*% U,
                 M = M,
                 iters = m,
                 loss_trace = loss_trace[1:(m + 1)],
                 prop_deviance_expl = 1 - loglike / null_loglike)
  class(object) <- "lpca"
  object
}


#' @title Predict Logistic PCA scores or reconstruction on new data
#' 
#' @param object logistic PCA object
#' @param newdata matrix with all binary entries. If missing, will use the 
#'  data that \code{object} was fit on
#' @param type the type of fitting required. \code{type = "PCs"} gives the PC scores, 
#'  \code{type = "link"} gives matrix on the logit scale and \code{type = "response"} 
#'  gives matrix on the probability scale
#' @param ... Additional arguments
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
#' lpca = logisticPCA(mat, k = 1, M = 4, main_effects = FALSE)
#' 
#' PCs = predict(lpca, mat_new)
#' @export
predict.lpca <- function(object, newdata, type = c("PCs", "link", "response"), ...) {
  type = match.arg(type)
  
  if (missing(newdata)) {
    PCs = object$PCs
  } else {
    eta = ((as.matrix(newdata) * 2) - 1) * object$M
    PCs = scale(eta, center = object$mu, scale = FALSE) %*% object$U
  }
  
  if (type == "PCs") {
    PCs
  } else {
    object$PCs = PCs
    fitted(object, type, ...)
  }
}

#' @title Fitted values using logistic PCA
#' 
#' @description 
#' Fit a lower dimentional representation of the binary matrix using logistic PCA
#' 
#' @param object logistic PCA object
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
#' # run logistic PCA on it
#' lpca = logisticPCA(mat, k = 1, M = 4, main_effects = FALSE)
#' 
#' # construct fitted probability matrix
#' fit = fitted(lpca, type = "response")
#' @export
fitted.lpca <- function(object, type = c("link", "response"), ...) {
  type = match.arg(type)
  n = nrow(object$PCs)
  
  theta = outer(rep(1, n), object$mu) + tcrossprod(object$PCs, object$U)
  
  if (type == "link") {
    return(theta)
  } else if (type == "response") {
    return(inv.logit.mat(theta))
  }
}

#' @title Plot logistic PCA
#' 
#' @description 
#' Plots the results of a logistic PCA
#' 
#' @param object logistic PCA object
#' @param type the type of plot \code{type = "trace"} plots the algorithms progress by
#' iteration, \code{type = "loadings"} plots the loadings first 2 principal components
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
#' # run logistic PCA on it
#' lpca = logisticPCA(mat, k = 2, M = 4, main_effects = FALSE)
#' 
#' \dontrun{
#' plot(lpca)
#' }
#' @export
plot.lpca <- function(object, type = c("trace", "loadings"), ...) {
  library("ggplot2")
  type = match.arg(type)
  
  if (type == "trace") {
    df = data.frame(Iteration = 0:object$iters,
                    NegativeLogLikelihood = object$loss_trace)
    p <- ggplot2::ggplot(df, aes(Iteration, NegativeLogLikelihood)) + 
      geom_line()
  } else if (type == "loadings") {
    df = data.frame(object$U)
    colnames(df) <- paste0("PC", 1:ncol(df))
    if (ncol(df) == 1) {
      p <- ggplot2::qplot(PC1, 0, data = df, ylab = NULL)
    } else {
      p <- ggplot2::ggplot(df, aes(PC1, PC2)) + geom_point()
    }
  }
  
  return(p)
}

#' @title Print logistic PCA object
#' 
#' @param x logistic PCA object
#' @param ... Additional arguments
#' 
#' @export
print.lpca <- function(x, ...) {
  cat(nrow(x$PCs), "rows and ")
  cat(nrow(x$U), "columns\n")
  cat("Rank", ncol(x$U), "solution with M =", x$M, "\n")
  cat("\n")
  cat(round(x$prop_deviance_expl * 100, 1), "% of deviance explained\n", sep = "")
  cat(x$iters, "iterations to converge\n")
  
  invisible(x)
}

#' @title CV for logistic PCA
#' 
#' @description 
#' Run cross validation on dimension and \code{M} for logistic PCA
#' 
#' @param x matrix with all binary entries
#' @param ks the different dimensions \code{k} to try
#' @param Ms the different approximations to the saturated model \code{M} to try
#' @param folds if \code{folds} is a scalar, then it is the number of folds. If 
#'  it is a vector, it should be the same length as the number of rows in \code{x}
#' @param quiet logical; whether the function should display progress
#' @param ... Additional arguments passed to logisticPCA
#' 
#' @return A matrix of the CV log likelihood with \code{k} in rows and 
#'  \code{M} in columns
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
#' loglikes = cv.lpca(mat, ks = 1:9, Ms = 3:6)
#' plot(loglikes)
#' }
#' @export
cv.lpca <- function(x, ks, Ms = seq(2, 10, by = 2), folds = 5, quiet = TRUE, ...) {
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
  
  log_likes = matrix(0, length(ks), length(Ms),
                     dimnames = list(k = ks, M = Ms))
  for (k in ks) {
    for (M in Ms) {
      if (!quiet) {
        cat("k =", k, "M =", M, "")
      }
      for (c in unique(cv)) {
        if (!quiet) {
          cat(".")
        }
        lpca = logisticPCA(x[c != cv, ], k = k, M = M, ...)
        pred_theta = predict(lpca, newdat = x[c == cv, ], type = "link")
        log_likes[k == ks, M == Ms] = log_likes[k == ks, M == Ms] + 
          log_like_Bernoulli(q = q[c == cv, ], theta = pred_theta)
        #         log_likes[k == ks, M == Ms] = log_likes[k == ks, M == Ms] + 
        #           sum(log(inv.logit.mat(q[c == cv, ] * pred_theta)))
      }
      if (!quiet) {
        cat("", log_likes[k == ks, M == Ms], "\n")
      }
    }
  }
  class(log_likes) <- c("matrix", "cv.lpca")
  which_min = which(log_likes == max(log_likes), arr.ind = TRUE)
  if (!quiet) {
    cat("Best: k =", ks[which_min[1]], "M =", Ms[which_min[2]], "\n")
  }
  
  return(log_likes)
}

#' @title Plot CV for logistic PCA
#' 
#' @description 
#' Plot cross validation results logistic PCA
#' 
#' @param object a \code{cv.lpca} object
#' @param ... Additional arguments
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
#' loglikes = cv.lpca(dat, ks = 1:9, Ms = 3:6)
#' plot(loglikes)
#' }
#' @export
plot.cv.lpca <- function(object, ...) {
  library(ggplot2)
  library(reshape2)
  df = melt(-object, value.name = "NegLogLikelihood")
  if (ncol(object) == 1) {
    df$M = factor(df$M)
    p <- ggplot(df, aes(k, NegLogLikelihood, colour = M)) + 
      geom_line()
  } else {
    df$k = factor(df$k)
    p <- ggplot(df, aes(M, NegLogLikelihood, colour = k)) + 
      geom_line()
  }
  return(p)
}
