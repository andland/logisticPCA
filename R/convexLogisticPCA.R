#' @title Convex Logistic Principal Component Analysis
#' 
#' @description 
#' Dimension reduction for binary data by extending Pearson's
#' PCA formulation to minimize Binomial deviance. The convex relaxation
#' to projection matrices, the Fantope, is used.
#' 
#' @param x matrix with all binary entries
#' @param k number of principal components to return
#' @param M value to approximate the saturated model
#' @param quiet logical; whether the calculation should give feedback
#' @param use_irlba logical; if \code{TRUE}, the function uses the irlba package 
#'   to more quickly calculate the eigen-decomposition
#' @param max_iters number of maximum iterations
#' @param conv_criteria convergence criteria. The difference between average deviance
#'   in successive iterations
#' @param random_start logical; whether to randomly inititalize the parameters. If \code{FALSE},
#'   function will use an eigen-decomposition as starting value
#' @param start_H starting value for the Fantope matrix
#' @param mu main effects vector. Only used if \code{main_effects = TRUE}
#' @param main_effects logical; whether to include main effects in the model
#' @param ss_factor step size multiplier. Amount by which to multiply the step size. Quadratic 
#'   convergence can be proven for \code{ss_factor = 1}, but I have found higher values sometimes work 
#'   better. The default is \code{ss_factor = 4}. If it is not converging, try \code{ss_factor = 1}.
#' 
#' @return An S3 object of class \code{clpca} which is a list with the
#' following components:
#' \item{mu}{the main effects}
#' \item{H}{a rank \code{k} Fantope matrix}
#' \item{U}{a \code{k}-dimentional orthonormal matrix with the loadings}
#' \item{PCs}{the princial components}
#' \item{M}{the parameter inputed}
#' \item{iters}{number of iterations required for convergence}
#' \item{loss_trace}{the trace of the average negative log likelihood using the Fantope matrix}
#' \item{proj_loss_trace}{the trace of the average negative log likelihood using the projection matrix}
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
#' # run convex logistic PCA on it
#' clpca = convexLogisticPCA(mat, k = 1, M = 4)
#' @export
convexLogisticPCA <- function(x, k = 2, M = 4, quiet = TRUE, use_irlba = FALSE, 
                              max_iters = 1000, conv_criteria = 1e-7, random_start = FALSE,
                              start_H, mu, main_effects = TRUE, ss_factor = 4) {
  x = as.matrix(x)
  n = nrow(x)
  d = ncol(x)
  q = 2 * x - 1
  q[is.na(q)] <- 0
  eta = q * abs(M)
  
  if (main_effects) {
    if (missing(mu)) {
      if (any(colMeans(x, na.rm = TRUE) == 0 | colMeans(x, na.rm = TRUE) == 1)) {
        stop("Some column(s) all 0 or 1. Remove and try again.")
      }
      mu = log(colMeans(x)) / log(1 - colMeans(x))
    }
  } else {
    mu = rep(0, d)
  }
  
  if (!missing(start_H)) {
    HU = project.Fantope(start_H, k)
    H = HU$H
  } else if (random_start) {
    U = matrix(rnorm(d * d), d, d)
    U = qr.Q(qr(U))
    HU = project.Fantope(U %*% t(U), k)
    H = HU$H
  } else {
    if (use_irlba) {
      udv = irlba(scale(q, center = main_effects, scale = F), nu = k, nv = k)
    } else {
      udv = svd(scale(q, center = main_effects, scale = F), nu = k, nv = k)
    }
    HU = project.Fantope(udv$v[, 1:k] %*% t(udv$v[, 1:k]), k)
    H = HU$H
  }
  
  mu_mat = outer(rep(1, n), mu)
  eta_centered = scale(eta, mu, FALSE)
  etatX = t(eta_centered) %*% x
  theta = mu_mat + eta_centered %*% H
  
  loglike = sum(log(inv.logit.mat(q * theta))[q != 0])
  min_loss = -loglike / sum(q != 0)
  best_HU = HU
  best_loglike = loglike
  if (!quiet) {
    cat(0,"  ", min_loss, "\n")
  }
  
  loss_trace <- proj_loss_trace <- numeric(max_iters)
  loss_trace[1] <- proj_loss_trace[1] <- min_loss
  
  H_lag = H
  for (m in 1:max_iters) {
    y = H + (m - 2) / (m + 1) * (H - H_lag)
    y = H
    H_lag = H
    # y = H
    step = 2 / (M^2 * n * d) * ss_factor
    
    etatP = t(eta_centered) %*% inv.logit.mat(mu_mat + eta_centered %*% y)
    deriv = etatX - etatP
    deriv = deriv + t(deriv) - diag(diag(deriv))
    
    H = y + step * deriv
    HU = project.Fantope(H, k)
    H = HU$H
    
    theta = mu_mat + eta_centered %*% H
    loglike = log_like_Bernoulli(q = q, theta = theta)
    loss_trace[m + 1] = -loglike / sum(q != 0)
    
    proj_theta = mu_mat + eta_centered %*% HU$U %*% t(HU$U)
    proj_loglike = log_like_Bernoulli(q = q, theta = proj_theta)
    proj_loss_trace[m + 1] = -proj_loglike / sum(q!=0)
    
    if (!quiet) {
      cat(m,"  ",loss_trace[m + 1],"  ",proj_loss_trace[m+1],"\n")
    }
    if (loss_trace[m+1] < min_loss) {
      min_loss = loss_trace[m + 1]
      best_HU = HU
      best_loglike = loglike
    }
    if (abs(loss_trace[m+1]-loss_trace[m]) < conv_criteria | min_loss == 0) {
      break
    }
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
  
  object = list(mu = mu,
                H = best_HU$H,
                U = best_HU$U,
                PCs = eta_centered %*% best_HU$U,
                M = M,
                iters = m,
                loss_trace = loss_trace[1:(m + 1)],
                proj_loss_trace = proj_loss_trace[1:(m + 1)],
                prop_deviance_expl = 1 - best_loglike / null_loglike)
  class(object) <- "clpca"
  return(object)
}

#' @title Project onto the Fantope
#' 
#' @description 
#' Project a symmetric matrix onto the convex set of the rank k Fantope
#' 
#' @param x a symmetric matrix
#' @param k the rank of the Fantope desired
#' 
#' @return 
#' \item{H}{a rank \code{k} Fantope matrix}
#' \item{U}{a \code{k}-dimentional orthonormal matrix with the first \code{k} eigenvectors of \code{H}}
#' @export
project.Fantope <- function(x, k) {
  eig = eigen(x, symmetric = TRUE)
  vals = eig$values
  lower = vals[length(vals)] - k / length(vals)
  upper = max(vals)
  while(TRUE) {
    theta = (lower+upper) / 2
    sum.eig.vals = sum(pmin(pmax(vals - theta, 0), 1))
    if (abs(sum.eig.vals-k) < 1e-10) {
      break
    } else if (sum.eig.vals>k) {
      lower = theta
    } else {
      upper = theta
    }
  }
  vals = pmin(pmax(vals - theta, 0), 1)
  return(list(H = eig$vectors %*% diag(vals) %*% t(eig$vectors),
              U = matrix(eig$vectors[, 1:k], nrow(x), k)))
}

#' @title Predict Convex Logistic PCA scores or reconstruction on new data
#' 
#' @param object convex logistic PCA object
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
#' clpca = convexLogisticPCA(mat, k = 1, M = 4, main_effects = FALSE)
#' 
#' PCs = predict(clpca, mat_new)
#' @export
predict.clpca <- function(object, newdata, type = c("PCs", "link", "response"), ...) {
  type = match.arg(type)
  
  if (type == "PCs") {
    if (missing(newdata)) {
      PCs = object$PCs
    } else {
      eta = ((as.matrix(newdata) * 2) - 1) * object$M
      PCs = scale(eta, center = object$mu, scale = FALSE) %*% object$U
    }
    return(PCs)
  } else {
    eta = object$M * (2 * as.matrix(newdata) - 1)
    theta = outer(rep(1, nrow(eta)), object$mu) + scale(eta, object$mu, FALSE) %*% object$H
    if (type == "link") {
      return(theta)
    } else {
      return(inv.logit.mat(theta))
    }
  }
}

#' @title Plot convex logistic PCA
#' 
#' @description 
#' Plots the results of a convex logistic PCA
#' 
#' @param object convex logistic PCA object
#' @param type the type of plot \code{type = "trace"} plots the algorithms progress by
#' iteration, \code{type = "loadings"} plots the loadings first 2 PC loadings
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
#' # run convex logistic PCA on it
#' clpca = convexLogisticPCA(mat, k = 2, M = 4, main_effects = FALSE)
#' 
#' \dontrun{
#' plot(clpca)
#' }
#' @export
plot.clpca <- function(object, type = c("trace", "loadings"), ...) {
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

#' @title Print convex logistic PCA object
#' 
#' @param x convex logistic PCA object
#' @param ... Additional arguments
#' 
#' @export
print.clpca <- function(x, ...) {
  cat(nrow(x$PCs), "rows and ")
  cat(nrow(x$H), "columns\n")
  cat("Rank", ncol(x$U), "Fantope solution with M =", x$M, "\n")
  cat("\n")
  cat(round(x$prop_deviance_expl * 100, 1), "% of deviance explained\n", sep = "")
  cat(x$iters, "iterations to converge\n")
  
  invisible(x)
}
