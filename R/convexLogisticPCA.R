#' @title Convex Logistic Principal Component Analysis
#'
#' @description
#' Dimensionality reduction for binary data by extending Pearson's
#' PCA formulation to minimize Binomial deviance. The convex relaxation
#' to projection matrices, the Fantope, is used.
#'
#' @param x matrix with all binary entries
#' @param k number of principal components to return
#' @param m value to approximate the saturated model
#' @param quiet logical; whether the calculation should give feedback
#' @param partial_decomp logical; if \code{TRUE}, the function uses the RSpectra package
#'   to quickly initialize \code{H} and project onto the Fantope when \code{ncol(x)} 
#'   is large and \code{k} is small
#' @param max_iters number of maximum iterations
#' @param conv_criteria convergence criteria. The difference between average deviance
#'   in successive iterations
#' @param random_start logical; whether to randomly inititalize the parameters. If \code{FALSE},
#'   function will use an eigen-decomposition as starting value
#' @param start_H starting value for the Fantope matrix
#' @param mu main effects vector. Only used if \code{main_effects = TRUE}
#' @param main_effects logical; whether to include main effects in the model
#' @param ss_factor step size multiplier. Amount by which to multiply the step size. Quadratic
#'   convergence rate can be proven for \code{ss_factor = 1}, but I have found higher values
#'   sometimes work better. The default is \code{ss_factor = 4}.
#'   If it is not converging, try \code{ss_factor = 1}.
#' @param weights an optional matrix of the same size as the \code{x} with non-negative weights
#' @param M depricated. Use \code{m} instead
#'
#' @return An S3 object of class \code{clpca} which is a list with the
#' following components:
#' \item{mu}{the main effects}
#' \item{H}{a rank \code{k} Fantope matrix}
#' \item{U}{a \code{ceiling(k)}-dimentional orthonormal matrix with the loadings}
#' \item{PCs}{the princial component scores}
#' \item{m}{the parameter inputed}
#' \item{iters}{number of iterations required for convergence}
#' \item{loss_trace}{the trace of the average negative log likelihood using the Fantope matrix}
#' \item{proj_loss_trace}{the trace of the average negative log likelihood using the projection matrix}
#' \item{prop_deviance_expl}{the proportion of deviance explained by this model.
#'    If \code{main_effects = TRUE}, the null model is just the main effects, otherwise
#'    the null model estimates 0 for all natural parameters.}
#' \item{rank}{the rank of the Fantope matrix \code{H}}
#' 
#' @references 
#' Landgraf, A.J. & Lee, Y., 2015. Dimensionality reduction for binary data through 
#' the projection of natural parameters. arXiv preprint arXiv:1510.06112.
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
#' clpca = convexLogisticPCA(mat, k = 1, m = 4)
#' @export
convexLogisticPCA <- function(x, k = 2, m = 4, quiet = TRUE, partial_decomp = FALSE,
                              max_iters = 1000, conv_criteria = 1e-6, random_start = FALSE,
                              start_H, mu, main_effects = TRUE, ss_factor = 4, weights, M) {
  if (!missing(M)) {
    m = M
    warning("M is depricated. Use m instead. ",
            "Using m = ", m)
  }
  if (partial_decomp) {
    if (!requireNamespace("RSpectra", quietly = TRUE)) {
      message("RSpectra must be installed to use partial_decomp")
      partial_decomp = FALSE
    }
  }
  
  x = as.matrix(x)
  n = nrow(x)
  d = ncol(x)
  q = 2 * x - 1
  q[is.na(q)] <- 0
  eta = q * abs(m)
  
  if (k >= d & partial_decomp) {
    message("k >= dimension. Setting partial_decomp = FALSE")
    partial_decomp = FALSE
    k = d
  }

  if (missing(weights) || length(weights) == 1) {
    weights = 1.0
    sum_weights = sum(!is.na(x))
  } else {
    if (!all(dim(weights) == dim(x))) {
      stop("x and weights are not the same dimension")
    }
    weights[is.na(x)] <- 0
    if (any(is.na(weights))) {
      stop("Can't have NA in weights")
    }
    if (any(weights < 0)) {
      stop("weights must be non-negative")
    }
    sum_weights = sum(weights)
  }

  if (main_effects) {
    if (missing(mu)) {
      if (length(weights) == 1) {
        x_bar = colMeans(x, na.rm = TRUE)
      } else {
        x_bar = colSums(x * weights, na.rm = TRUE) / colSums(weights, na.rm = TRUE)
      }
      if (any(x_bar == 0 | x_bar == 1)) {
        stop("Some column(s) all 0 or 1. Remove and try again.")
      }
      mu = log(x_bar) - log(1 - x_bar)
    }
  } else {
    mu = rep(0, d)
  }

  if (!missing(start_H)) {
    HU = project.Fantope(start_H, k, partial_decomp = partial_decomp)
    H = HU$H
  } else if (random_start) {
    U = matrix(rnorm(d * d), d, d)
    U = qr.Q(qr(U))
    HU = project.Fantope(U %*% t(U), k, partial_decomp = partial_decomp)
    H = HU$H
  } else {
    if (partial_decomp) {
      udv = RSpectra::svds(scale(q, center = main_effects, scale = F), k = ceiling(k))
    } else {
      udv = svd(scale(q, center = main_effects, scale = F), nu = ceiling(k), nv = ceiling(k))
    }
    H = tcrossprod(udv$v[, 1:ceiling(k), drop = FALSE])
    HU = list(H = H,
         U = udv$v[, 1:ceiling(k), drop = FALSE],
         rank = ceiling(k))
  }

  mu_mat = outer(rep(1, n), mu)
  eta_centered = scale(eta, mu, FALSE)

  # when x is missing eta = mu. So eta_centered is 0
  eta_centered[q == 0] <- 0

  # Lipschitz constant
  L = sum(eta_centered^2) * max(weights) / 2

  # only sum over non-missing x. Equivalent to replacing missing x with 0
  x[q == 0] <- 0

  etatX = t(eta_centered) %*% (weights * x)
  theta = mu_mat + eta_centered %*% H

  loglike = log_like_Bernoulli_weighted(q = q, theta = theta, weights)
  min_loss = -loglike / sum_weights
  best_HU = HU
  best_loglike = loglike
  if (!quiet) {
    cat(0,"  ", min_loss, "\n")
  }

  loss_trace <- proj_loss_trace <- numeric(max_iters + 1)
  loss_trace[1] <- proj_loss_trace[1] <- min_loss

  H_lag = H
  for (i in 1:max_iters) {
    y = H + (i - 2) / (i + 1) * (H - H_lag)
    # y = H
    H_lag = H
    # y = H
    step = ss_factor / L

    Phat = inv.logit.mat(mu_mat + eta_centered %*% y)
    Phat[q == 0] <- 0
    etatP = t(eta_centered) %*% (Phat * weights)
    deriv = etatX - etatP
    deriv = deriv + t(deriv) - diag(diag(deriv))

    H = y + step * deriv
    HU = project.Fantope(H, k, partial_decomp = partial_decomp)
    H = HU$H

    theta = mu_mat + eta_centered %*% H
    loglike = log_like_Bernoulli_weighted(q = q, theta = theta, weights)
    loss_trace[i + 1] = -loglike / sum_weights

    proj_theta = mu_mat + eta_centered %*% tcrossprod(HU$U)
    proj_loglike = log_like_Bernoulli_weighted(q = q, theta = proj_theta, weights)
    proj_loss_trace[i + 1] = -proj_loglike / sum_weights

    if (!quiet) {
      cat(i,"  ",loss_trace[i + 1],"  ",proj_loss_trace[i + 1],"\n")
    }
    if (loss_trace[i + 1] < min_loss) {
      min_loss = loss_trace[i + 1]
      best_HU = HU
      best_loglike = loglike
    }
    if (abs(loss_trace[i + 1]-loss_trace[i]) < conv_criteria | min_loss == 0) {
      break
    }
  }

  # calculate the null log likelihood for % deviance explained
  # assumes no missing data
  # if (main_effects) {
  #   null_proportions = x_bar
  # } else {
  #   null_proportions = rep(0.5, d)
  # }
  # null_loglikes <- null_proportions * log(null_proportions) +
  #   (1 - null_proportions) * log(1 - null_proportions)
  # null_loglike = sum((null_loglikes * colSums(q!=0))[!(null_proportions %in% c(0, 1))])
  null_loglike = log_like_Bernoulli_weighted(q = q, theta = mu_mat, weights)

  object = list(mu = mu,
                H = best_HU$H,
                U = best_HU$U,
                PCs = eta_centered %*% best_HU$U,
                m = m,
                M = m, # need to depricate after 0.1.1
                iters = i,
                loss_trace = loss_trace[1:(i + 1)],
                proj_loss_trace = proj_loss_trace[1:(i + 1)],
                prop_deviance_expl = 1 - best_loglike / null_loglike,
                rank = best_HU$rank)
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
#' @param partial_decomp logical; if \code{TRUE}, the function uses the RSpectra package
#'   to quickly calculate the eigendecomposition when \code{ncol(x)} is large and \code{k} is small
#'
#' @return
#' \item{H}{a rank \code{k} Fantope matrix}
#' \item{U}{a \code{k}-dimentional orthonormal matrix with the first \code{k} eigenvectors of \code{H}}
#' \item{rank}{the rank of the Fantope matrix \code{H}}
#' @export
project.Fantope <- function(x, k, partial_decomp = FALSE) {
  if (partial_decomp) {
    if (!requireNamespace("RSpectra", quietly = TRUE)) {
      message("RSpectra must be installed to use partial_decomp")
      partial_decomp = FALSE
    }
  }
  
  d = ncol(x)
  
  if (partial_decomp) {
    eig = eigen(x, symmetric = TRUE, only.values = TRUE)
  } else {
    eig = eigen(x, symmetric = TRUE)
  }
  vals = eig$values
  
  if (sum(pmin(1, eig$values)) < k) {
    warning("k is larger than the rank of x")
  }
  lower = vals[length(vals)] - k / length(vals)
  upper = max(vals)
  while(TRUE) {
    theta = (lower + upper) / 2
    sum.eig.vals = sum(pmin(pmax(vals - theta, 0), 1))
    if (abs(sum.eig.vals - k) < 1e-10) {
      break
    } else if (sum.eig.vals > k) {
      lower = theta
    } else {
      upper = theta
    }
  }
  vals = pmin(pmax(vals - theta, 0), 1)
  num_vals = sum(vals > .Machine$double.eps * 10)
  
  if (partial_decomp & num_vals < d) {
    # do a partial decomposition of between k and number of non-zero values (plus a buffer)
    eig = RSpectra::eigs_sym(x, k = max(min(num_vals + 1, d), ceiling(k)))
  }
  # if RSpectra gave a bad result or the number of non-zero e-vals = d
  if (!partial_decomp || num_vals == d || any(eig$values[1:num_vals] < 0)) {
    if (partial_decomp) {
      eig = eigen(x, symmetric = TRUE)
    }
    H = eig$vectors %*% diag(vals) %*% t(eig$vectors)
  } else {
    H = eig$vectors[, 1:num_vals, drop = FALSE] %*% 
      diag(vals[1:num_vals], num_vals, num_vals) %*% 
      t(eig$vectors[, 1:num_vals, drop = FALSE])
  }
  
  return(list(H = H,
              U = matrix(eig$vectors[, 1:ceiling(k)], nrow(x), ceiling(k)),
              rank = num_vals))
}

#' @title Predict Convex Logistic PCA scores or reconstruction on new data
#'
#' @description Predict Convex Logistic PCA scores or reconstruction on new data
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
#' clpca = convexLogisticPCA(mat, k = 1, m = 4, main_effects = FALSE)
#'
#' PCs = predict(clpca, mat_new)
#' @export
predict.clpca <- function(object, newdata, type = c("PCs", "link", "response"), ...) {
  type = match.arg(type)

  if (type == "PCs") {
    if (missing(newdata)) {
      PCs = object$PCs
    } else {
      eta = ((as.matrix(newdata) * 2) - 1) * object$m
      eta_centered = scale(eta, center = object$mu, scale = FALSE)
      eta_centered[is.na(newdata)] <- 0
      PCs = eta_centered %*% object$U
    }
    return(PCs)
  } else {
    eta = ((as.matrix(newdata) * 2) - 1) * object$m
    eta_centered = scale(eta, center = object$mu, scale = FALSE)
    eta_centered[is.na(newdata)] <- 0
    theta = outer(rep(1, nrow(eta)), object$mu) + eta_centered %*% object$H
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
#' @param x convex logistic PCA object
#' @param type the type of plot \code{type = "trace"} plots the algorithms progress by
#' iteration, \code{type = "loadings"} plots the first 2 PC loadings,
#' \code{type = "scores"} plots the first 2 PC scores
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
#' clpca = convexLogisticPCA(mat, k = 2, m = 4, main_effects = FALSE)
#'
#' \dontrun{
#' plot(clpca)
#' }
#' @export
plot.clpca <- function(x, type = c("trace", "loadings", "scores"), ...) {
  type = match.arg(type)

  if (type == "trace") {
    df = data.frame(Iteration = 0:x$iters,
                    NegativeLogLikelihood = x$loss_trace)
    p <- ggplot2::ggplot(df, ggplot2::aes_string("Iteration", "NegativeLogLikelihood")) +
      ggplot2::geom_line()
  } else if (type == "loadings") {
    df = data.frame(x$U)
    colnames(df) <- paste0("PC", 1:ncol(df))
    if (ncol(df) == 1) {
      df$PC2 = 0
      p <- ggplot2::ggplot(df, ggplot2::aes_string("PC1", "PC2")) + ggplot2::geom_point() + 
        ggplot2::labs(y = NULL)
    } else {
      p <- ggplot2::ggplot(df, ggplot2::aes_string("PC1", "PC2")) + ggplot2::geom_point()
    }
  } else if (type == "scores") {
    df = data.frame(x$PCs)
    colnames(df) <- paste0("PC", 1:ncol(df))
    if (ncol(df) == 1) {
      df$PC2 = 0
      p <- ggplot2::ggplot(df, ggplot2::aes_string("PC1", "PC2")) + ggplot2::geom_point() + 
        ggplot2::labs(y = NULL)
    } else {
      p <- ggplot2::ggplot(df, ggplot2::aes_string("PC1", "PC2")) + ggplot2::geom_point()
    }
  }

  return(p)
}

#' @export
print.clpca <- function(x, ...) {
  cat(nrow(x$PCs), "rows and ")
  cat(nrow(x$H), "columns\n")
  cat("Rank", ncol(x$U), "Fantope solution with m =", x$m, "\n")
  cat("\n")
  cat(round(x$prop_deviance_expl * 100, 1), "% of deviance explained\n", sep = "")
  cat(x$iters, "iterations to converge\n")

  invisible(x)
}

#' @title CV for convex logistic PCA
#'
#' @description
#' Run cross validation on dimension and \code{m} for convex logistic PCA
#'
#' @param x matrix with all binary entries
#' @param ks the different dimensions \code{k} to try
#' @param ms the different approximations to the saturated model \code{m} to try
#' @param folds if \code{folds} is a scalar, then it is the number of folds. If
#'  it is a vector, it should be the same length as the number of rows in \code{x}
#' @param quiet logical; whether the function should display progress
#' @param Ms depricated. Use \code{ms} instead
#' @param ... Additional arguments passed to convexLogisticPCA
#'
#' @return A matrix of the CV negative log likelihood with \code{k} in rows and
#'  \code{m} in columns
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
#' negloglikes = cv.clpca(mat, ks = 1:9, ms = 3:6)
#' plot(negloglikes)
#' }
#' @export
cv.clpca <- function(x, ks, ms = seq(2, 10, by = 2), folds = 5, quiet = TRUE, Ms, ...) {
  if (!missing(Ms)) {
    ms = Ms
    warning("Ms is depricated. Use ms instead.\n", 
            "Using ms in ", paste(ms, collapse = ","))
  }
  # TODO: does not support weights
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

  log_likes = matrix(0, length(ks), length(ms),
                     dimnames = list(k = ks, m = ms))
  for (k in ks) {
    for (m in ms) {
      if (!quiet) {
        cat("k =", k, "m =", m, "")
      }
      for (c in unique(cv)) {
        if (!quiet) {
          cat(".")
        }
        clpca = convexLogisticPCA(x[c != cv, ], k = k, m = m, ...)
        pred_theta = predict(clpca, newdat = x[c == cv, ], type = "link")
        log_likes[k == ks, m == ms] = log_likes[k == ks, m == ms] +
          log_like_Bernoulli(q = q[c == cv, ], theta = pred_theta)
      }
      if (!quiet) {
        cat("", -log_likes[k == ks, m == ms], "\n")
      }
    }
  }
  class(log_likes) <- c("matrix", "cv.lpca")
  which_min = which(log_likes == max(log_likes), arr.ind = TRUE)
  if (!quiet) {
    cat("Best: k =", ks[which_min[1]], "m =", ms[which_min[2]], "\n")
  }

  return(-log_likes)
}
