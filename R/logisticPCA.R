#' @title Logistic Principal Component Analysis
#'
#' @description
#' Dimensionality reduction for binary data by extending Pearson's
#' PCA formulation to minimize Binomial deviance
#'
#' @param x matrix with all binary entries
#' @param k number of principal components to return
#' @param m value to approximate the saturated model. If \code{m = 0}, m is solved for
#' @param quiet logical; whether the calculation should give feedback
#' @param use_irlba logical; if \code{TRUE}, the function uses the irlba package
#'   to more quickly calculate the eigen-decomposition. This is usually faster when 
#'   \code{ncol(x) > 100} and \code{k} is small
#' @param max_iters number of maximum iterations
#' @param conv_criteria convergence criteria. The difference between average deviance
#'   in successive iterations
#' @param random_start logical; whether to randomly inititalize the parameters. If \code{FALSE},
#'   function will use an eigen-decomposition as starting value
#' @param start_U starting value for the orthogonal matrix
#' @param start_mu starting value for mu. Only used if \code{main_effects = TRUE}
#' @param main_effects logical; whether to include main effects in the model
#' @param validation optional validation matrix. If supplied and \code{m = 0}, the
#'   validation data is used to solve for \code{m}
#' @param M depricated. Use \code{m} instead
#'
#' @return An S3 object of class \code{lpca} which is a list with the
#' following components:
#' \item{mu}{the main effects}
#' \item{U}{a \code{k}-dimentional orthonormal matrix with the loadings}
#' \item{PCs}{the princial component scores}
#' \item{m}{the parameter inputed or solved for}
#' \item{iters}{number of iterations required for convergence}
#' \item{loss_trace}{the trace of the average negative log likelihood of the algorithm.
#'    Should be non-increasing}
#' \item{prop_deviance_expl}{the proportion of deviance explained by this model.
#'    If \code{main_effects = TRUE}, the null model is just the main effects, otherwise
#'    the null model estimates 0 for all natural parameters.}
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
#' # run logistic PCA on it
#' lpca = logisticPCA(mat, k = 1, m = 4, main_effects = FALSE)
#'
#' # Logistic PCA likely does a better job finding latent features
#' # than standard PCA
#' plot(svd(mat_logit)$u[, 1], lpca$PCs[, 1])
#' plot(svd(mat_logit)$u[, 1], svd(mat)$u[, 1])
#' @export
logisticPCA <- function(x, k = 2, m = 4, quiet = TRUE, use_irlba = FALSE,
                        max_iters = 1000, conv_criteria = 1e-5, random_start = FALSE,
                        start_U, start_mu, main_effects = TRUE, validation, M) {
  if (!missing(M)) {
    m = M
    warning("M is depricated. Use m instead. ",
            "Using m = ", m)
  }
  use_irlba = use_irlba && requireNamespace("irlba", quietly = TRUE)
  q = as.matrix(2 * x - 1)
  missing_mat = is.na(q)
  q[is.na(q)] <- 0 # forces Z to be equal to theta when data is missing
  n = nrow(q)
  d = ncol(q)
  if (m == 0) {
    m = 4
    solve_M = TRUE
    if (!missing(validation)) {
      if (ncol(validation) != ncol(x)) {
        stop("validation does not have the same variables as x")
      }
      validation = as.matrix(validation)
      q_val = 2 * validation - 1
      q_val[is.na(q_val)] <- 0
    }
  } else {
    solve_M = FALSE
  }

  if (main_effects) {
    if (!missing(start_mu)) {
      mu = start_mu
    } else {
      mu = colMeans(m * q)
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
  eta = m * q + missing_mat * outer(rep(1, n), mu)
  theta = outer(rep(1, n), mu) + scale(eta, center = mu, scale = FALSE) %*% tcrossprod(U)
  loglike <- log_like_Bernoulli(q = q, theta = theta)
  loss_trace[1] = (-loglike) / sum(q!=0)
  ptm <- proc.time()

  if (!quiet) {
    cat(0, "  ", loss_trace[1], "")
    cat("0 hours elapsed\n")
  }

  for (i in 1:max_iters) {
    last_U = U
    last_m = m
    last_mu = mu

    if (solve_M) {
      if (missing(validation)) {
        Phat = inv.logit.mat(theta)
        M_slope = sum(((Phat - x) * (q %*% tcrossprod(U)))[q != 0])
        M_curve = sum((Phat * (1 - Phat) * (q %*% tcrossprod(U))^2)[q != 0])
      } else {
        lpca_obj = structure(list(mu = mu, U = U, m = m),
                             class = "lpca")
        Phat = predict(lpca_obj, newdata = validation, type = "response")
        M_slope = sum(((Phat - validation) * (q_val %*% tcrossprod(U)))[q_val != 0])
        M_curve = sum((Phat * (1 - Phat) * (q_val %*% tcrossprod(U))^2)[q_val != 0])
      }
      m = max(m - M_slope / M_curve, 0)

      eta = m * q + missing_mat * outer(rep(1, n), mu)
      theta = outer(rep(1, n), mu) + scale(eta, center = mu, scale = FALSE) %*% tcrossprod(U)
    }

    Z = as.matrix(theta + 4 * q * (1 - inv.logit.mat(q * theta)))
    if (main_effects) {
      mu = as.numeric(colMeans(Z - eta %*% tcrossprod(U)))
    }

    eta = m * q + missing_mat * outer(rep(1, n), mu)

    mat_temp = crossprod(scale(eta, center = mu, scale = FALSE), Z)
    mat_temp = mat_temp + t(mat_temp) - crossprod(eta) + n * outer(mu, mu)

    # irlba sometimes gives poor estimates of e-vectors
    # so I switch to standard eigen if it does
    repeat {
      if (use_irlba) {
        if (packageVersion("irlba") < "2.0.0") {
          udv = irlba::irlba(mat_temp, nu=k, nv=k)
          U = matrix(udv$u[, 1:k], d, k)
        } else {
          eig = irlba::partial_eigen(mat_temp, n = k)
          U = matrix(eig$vectors[, 1:k], d, k)
        }
      } else {
        eig = eigen(mat_temp, symmetric=TRUE)
        U = matrix(eig$vectors[, 1:k], d, k)
      }

      theta = outer(rep(1, n), mu) + scale(eta, center = mu, scale = FALSE) %*% tcrossprod(U)
      this_loglike <- log_like_Bernoulli(q = q, theta = theta)

      if (!use_irlba | this_loglike>=loglike) {
        loglike = this_loglike
        break
      } else {
        use_irlba = FALSE
        warning("irlba::partial_eigen was too inaccurate. Switched to base::eigen")
      }
    }

    loss_trace[i + 1] = (-loglike) / sum(q!=0)

    if (!quiet) {
      time_elapsed = as.numeric(proc.time() - ptm)[3]
      tot_time = max_iters / i * time_elapsed
      time_remain = tot_time - time_elapsed
      cat(i, "  ", loss_trace[i + 1], "")
      cat(round(time_elapsed / 3600, 1), "hours elapsed. Max", round(time_remain / 3600, 1), "hours remain.\n")
    }
    if (i > 4) {
      # when solving for m, the monoticity does not apply
      if (solve_M) {
        if (abs(loss_trace[i] - loss_trace[i + 1]) < conv_criteria) {
          break
        }
      } else {
        if ((loss_trace[i] - loss_trace[i + 1]) < conv_criteria) {
          break
        }
      }
    }
  }

  # test if loss function increases
  if ((loss_trace[i + 1] - loss_trace[i]) > (1e-10)) {
    U = last_U
    mu = last_mu
    m = last_m
    i = i - 1

    if (!solve_M) {
      warning("Algorithm stopped because deviance increased.\nThis should not happen!")
    }
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

  eta = m * q + missing_mat * outer(rep(1, n), mu)

  object <- list(mu = mu,
                 U = U,
                 PCs = scale(eta, center = mu, scale = FALSE) %*% U,
                 m = m,
                 M = m, # need to depricate after 0.1.1
                 iters = i,
                 loss_trace = loss_trace[1:(i + 1)],
                 prop_deviance_expl = 1 - loglike / null_loglike)
  class(object) <- "lpca"
  object
}


#' @title Predict Logistic PCA scores or reconstruction on new data
#'
#' @description Predict Logistic PCA scores or reconstruction on new data
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
#' lpca = logisticPCA(mat, k = 1, m = 4, main_effects = FALSE)
#'
#' PCs = predict(lpca, mat_new)
#' @export
predict.lpca <- function(object, newdata, type = c("PCs", "link", "response"), ...) {
  type = match.arg(type)

  if (missing(newdata)) {
    PCs = object$PCs
  } else {
    q = as.matrix(newdata) * 2 - 1
    q[is.na(q)] <- 0
    eta = object$m * q + is.na(q) * outer(rep(1, nrow(newdata)), object$mu)
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
#' lpca = logisticPCA(mat, k = 1, m = 4, main_effects = FALSE)
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
#' @param x logistic PCA object
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
#' # run logistic PCA on it
#' lpca = logisticPCA(mat, k = 2, m = 4, main_effects = FALSE)
#'
#' \dontrun{
#' plot(lpca)
#' }
#' @export
plot.lpca <- function(x, type = c("trace", "loadings", "scores"), ...) {
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
print.lpca <- function(x, ...) {
  cat(nrow(x$PCs), "rows and ")
  cat(nrow(x$U), "columns\n")
  cat("Rank", ncol(x$U), "solution with m =", x$m, "\n")
  cat("\n")
  cat(round(x$prop_deviance_expl * 100, 1), "% of deviance explained\n", sep = "")
  cat(x$iters, "iterations to converge\n")

  invisible(x)
}

#' @title CV for logistic PCA
#'
#' @description
#' Run cross validation on dimension and \code{m} for logistic PCA
#'
#' @param x matrix with all binary entries
#' @param ks the different dimensions \code{k} to try
#' @param ms the different approximations to the saturated model \code{m} to try
#' @param folds if \code{folds} is a scalar, then it is the number of folds. If
#'  it is a vector, it should be the same length as the number of rows in \code{x}
#' @param quiet logical; whether the function should display progress
#' @param Ms depricated. Use \code{ms} instead
#' @param ... Additional arguments passed to \code{logisticPCA}
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
#' negloglikes = cv.lpca(mat, ks = 1:9, ms = 3:6)
#' plot(negloglikes)
#' }
#' @export
cv.lpca <- function(x, ks, ms = seq(2, 10, by = 2), folds = 5, quiet = TRUE, Ms, ...) {
  if (!missing(Ms)) {
    ms = Ms
    warning("Ms is depricated. Use ms instead.\n", 
            "Using ms in ", paste(ms, collapse = ","))
  }
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
        lpca = logisticPCA(x[c != cv, ], k = k, m = m, ...)
        pred_theta = predict(lpca, newdat = x[c == cv, ], type = "link")
        log_likes[k == ks, m == ms] = log_likes[k == ks, m == ms] +
          log_like_Bernoulli(q = q[c == cv, ], theta = pred_theta)
        #         log_likes[k == ks, m == ms] = log_likes[k == ks, m == ms] +
        #           sum(log(inv.logit.mat(q[c == cv, ] * pred_theta)))
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

#' @title Plot CV for logistic PCA
#'
#' @description
#' Plot cross validation results logistic PCA
#'
#' @param x a \code{cv.lpca} object
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
#' negloglikes = cv.lpca(dat, ks = 1:9, ms = 3:6)
#' plot(negloglikes)
#' }
#' @export
plot.cv.lpca <- function(x, ...) {
  # replaces reshape2::melt(-x, value.name = "NegLogLikelihood")
  ms = type.convert(colnames(x))
  ks = type.convert(rownames(x))
  df = data.frame(k = rep(ks, times = length(ms)),
                  m = rep(ms, each = length(ks)),
                  NegLogLikelihood = as.vector(x))
  
  if (ncol(x) == 1) {
    df$m = factor(df$m)
    p <- ggplot2::ggplot(df, ggplot2::aes_string("k", "NegLogLikelihood", colour = "m")) +
      ggplot2::geom_line()
  } else {
    df$k = factor(df$k)
    p <- ggplot2::ggplot(df, ggplot2::aes_string("m", "NegLogLikelihood", colour = "k")) +
      ggplot2::geom_line()
  }
  return(p)
}
