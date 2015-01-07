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

#' @title Logistic Principal Component Analysis
#' 
#' @description 
#' Dimension reduction for binary data by extending Pearson's
#' PCA formulation to minimize Binomial deviance
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
#' @param start_U starting value for the orthogonal matrix
#' @param start_mu starting value for mu. Only used if \code{main_effects = TRUE}
#' @param main_effects logical; whether to include main effects in the model
#' 
#' @return An S3 object of class \code{lpca} which is a list with the
#' following components:
#' \item{mu}{the main effects}
#' \item{U}{a \code{k}-dimentional orthonormal matrix with the loadings}
#' \item{PCs}{the princial components}
#' \item{M}{the same parameter as inputed}
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
  # better name for k and dat.
  use_irlba = use_irlba && requireNamespace("irlba", quietly = TRUE)
  q = as.matrix(2 * x - 1)
  q[is.na(q)] <- 0 # forces Z to be equal to theta when data is missing
  n = nrow(q)
  d = ncol(q)
  eta = q * abs(M)
  
  if (main_effects) {
    if (!missing(start_mu)) {
      mu = start_mu
    } else {
      mu = colMeans(eta)
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
  
  etaTeta = crossprod(eta)
  
  loss_trace = numeric(max_iters + 1)
  theta = outer(rep(1, n), mu) + scale(eta, center = mu, scale = FALSE) %*% U %*% t(U)
  loglike <- .Call(compute_loglik, q, theta)
  # loglike = sum(x * theta) - sum(pmax(0, theta))
  loss_trace[1] = (-loglike) / sum(q!=0)
  ptm <- proc.time()
  
  if (!quiet) {
    cat(0, "  ", loss_trace[1], "")
    cat("0 hours elapsed\n")
  }
  
  for (m in 1:max_iters) {
    last_U = U
    last_mu = mu
    
    Z = as.matrix(theta + 4 * q * (1 - inv.logit.mat(q * theta)))
    if (main_effects) {
      mu = as.numeric(colMeans(Z - eta %*% U %*% t(U)))
    }
    
    mat_temp = t(scale(eta, center = mu, scale = FALSE)) %*% Z
    mat_temp = mat_temp + t(mat_temp) - etaTeta + n * outer(mu, mu)
    
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
      
      theta = outer(rep(1, n), mu) + scale(eta, center = mu, scale = FALSE) %*% tcrossprod(U)
      this_loglike <- .Call(compute_loglik, q, theta)
      # this_loglike = sum(x * theta) - sum(pmax(0, theta))
      
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
      if ((loss_trace[m] - loss_trace[m+1]) < conv_criteria)
        break
    }
  }
  
  # test if loss function increases
  if ((loss_trace[m + 1] - loss_trace[m]) > (1e-10)) {
    U = last_U
    mu = last_mu
    m = m - 1
    
    theta = outer(rep(1, n), mu) + scale(eta, center = mu, scale = FALSE) %*% U %*% t(U)
    loglike = sum(log(inv.logit.mat(q * theta))[q!=0])
    warning("Algorithm stopped because deviance increased.\nThis should not happen!")
  }
  
  # calculate the null log likelihood for % deviance explained
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
                 PCs = scale(eta, center = mu, scale = FALSE) %*% U,
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
#' # plot(lpca)
#' @export
plot.lpca <- function(object, type = c("trace", "loadings"), ...) {
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