# Logistic PCA

[![Build Status](https://travis-ci.org/andland/logisticPCA.png?branch=master)](https://travis-ci.org/andland/logisticPCA) [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/logisticPCA)](http://cran.r-project.org/web/packages/logisticPCA)

`logisticPCA` is an R package for dimensionality reduction of binary data. Please note that it is still in the very early stages of development and the conventions will possibly change in the future. A manuscript describing logistic PCA can be found [here](http://www.stat.osu.edu/~yklee/mss/tr890.pdf).

## Installation

To install R, visit [r-project.org/](http://www.r-project.org/).

The package can be installed by downloading from CRAN.
```R
install.packages("logisticPCA")
```

To install the development version, first install `devtools` from CRAN. Then run the following commands.
```R
# install.packages("devtools")
library("devtools")
install_github("andland/logisticPCA")
```

## Classes
Three types of dimensionality reduction are given. For all the functions, the user must supply the desired dimension `k`. The data must be an `n x d` matrix comprised of binary variables (i.e. all `0`'s and `1`'s).

### Logistic PCA
`logisticPCA()` estimates the natural parameters of a Bernoulli distribution in a lower dimensional space. This is done by projecting the natural parameters from the saturated model. A rank-`k` projection matrix, or equivalently a `d x k` orthogonal matrix `U`, is solved for to minimize the Bernoulli deviance. Since the natural parameters from the saturated model are either negative or positive infinity, an additional tuning parameter `M` is needed to approximate them. You can use `cv.lpca()` to select `M` by cross validation. Typical values are in the range of `3` to `10`. 

`mu` is a main effects vector of length `d` and `U` is the `d x k` loadings matrix.

### Logistic SVD
`logisticSVD()` estimates the natural parameters by a matrix factorization. `mu` is a main effects vector of length `d`, `B` is the `d x k` loadings matrix, and `A` is the `n x k` principal component score matrix.

### Convex Logistic PCA
`convexLogisticPCA()` relaxes the problem of solving for a projection matrix to solving for a matrix in the `k`-dimensional Fantope, which is the convex hull of rank-`k` projection matrices. This has the advantage that the global minimum can be obtained efficiently. The disadvantage is that the `k`-dimensional Fantope solution may have a rank much larger than `k`, which reduces interpretability. It is also necessary to specify `M` in this function.

`mu` is a main effects vector of length `d`, `H` is the `d x d` Fantope matrix, and `U` is the `d x k` loadings matrix, which are the first `k` eigenvectors of `H`.

## Methods
Each of the classes has associated methods to make data analysis easier.

* `print()`: Prints a summary of the fitted model.
* `fitted()`: Fits the low dimensional matrix of either natural parameters or probabilities.
* `predict()`: Predicts the PCs on new data. Can also predict the low dimensional matrix of natural parameters or probabilities on new data.
* `plot()`: Either plots the deviance trace, the first two PC loadings, or the first two PC scores using the package `ggplot2`.

In addition, there are functions for performing cross validation.

* `cv.lpca()`, `cv.lsvd()`, `cv.clpca()`: Run cross validation over the rows of the matrix to assess the fit of `M` and/or `k`.
* `plot.cv()`: Plots the results of the `cv()` method.
