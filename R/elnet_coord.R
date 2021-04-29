#' Title
#'
#' @param x Predictor matrix x, with n rows and p columns
#' @param y Response vector y, with n elements
#' @param lambda Lasso tuning parameter
#' @param alpha Elastic net alpha parameter. elastic net simplifies to Ridge
#'  regression under alpha = 0, and to Lasso regression under alpha = 1
#' @param tol Tolerance for convergence of betas. Algorithm will complete when
#' the maximum change in beta values is smaller than the tol value
#'
#' @return Matrix where each column is the p-length vector of beta estimates
#' for one iteration
#' @export
#'
#' @examples
#' beta = c(2,0,-2,0,1,0,-1,rep(0,13))
#' covar = diag(1,20)
#' covar[1,2] = covar[2,1] = covar[5,6] = covar[6,5] = .8
#' x = MASS::mvrnorm(20, rep(0, 20), covar)
#' y = x %*% beta + rnorm(20s, 0, 1)
#'
elnet_coord <- function(x, y, lambda, alpha, tol){
  x = scale(x)
  la = lambda*(1- alpha) + 1
  betahat = rep(0,20)
  betamat = matrix(0, 20, 1)
  maxdiff = 1
  i = 2
  while(maxdiff > tol){
    for(j in 1:ncol(x)){
      r = y - x[,-j] %*% betahat[-j]
      rx = 1/nrow(x) * t(r) %*% x[,j]
      betahat[j] = soft(rx, lambda, alpha)/la
    }
    betamat = cbind(betamat, betahat)
    maxdiff = max(abs(betamat[, i] - betamat[, i-1]))
    i = i + 1
  }
  return(betamat)
}
