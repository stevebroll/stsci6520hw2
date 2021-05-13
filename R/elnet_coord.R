#' Title
#'
#' @param x Predictor matrix x, with n rows and p columns
#' @param y Response vector y, with n elements
#' @param alpha Elastic net alpha parameter. elastic net simplifies to Ridge
#'  regression under alpha = 0, and to Lasso regression under alpha = 1
#' @param lambda Lasso tuning parameter
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
#' y = x %*% beta + rnorm(20, 0, 1)
#'
elnet_coord <- function(x, y, alpha, lambda = .1, tol = .0001){
  if(alpha < 0 || alpha > 1){
    stop("alpha should be in [0,1]")
  }
  x = as.matrix(scale(x))
  p = ncol(x)
  la = lambda*(1- alpha) + 1
  betahat = rep(0,p)
  maxdiff = 1
  while(maxdiff > tol){
    old_betahat = betahat
    for(j in 1:p){
      r = y - x[,-j] %*% betahat[-j]
      rx = 1/nrow(x) * t(r) %*% x[,j]
      betahat[j] = soft(rx, lambda, alpha)/la
    }
    maxdiff = max(abs(betahat - old_betahat))
  }
  betahat = data.frame(Beta = betahat)
  #cat(paste("Betas for alpha = ", alpha, " and lambda = ", lambda, ": \n"))
  return(betahat)
}
