#' Soft thresholding for elnet_coord
#'
#' @param rx rx object from elnet_coord function, product of partial residual
#' and a column of x
#' @param lambda lambda value specified in elnet_coord
#' @param alpha alpha value specified in elnet_coord
#'
#' @return Soft thresholded rx value
#'
#'
soft <- function(rx, lambda, alpha){
  return(sign(rx)* max(0, (abs(rx) - lambda * alpha /2)))
}
