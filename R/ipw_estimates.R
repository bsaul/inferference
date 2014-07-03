#' Create IPW estimates 
#'  
#'  
#' @param y vector of unweighted group means (i.e. output from \code{\link{group_means}}
#' @param a treatment level (0,1) for which to compute weights. Defaults to NULL which returns marginal.
#' @param weights weight matrix to use
#' @param rescale.factor factor by which to rescale values
#' @param na.rm exclude groups missing weights? defaults to FALSE.
#' @return length(alpha) vector of IPW estimates
#' @export

ipw_estimates <- function(y, G, A, B, data, weights, weight_dervs, predictors,
                          rescale.factor, na.rm = FALSE){
  
  # Point Estimates
  points <- ipw_point_estimates(y = y, A = A, G = G, B = B, data = data,
                                weights = weights, rescale.factor = rescale.factor)

  # Parts to Variance Estimates
  bscores <- Bscore(predictors, B, G, theta, data)
  #V       <- V_matrices(bscores, points)
  Upart   <- Upartial(y = y, A = A, G = G, B = B, data = data,
                      weights = weights, weight_dervs = weight_dervs, 
                      rescale.factor = rescale.factor, na.rm = na.rm)

  out <- list(point_estimates = points, bscores = bscores, Upart = Upart)
  return(out)
}

