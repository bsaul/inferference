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

ipw_estimates <- function(y, G, A, data, weights, weight_dervs,
                          predictors, rescale.factor, set.NA.to.0 = TRUE){
  
  # Point Estimates
  points <- ipw_point_estimates(y = y, A = A, G = G, data = data,
                                weights = weights, 
                                rescale.factor = rescale.factor,
                                set.NA.to.0 = set.NA.to.0)

  # Parts to Variance Estimates
  Upart   <- ipw_point_estimates(y = y, A = A, G = G, data = data,
                                 weights = weight_dervs, 
                                 rescale.factor = rescale.factor, 
                                 set.NA.to.0 = set.NA.to.0)

  out <- list(point_estimates = points, 
              Upart = Upart)
  return(out)
}

