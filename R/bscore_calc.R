#-----------------------------------------------------------------------------#
#' Compute psi_xp(B_i, X_i, theta)    
#'
#' @param predictors character vector of predictors from model
#' @param B character naming B variable in data
#' @param G character vector (length == 1) naming group variable in data
#' @param theta
#' @param data dataframe
#' @return N X length(theta) matrix of scores
#' @export
#-----------------------------------------------------------------------------#

bscore_calc <- function(predictors, B, G, theta, data){
  N <- length(unique(data[, G]))
  
  if(length(theta) != (length(predictors) + 2)){
    stop("The length of theta is not equal to the number of predictors + 2 ")
  }
  
  out <- matrix(nrow = N, ncol = length(theta))  
  for(ii in 1:N){
    grp <- data[data$group == ii, ]
    X <- as.matrix(cbind(1, grp[, predictors]))
    out[ii, ] <- cmp_scores_ll(theta = theta, B = grp[ , B], X = X)
  }
  
  return(out)
}