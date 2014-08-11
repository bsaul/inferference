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

bscore_calc <- function(predictors, 
                        B, 
                        G, 
                        theta, 
                        data,
                        set.NA.to.0 = TRUE){

  if(length(theta) != (length(predictors) + 2)){
    stop("The length of theta is not equal to the number of predictors + 2 ")
  }
  
  dd <- data[order(as.numeric(data[, G])), ]  
  X <- cbind(1, dd[, c(predictors, B)])
  p <- length(theta) - 1
  
  b.list <- by(X, INDICES = dd[, G], simplify = TRUE, 
               FUN = function(x) {
               x <- as.matrix(x) # PrAX expects a matrix
               cmp_scores_ll(theta = theta, 
                             B = x[ , (p + 1)], 
                             X = x[ , 1:p]) })

  out <- matrix(unlist(b.list), ncol = p + 1, byrow = TRUE,
                     dimnames = list(sort(unique(data[ , G])), 
                                     names(theta)))
  ## replace any Bscores with 0 ##
  if(set.NA.to.0 == TRUE) {
    out[is.na(out)] <- 0
  }
  
  return(out)
}