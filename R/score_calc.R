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
                        groups, 
                        theta, 
                        data,
                        set.NA.to.0 = TRUE){
  ## Warnings ##
  if(length(theta) != (length(predictors) + 2)){
    stop("The length of theta is not equal to the number of predictors + 2 ")
  }
  
  ## Necessary bits ##
  G  <- data[, groups]
  X  <- cbind(1, data[, c(predictors, B)])
  p  <- ncol(X) - 1
  gg <- unique(G)
  
  ## Compute score for each group and parameter ##
  b.list <- by(X, INDICES = G, simplify = TRUE, 
               FUN = function(x) {
               x <- as.matrix(x) # PrAX expects a matrix
               cmp_scores_ll(theta = theta, 
                             B = x[ , (p + 1)], 
                             X = x[ , 1:p]) })
  
  ## Reshape list into matrix ##
  out <- matrix(unlist(b.list, use.names = FALSE), 
                ncol = p + 1, 
                byrow = TRUE,
                dimnames = list(gg, names(theta)))
  
  ## replace any Bscores with 0 ##
  if(set.NA.to.0 == TRUE) {
    out[is.na(out)] <- 0
  }
  
  return(out)
}