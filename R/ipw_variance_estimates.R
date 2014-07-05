#-----------------------------------------------------------------------------#
# Compute psi_xp(B_i, X_i, theta)    
#
# @param predictors character vector of predictors from model
# @param B character naming B variable in data
# @param G character vector (length == 1) naming group variable in data
# @param theta
# @param data dataframe
# @return N X length(theta) matrix of scores
#-----------------------------------------------------------------------------#

Bscore <- function(predictors, B, G, theta, data){
  N <- length(unique(data[, G]))
  
  out <- matrix(nrow = N, ncol = length(theta))  
  for(ii in 1:N){
    grp <- data[data$group == ii, ]
    X <- as.matrix(cbind(1, grp[, predictors]))
    out[ii, ] <- cmp_scores_ll(theta = theta, B = grp[ , B], X = X)
  }
  
  return(out)
}

#-----------------------------------------------------------------------------#
#' V Matrix
#'  
#' @param predictors
#' @return V matrix
#' @export
#' 
#-----------------------------------------------------------------------------#
V_matrix_effect <- function(Bscores, ipw_obj, alpha1, alpha2, trt.lvl1, trt.lvl2){

  N <- ipw_obj$summary$ngroups
  p <- dim(Bscores)[2]
  
  tmp <- array(dim = c(p+1, p+1, N))

  for(ii in 1:N){
    hold <- c(Bscores[ii, ], ipw_obj$ipw$contrasts$group_resid[alpha1, trt.lvl1, alpha2, trt.lvl2, ii, alpha1, trt.lvl1, alpha2, trt.lvl2])
    tmp[ , , ii] <- hold %*% t(hold)
  }
  
  V <- apply(tmp, 1:2, sum, na.rm = T)/N

  return(V)
}


