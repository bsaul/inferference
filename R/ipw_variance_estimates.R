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
# Compute V11 
#
# @param 
# @return 
#-----------------------------------------------------------------------------#

V11 <- function(bscores, na.rm = FALSE){  

  if(na.rm == FALSE) {
    # Replace missing values with 0
    bscores[is.na(bscores)] <- 0
  } else {
    # Remove any rows with NA values
    bscores <- bscores[apply(bscores, 1, function(x) sum(is.na(x))) == 0, ]
  }
  
  N <- dim(bscores)[1]
  p <- dim(bscores)[2]
  
  tmp <- array(dim = c(p, p, N))
  
  for(ii in 1:N){
    tmp[ , , ii] <- bscores[ii, ] %*% t(bscores[ii, ])
  }
  
  out <- apply(tmp, 1:2, mean, na.rm = T)
    
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
V_matrix_effects <- function(Bscores, ipw_obj, alpha1, alpha2, 
                             trt.lvl1, trt.lvl2, effect, marginal,
                             na.rm = FALSE){

  N <- dim(Bscores)[1]
  p <- dim(Bscores)[2]
  a1 <- alpha1
  a2 <- alpha2
  t1 <- trt.lvl1
  t2 <- trt.lvl2
  
  ## replace any Bscores with 0 ##
  if(na.rm == FALSE) {
    Bscores[is.na(Bscores)] <- 0
  }
  
  fff <- ifelse(marginal == TRUE, 'marginal_outcomes', 'outcomes')
  
  if(effect == 'contrast'){   
    x <- ipw_obj[[fff]]$contrasts$group_resid
    if(marginal == TRUE){
      x <- x[a1, a2, ]
    } else {
      x <- x[a1, t1, a2, t2, , a1, t1, a2, t2]
    }
  } 
  else if(effect == 'outcome'){
    x <- ipw_obj[[fff]]$group_resid
    if(marginal == TRUE){
      x <- x[ , a1]
    } else {
      x <- x[  ,a1, t1]
    }
  }
  
  tmp <- array(dim = c(p+1, p+1, N))
  
  for(ii in 1:N){
    hold <- c(Bscores[ii, ], x[ii])
    tmp[ , , ii] <- hold %*% t(hold)
  }
  
  V <- apply(tmp, 1:2, mean, na.rm = T)

  return(V)
}
