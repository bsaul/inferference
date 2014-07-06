#-----------------------------------------------------------------------------#
#' V Matrix
#'  
#' @param predictors
#' @return V matrix
#' @export
#' 
#-----------------------------------------------------------------------------#
V_matrix <- function(Bscores, ipw_obj, alpha1, alpha2, 
                             trt.lvl1, trt.lvl2, effect, marginal,
                             set.NA.to.0  = TRUE){

  N <- dim(Bscores)[1]
  p <- dim(Bscores)[2]
  a1 <- alpha1
  a2 <- alpha2
  t1 <- trt.lvl1
  t2 <- trt.lvl2
  
  ## replace any Bscores with 0 ##
  if(set.NA.to.0 == TRUE) {
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
