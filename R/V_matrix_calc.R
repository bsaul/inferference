#-----------------------------------------------------------------------------#
#' V Matrix
#'  
#' @param scores
#' @return V matrix
#' @export
#' 
#-----------------------------------------------------------------------------#
V_matrix <- function(scores, 
                     point_estimates, 
                     alpha1, 
                     trt.lvl1, 
                     alpha2   = NA, 
                     trt.lvl2 = NA, 
                     effect, 
                     marginal){
  ## Necessary bits ##
  N <- dim(scores)[1]
  p <- dim(scores)[2]
  a1 <- alpha1
  a2 <- alpha2
  t1 <- trt.lvl1
  t2 <- trt.lvl2
  
  ## Grab the last element of the psi(O, theta) vector: psi_a, alpha ##
  fff <- ifelse(marginal == TRUE, 'marginal_outcomes', 'outcomes')
    
  if(effect == 'contrast'){   
    hold_oal <- point_estimates[[fff]]$overall 
    hold_grp <- point_estimates[[fff]]$groups

    if(marginal == TRUE){
      pe <- hold_oal[a1] - hold_oal[a2]
      grp.pe <- hold_grp[ , a1] - hold_grp[, a2]
      xx <- grp.pe - pe
    } else {
      pe <- hold_oal[a1, t1] - hold_oal[a2, t2]
      grp.pe <- hold_grp[ ,a1, t1] - hold_grp[, a2, t2]
      xx <- grp.pe - pe
    }
  } 
  else if(effect == 'outcome'){
    if(marginal == TRUE){
      xx <- hold_grp[ , a1] - hold_oal[a1]
    } else {
      xx <- hold_grp[  , a1, t1] - hold_oal[a1, t1]
    }
  }
  
  ## Create the V matrix for each group ##
  V_i <- array(dim = c(p+1, p+1, N))
  
  for(ii in 1:N){
    hold <- c(scores[ii, ], xx[ii])
    V_i[ , , ii] <- hold %*% t(hold)
  }
  
  ## Create final V matrix by taking the mean of each element across groups ##
  V <- apply(V_i, 1:2, sum, na.rm = T)/N

  return(V)
}
