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

  N <- dim(scores)[1]
  p <- dim(scores)[2]
  a1 <- alpha1
  a2 <- alpha2
  t1 <- trt.lvl1
  t2 <- trt.lvl2
  
  fff <- ifelse(marginal == TRUE, 'marginal_outcomes', 'outcomes')
  
  if(effect == 'contrast'){   
    hold_oal <- point_estimates[[fff]]$overall 
    hold_grp <- point_estimates[[fff]]$groups

    if(marginal == TRUE){
      pe <- hold_oal[a1] - hold_oal[a2]
      grp.pe <- hold_grp[ , a1] - hold_grp[, a2]
      x <- grp.pe - pe
    } else {
      pe <- hold_oal[a1, t1] - hold_oal[a2, t2]
      grp.pe <- hold_grp[ ,a1, t1] - hold_grp[, a2, t2]
      x <- grp.pe - pe
    }
  } 
  else if(effect == 'outcome'){
    x <- point_estimates[[fff]]$group_resid
    if(marginal == TRUE){
      x <- x[ , a1]
    } else {
      x <- x[  ,a1, t1]
    }
  }
  
  tmp <- array(dim = c(p+1, p+1, N))
  
  for(ii in 1:N){
    hold <- c(scores[ii, ], x[ii])
    tmp[ , , ii] <- hold %*% t(hold)
  }
  
  V <- apply(tmp, 1:2, mean, na.rm = T)

  return(V)
}
