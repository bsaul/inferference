#-----------------------------------------------------------------------------#
#' V Matrix
#'  
#' @param predictors
#' @return V matrix
#' @export
#' 
#-----------------------------------------------------------------------------#
V_matrix <- function(Bscores, ipw_obj, 
                     alpha1, alpha2, 
                     trt.lvl1 = NULL, trt.lvl2 = NULL, 
                     effect, marginal,
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
    hold_oal <- ipw_obj[[fff]]$overall 
    hold_grp <- ipw_obj[[fff]]$groups

    #x <- grp.pe - pe
    #x <- ipw_obj[[fff]]$contrasts$group_resid
    if(marginal == TRUE){
      pe <- hold_oal[a1] - hold_oal[a2]
      grp.pe <- hold_grp[ ,a1] - hold_grp[, a2]
      x <- grp.pe - pe
      #x <- x[a1, a2, ]
    } else {
      pe <- hold_oal[a1, t1] - hold_oal[a2, t2]
      grp.pe <- hold_grp[ ,a1, t1] - hold_grp[, a2, t2]
      x <- grp.pe - pe
      #x <- x[a1, t1, a2, t2, , a1, t1, a2, t2]
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
