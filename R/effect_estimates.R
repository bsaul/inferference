#-----------------------------------------------------------------------------#
#' Calculate Effect estimates
#'  
#' @param TBD
#' @return TBD
#' @export
#-----------------------------------------------------------------------------#

calc_effect <- function(obj, alpha1, alpha2, trt.lvl1, trt.lvl2){
  N <- obj[[1]]$summary$ngroups
  a1 <- alpha1
  a2 <- alpha2
  t1 <- trt.lvl1
  t2 <- trt.lvl2
  
  pe <- obj$point_estimates$outcomes$contrasts$overall[a1, t1, a2, t2]
  
  ## VARIANCE ESTIMATION ####
  bscores <- obj$bscores
  
  V <- V_matrix(bscores, obj[[1]], a1, a2, t1, t2, 
                effect = 'contrast', marginal = FALSE )
  
  U21 <- apply(obj[[3]]$outcomes$contrasts$group_resid[a1, t1, a2, t2, , a1, t1, a2, t2, ], 
               2, mean, na.rm = T)
  
  V21 <- V[dim(V)[1], 1:(dim(V)[2] - 1)]
  V11 <- V[1:(dim(V)[1] - 1), 1:(dim(V)[2] - 1)]
  V22 <- V[dim(V)[1], dim(V)[2]]
  ## VARIANCE Estimate
  ve <- ((U21 - 2*V21) %*% solve(V11) %*% U21 + V22)/N
  me <- 1.96 * sqrt(ve)
  
  toprint <- paste0('Estimate: ', round(pe, 2), ' 95% CI: (', 
                round(pe - me, 2), ', ', round(pe + me, 2), ')' )
  
  print(toprint)
  
  out <- c(pe, ve, pe-me, pe+me)
  return(out)
}

#-----------------------------------------------------------------------------#
#' Calculate Direct effect estimates
#'  
#' @param TBD
#' @return TBD
#' @export
#-----------------------------------------------------------------------------#

direct_effect <- function(obj, alpha, trt.lvl1 = '0', trt.lvl2 = '1'){
  out <- calc_effect(obj, alpha, alpha, trt.lvl1, trt.lvl2)
  return(out)
}

#-----------------------------------------------------------------------------#
#' Calculate indirect effect estimates
#'  
#' @param TBD
#' @return TBD
#' @export
#-----------------------------------------------------------------------------#

indirect_effect <- function(obj, alpha1, alpha2, trt.lvl = '0'){
  out <- calc_effect(obj, alpha1, alpha2, trt.lvl, trt.lvl)
  return(out)
}


#-----------------------------------------------------------------------------#
#' Calculate indirect effect estimates
#'  
#' @param TBD
#' @return TBD
#' @export
#-----------------------------------------------------------------------------#

total_effect <- function(obj, alpha1, alpha2, trt.lvl1 = '0', trt.lvl2 = '1'){
  out <- calc_effect(obj, alpha1, alpha2, trt.lvl1, trt.lvl2)
  return(out)
}