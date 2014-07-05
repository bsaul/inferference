#' Calculate IPW direct effect estimates
#'  
#' @param TBD
#' @return TBD
#' @export


ipw_de <- function(obj, alpha, trt.lvl1, trt.lvl2){
  N <- obj[[1]]$summary$ngroups
  a <- alpha
  t1 <- trt.lvl1
  t2 <- trt.lvl2
    
  pe <- obj$point_estimates$outcomes$contrasts$overall[a, t1, a, t2]
  bscores <- obj$bscores

  V <- V_matrix(bscores, obj[[1]], a, a, t1, t2, 
                effect = 'contrast', marginal = FALSE )
  
  U21 <- apply(obj[[3]]$outcomes$contrasts$group_resid[a, t1, a, t2, , a, t1, a, t2, ], 
                2, mean, na.rm = T)
  
  V21 <- V[dim(V)[1], 1:(dim(V)[2] - 1)]
  V11 <- V[1:(dim(V)[1] - 1), 1:(dim(V)[2] - 1)]
  V22 <- V[dim(V)[1], dim(V)[2]]
  ## VARIANCE Estimate
  ve <- ((U21 - 2*V21) %*% solve(V11) %*% U21 + V22)/N
  me <- 1.96 * sqrt(ve)
  
  out <- paste0('Estimate: ', round(pe, 2), ' 95% CI: (', 
                round(pe - me, 2), ', ', round(pe + me, 2), ')' )
   
  return(print(out))
}


