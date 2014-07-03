#' Calculate IPW direct effect estimates
#'  
#' @param TBD
#' @return TBD
#' @export


ipw_de <- function(obj, alpha, trt1, trt2){
  N <- obj[[1]]$summary$ngroups
  
  point_estimate <- obj[[1]]$ipw$contrasts$overall[alpha, trt1, alpha, trt2]
  V <- V_matrix_effect(obj[[2]], obj[[1]], alpha, alpha, trt1, trt2)
  U21 <- apply(obj[[3]]$ipw$contrasts$group_resid[alpha, trt1, alpha, trt2, , alpha, trt1, alpha, trt2, ], 
                2, sum, na.rm = T)/N
  
  V21 <- V[dim(V)[1], 1:(dim(V)[2] - 1)]
  V11 <- V[1:(dim(V)[1] - 1), 1:(dim(V)[2] - 1)]
  V22 <- V[dim(V)[1], dim(V)[2]]
  ## VARIANCE Estimate
  variance <- ((U21 - 2*V21) %*% solve(V11) %*% U21 + V22)/N
  
  print(point_estimate + 1.96 * sqrt(variance))
}


