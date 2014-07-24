#-----------------------------------------------------------------------------#
#' Calculate Effect estimates
#'  
#' @param TBD
#' @return TBD
#' @export
#-----------------------------------------------------------------------------#

calc_effect <- function(obj, 
                        alpha1, trt.lvl1, 
                        alpha2, trt.lvl2,
                        effect,
                        marginal,
                        print = FALSE)
{
  N <- obj[[1]]$summary$ngroups
  a1 <- alpha1
  a2 <- alpha2
  t1 <- trt.lvl1
  t2 <- trt.lvl2
  
  fff <- ifelse(marginal == TRUE, 'marginal_outcomes', 'outcomes')
  
  hold_oal <- obj$point_estimates[[fff]]$overall
  U_hold_oal <- obj$Upart[[fff]]$overall 
  U_hold_grp <- obj$Upart[[fff]]$groups
  
  if(effect == 'contrast'){
    if(marginal == TRUE){
      pe <- hold_oal[a1] - hold_oal[a2]
      U_pe <- U_hold_oal[ ,a1] - U_hold_oal[ ,a2]
      U_pe_grp <- U_hold_grp[ , ,a1] - U_hold_grp[ , ,a2]  
      U_grp_diff <- t(U_pe - t(U_pe_grp))
    } else {
      pe <- hold_oal[a1, t1] - hold_oal[a2, t2]
      U_pe <- U_hold_oal[ ,a1, t1] - U_hold_oal[ ,a2, t2]
      U_pe_grp <- U_hold_grp[ , ,a1, t1] - U_hold_grp[ , ,a2, t2]  
      U_grp_diff <- t(U_pe - t(U_pe_grp))
    }
  } 
  else if(effect == 'outcome'){
    if(marginal == TRUE){
      pe <- hold_oal[a1]
      U_pe <- U_hold_oal[ ,a1]
      U_pe_grp <- U_hold_grp[ , ,a1]
      U_grp_diff <- t(U_pe - t(U_pe_grp))
    } else {
      pe <- hold_oal[a1, t1]
      U_pe <- U_hold_oal[ ,a1, t1]
      U_pe_grp <- U_hold_grp[ , ,a1, t1]
      U_grp_diff <- t(U_pe - t(U_pe_grp))
    }
  }
  
  ## VARIANCE ESTIMATION ####
  
  # V matrix
  V <- V_matrix(Bscores = obj$bscores, ipw_obj = obj[[1]], 
                alpha1 = a1, alpha2 = a2, 
                trt.lvl1 = t1, trt.lvl2 = t2, 
                effect = effect, marginal = marginal)

  
  U21 <- apply(U_grp_diff, 2, mean, na.rm = T)
  
  V21 <- V[dim(V)[1], 1:(dim(V)[2] - 1)]
  V11 <- V[1:(dim(V)[1] - 1), 1:(dim(V)[2] - 1)]
  V22 <- V[dim(V)[1], dim(V)[2]]
 
  ## VARIANCE Estimate
  ve <- ((U21 - 2*V21) %*% solve(V11) %*% U21 + V22)/N
  me <- 1.96 * sqrt(ve)
  
  if(print == TRUE){
    toprint <- paste0('Estimate: ', round(pe, 2), ' 95% CI: (', 
                round(pe - me, 2), ', ', round(pe + me, 2), ')' )
    print(toprint)
  }
  
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

direct_effect <- function(obj, alpha, 
                          trt.lvl1 = '0', trt.lvl2 = '1',
                          print = FALSE)
{
  out <- calc_effect(obj, alpha, trt.lvl1, alpha, trt.lvl2,
                     effect = 'contrast', marginal = FALSE,
                     print = print)
  return(out)
}

#-----------------------------------------------------------------------------#
#' Calculate indirect effect estimates
#'  
#' @param TBD
#' @return TBD
#' @export
#-----------------------------------------------------------------------------#

indirect_effect <- function(obj, alpha1, alpha2, 
                            trt.lvl = '0', 
                            print = FALSE)
{
  out <- calc_effect(obj, alpha1, trt.lvl, alpha2, trt.lvl,
                     effect = 'contrast', marginal = FALSE,
                     print = print)
  return(out)
}


#-----------------------------------------------------------------------------#
#' Calculate Total effect estimates
#'  
#' @param TBD
#' @return TBD
#' @export
#-----------------------------------------------------------------------------#

total_effect <- function(obj, alpha1, alpha2, 
                         trt.lvl1 = '0', trt.lvl2 = '1', 
                         print = FALSE)
{
  out <- calc_effect(obj, alpha1, trt.lvl1, alpha2, trt.lvl2,
                     effect = 'contrast', marginal = FALSE,
                     print = print)
  return(out)
}

#-----------------------------------------------------------------------------#
#' Calculate Overall effect estimates
#'  
#' @param TBD
#' @return TBD
#' @export
#-----------------------------------------------------------------------------#

overall_effect <- function(obj, alpha1, alpha2, 
                           print = FALSE)
{
  out <- calc_effect(obj, alpha1, trt.lvl1 = NA, alpha2, trt.lvl2 = NA,
                     effect = 'contrast', marginal = TRUE,
                     print = print)
  return(out)
}