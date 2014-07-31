#-----------------------------------------------------------------------------#
#' Calculate IPW estimates
#'
#'  Computes either outcome or effect estimates from the object output by 
#'  \code{\link{run_interference}}.  
#'  
#'  @details See \code{\link{direct_effect}}, \code{\link{indirect_effect}},
#'  \code{\link{total_effect}}, and \code{\link{overall_effect}} for convenient
#'  wrappers of \code{calc_effect} to compute common causal effects.
#'  
#'  This table summarizes the value that \code{calc_effect} returns.
#'  \tabular{llc}{
#'  Marginal \tab Effect    \tab Value returned \cr
#'  FALSE    \tab 'outcome' \tab \eqn{\hat{Y}(trt.lvl1, alpha1)}{Yhat(trt.lvl1, alpha1)} \cr
#'  TRUE     \tab 'outcome' \tab \eqn{\hat{Y}(alpha1)}{Yhat(alpha1)}  \cr
#'  FALSE    \tab 'contrast' \tab 
#'  \eqn{\hat{Y}(trt.lvl1, alpha1) - \hat{Y}(trt.lvl2, alpha2)}{Yhat(trt.lvl1, alpha1) - Yhat(trt.lvl2, alpha2)} \cr
#'  TRUE    \tab 'contrast' \tab 
#'  \eqn{\hat{Y}(alpha1) - \hat{Y}(alpha2)}{Yhat(alpha1) - Yhat(alpha2)} \cr
#' }
#'  
#' @param obj the name of the object created by \code{\link{run_interference}}
#' @param alpha1 the allocation scheme for the outcome of interest or the first
#' scheme in the constrast of interest. See details.
#' @param trt.lvl1 the treatment level for the outcome of interest or the first
#' treatment in the constrast of interest. If marginal = TRUE, this is ignored.
#' @param alpha2 the second allocation scheme for the contrast of interest.
#' Ignored if effect = 'outcome'.
#' @param trt.lvl2  the second treatment in the constrast of interest. 
#' If marginal = TRUE or effect = 'outcome', this is ignored.
#' @param effect either 'contrast' or 'outcome'
#' @param marginal TRUE or FALSE
#' @param conf.level Confidence level for confidence intervals. Defaults to 0.95.
#' @param print TRUE/FALSE. If TRUE, the point estimates and confidence interval
#' are printed to the console. 
#' @return A \code{data.frame} with 1 record and 4 variables: point (the point
#'  estimate), variance (the variance estimate), ll (the lower bound of the 
#'  confidence interval), and ul (the upper bound of the confidence interval).
#' @export
#-----------------------------------------------------------------------------#

calc_effect <- function(obj, 
                        alpha1, 
                        trt.lvl1, 
                        alpha2 = NA, 
                        trt.lvl2 = NA,
                        effect,
                        marginal,
                        rescale.factor = 1,
                        conf.level = 0.95,
                        print = FALSE){
  
  # Print error if either estimates with alpha1 have been computed 
  # or a constrast is being estimated when estimates for alpha2
  # have not been computed
  if (!(alpha1 %in% obj[[1]]$summary$alphas) | 
     (effect == 'contrast' & !(alpha2 %in% obj[[1]]$summary$alphas))){
    stop(paste('At least one of the chosen coverage levels has not been estimated.\n',
               'Select from the following: \n', 
               paste(obj[[1]]$summary$alphas, collapse = ' ')))
  }
  
  N <- obj[[1]]$summary$ngroups
  a1 <- as.character(alpha1)
  a2 <- as.character(alpha2)
  t1 <- as.character(trt.lvl1)
  t2 <- as.character(trt.lvl2)
  
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
  V <- V_matrix(Bscores = obj$bscores, 
                ipw_obj = obj$point_estimates, 
                alpha1 = a1, alpha2 = a2, 
                trt.lvl1 = t1, trt.lvl2 = t2, 
                effect = effect, marginal = marginal)
  
  U21 <- apply(U_grp_diff, 2, mean, na.rm = T)
  
  V21 <- V[dim(V)[1], 1:(dim(V)[2] - 1)]
  V11 <- V[1:(dim(V)[1] - 1), 1:(dim(V)[2] - 1)]
  V22 <- V[dim(V)[1], dim(V)[2]]
 
  ## VARIANCE Estimate
  ve <- ((U21 - 2*V21) %*% solve(V11) %*% U21 + V22)/N
  
  ## CONFIDENCE INTERVALS
  qq <- conf.level + (1 - conf.level)/2
  
  me <- qnorm(qq) * sqrt(ve)
  
  if(print == TRUE){
    toprint <- paste0('Estimate: ', round(pe, 2), conf.level*100, '% CI: (', 
                round(pe - me, 2), ', ', round(pe + me, 2), ')' )
    print(toprint)
  }
  
  out <- data.frame(point = pe,
                    variance = ve, 
                    ll = pe - me, 
                    ul = pe + me)
  return(out)
}

#-----------------------------------------------------------------------------#
#' Calculate Direct Effect estimates
#'  
#' @description By default, this function computes:
#' \eqn{\hat{Y}(0, alpha) - \hat{Y}(1, alpha)}{Yhat(0, alpha) - Yhat(1, alpha)}
#'  
#' @param obj the name of the object created by \code{\link{run_interference}}
#' @param alpha the allocation scheme for which to estimate direct effects
#' @param trt.lvl1 Defaults to 0.
#' @param trt.lvl2 Defaults to 1.
#' @param print see \code{\link{calc_effect}}
#' @param conf.level see \code{\link{calc_effect}}
#' @return See \code{\link{calc_effect}}.
#' @export
#-----------------------------------------------------------------------------#

direct_effect <- function(obj, 
                          alpha, 
                          trt.lvl1 = 0, 
                          trt.lvl2 = 1,
                          print = FALSE, 
                          conf.level = 0.95){
  out <- calc_effect(obj, alpha, trt.lvl1, alpha, trt.lvl2,
                     effect = 'contrast', marginal = FALSE,
                     print = print, conf.level = conf.level)
  return(out)
}

#-----------------------------------------------------------------------------#
#' Calculate indirect effect estimates
#'  
#' @description By default, this function computes:
#' \eqn{\hat{Y}(0, alpha1) - \hat{Y}(0, alpha2)}{Yhat(0, alpha1) - Yhat(0, alpha2)}
#'  
#'  
#' @param obj the name of the object created by \code{\link{run_interference}}
#' @param alpha1 the allocation scheme for which to estimate indirect effects
#' @param alpha2 the allocation scheme for which to estimate indirect effects
#' @param trt.lvl Defaults to 0.
#' @param print see \code{\link{calc_effect}}
#' @param conf.level see \code{\link{calc_effect}}
#' @return See \code{\link{calc_effect}}.
#' @export
#-----------------------------------------------------------------------------#

indirect_effect <- function(obj, 
                            alpha1, 
                            alpha2, 
                            trt.lvl = 0, 
                            print = FALSE, 
                            conf.level = 0.95){
  
  out <- calc_effect(obj, alpha1, trt.lvl, alpha2, trt.lvl,
                     effect = 'contrast', marginal = FALSE,
                     print = print, conf.level = conf.level)
  return(out)
}

#-----------------------------------------------------------------------------#
#' Calculate Total effect estimates
#'
#' @description By default, this function computes:
#' \eqn{\hat{Y}(0, alpha1) - \hat{Y}(1, alpha2)}{Yhat(0, alpha1) - Yhat(1, alpha2)}
#'  
#' @param obj the name of the object created by \code{\link{run_interference}}
#' @param alpha1 the allocation scheme for which to estimate total effects
#' @param alpha2 the allocation scheme for which to estimate total effects
#' @param trt.lvl1 Defaults to 0.
#' @param trt.lvl2 Defaults to 1.
#' @param print see \code{\link{calc_effect}}
#' @param conf.level see \code{\link{calc_effect}}
#' @return See \code{\link{calc_effect}}.
#' @export
#-----------------------------------------------------------------------------#

total_effect <- function(obj, 
                         alpha1, 
                         alpha2, 
                         trt.lvl1 = 0, 
                         trt.lvl2 = 1, 
                         print = FALSE, 
                         conf.level = 0.95){
  
  out <- calc_effect(obj, alpha1, trt.lvl1, alpha2, trt.lvl2,
                     effect = 'contrast', marginal = FALSE,
                     print = print, conf.level = conf.level)
  return(out)
}

#-----------------------------------------------------------------------------#
#' Calculate Overall effect estimates
#' 
#' @description By default, this function computes:
#' \eqn{\hat{Y}(alpha1) - \hat{Y}(lpha2)}{Yhat(alpha1) - Yhat(alpha2)}
#' 
#' @param obj the name of the object created by \code{\link{run_interference}}
#' @param alpha1 the allocation scheme for which to estimate overall effects
#' @param alpha2 the allocation scheme for which to estimate overall effects
#' @param print see \code{\link{calc_effect}}
#' @param conf.level see \code{\link{calc_effect}}
#' @return See \code{\link{calc_effect}}.
#' @export
#-----------------------------------------------------------------------------#

overall_effect <- function(obj, 
                           alpha1, 
                           alpha2, 
                           print = FALSE, 
                           conf.level = 0.95){
  
  out <- calc_effect(obj, alpha1, trt.lvl1 = NA, alpha2, trt.lvl2 = NA,
                     effect = 'contrast', marginal = TRUE,
                     print = print, conf.level = conf.level)
  return(out)
}