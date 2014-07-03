#----------------------------------------#
######       bscore()v              ######
# compute psi_xp(B_i, X_i, theta)        #
#----------------------------------------#

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


#' Partial U
#'  
#' @param predictors
#' @return length(alpha) vector of IPW estimates
#' @export
Upartial <- function(y, G, A, B, data, weights, weight_dervs,
                     rescale.factor, na.rm = FALSE){
  out <- list()
  
  groups <- as.numeric(dimnames(weight_dervs)[[1]])
  predictors <- dimnames(weight_dervs[[2]])
  alphas <- as.numeric(dimnames(weight_dervs)[[3]])
  trt_lvls <- sort(unique(data[, A]))
  
  N <- length(groups)
  k <- length(alphas)
  l <- length(trt_lvls)
  p <- length(theta) - 1
  
  # replace any missing weights with 0
  if(na.rm == FALSE) {
    weight_dervs[is.na(weight_dervs)] <- 0
  }
  
  ## CALCULATE MARGINAL ESTIMATES ####
  ybar <- group_means(Y = y, A = A, G = G, a = NA, data = data)
  grp_est <- apply(weight_dervs, 2:3, function(x) x * ybar) * rescale.factor
  
  out$ipw$marginal_outcomes$groups <- grp_est
  out$ipw$marginal_outcomes$overall <- ymarg <- apply(grp_est, 2, mean, na.rm = TRUE)
  
  #   ## CALCULATE OVERALL EFFECTS ####
  #   out$ipw$effects$oe$overall <- outer(ymarg, ymarg, '-')
  #   hold_oe_grp <- array(dim = c(k, k, N), dimnames = list(alphas, alphas, groups))
  #   
  #   for(ii in 1:N){
  #     hold_oe_grp[ , , ii] <- outer(grp_est[ii, ], grp_est[ii, ], '-')
  #   }
  #   
  #   out$ipw$effects$oe$groups <- hold_oe_grp
  
  ## CALCULATE OUTCOME ESTIMATES PER TREATMENT LEVEL####
  
  hold_grp <- hold_grp_diff <- array(dim = c(N, p + 1, k, l),
                                     dimnames = list(groups, predictors, alphas, trt_lvls))
  hold_oal <- array(dim = c(p+1, k, l),
                    dimnames = list(predictors, alphas, trt_lvls))
  
  for(ll in 1:l){    
    a <- trt_lvls[ll]
    
    # Compute means per group
    ybar_trt <- group_means(Y = y, A = A, G = G, a = a, data = data)
    
    # Modify weights per treatment level
    weights_d_trt <- array(dim= c(N, p + 1, k),
                           dimnames = list(1:N, predictors, alphas))
    for(pp in 1:(p + 1)){
      weights_d_trt[ , pp, ] <- t(t(weight_dervs[ , pp, ])/((alphas^a*(1-alphas)^(1-a))))
    }
    
    # Compute estimates
    grp_est <- apply(weights_d_trt, 2:3, function(x) x * ybar_trt) * rescale.factor
    oal_est <- apply(grp_est, 2:3, mean, na.rm = TRUE)
    
    hold_grp[ , , , ll] <- grp_est
    hold_oal[ , , ll]   <- oal_est
    
    for(alpha in 1:k){
      for(g in 1:N){
        hold_grp_diff[g, , alpha, ll] <- hold_grp[g,  , alpha, ll] - hold_oal[ , alpha, ll]
      }
    }
  }
  
  out$ipw$outcomes <- list(groups = hold_grp, 
                           overall = hold_oal, 
                           group_resid = hold_grp_diff)
  
  ## CALCULATE EFFECT CONTRASTS ####
  
  hold_contrasts_oal <- array(dim = c(k, l, k, l, p + 1),
                              dimnames = list(alphas, trt_lvls, alphas, 
                                              trt_lvls, predictors))
  
  hold_contrast_grp <- array(dim = c(k, l, k, l, N, p + 1),
                             dimnames = list(alphas, trt_lvls, alphas, trt_lvls, 
                                             groups, predictors))
  
  hold_resid <- array(dim = c(k, l, k, l, N, k, l, k, l, p + 1),
                      dimnames = list(alphas, trt_lvls, alphas, trt_lvls, 
                                      groups, alphas, trt_lvls, alphas, trt_lvls, 
                                      predictors))
  for(pp in 1:(p+1)){
    hold_contrasts_oal[ , , , , pp] <- outer(hold_oal[pp, , ], hold_oal[pp, , ], '-')
    
    for(ii in 1:N){
      hold_contrast_grp[ , , , , ii, pp] <- outer(hold_grp[ii, pp, , ], hold_grp[ii, pp, , ], '-')
    }
    
    hold_resid[ , , , , , , , , , pp] <- outer(hold_contrast_grp[ , , , , , pp],
                                               hold_contrasts_oal[ , , , , pp], '-')
  }
  
  
  out$ipw$contrasts <- list(overall = hold_contrasts_oal,
                            groups  = hold_contrast_grp,
                            group_resid = hold_resid)
  
  ## DONE ####
  return(out)
}




#' V Matrix
#'  
#' @param predictors
#' @return V matrix
#' @export
V_matrix_effect <- function(Bscores, ipw_obj, alpha1, alpha2, trt.lvl1, trt.lvl2){

  N <- ipw_obj$summary$ngroups
  p <- ipw_obj$summary$npredictors
  
  tmp <- array(dim = c(p+2, p+2, N))

  for(ii in 1:N){
    hold <- c(Bscores[ii, ], ipw_obj$ipw$contrasts$group_resid[alpha1, trt.lvl1, alpha2, trt.lvl2, ii, alpha1, trt.lvl1, alpha2, trt.lvl2])
    tmp[ , , ii] <- hold %*% t(hold)
  }
  
  V <- apply(tmp, 1:2, sum, na.rm = T)/N

  return(V)
}


