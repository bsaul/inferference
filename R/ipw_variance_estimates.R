#' B Score
#'  
#' @param predictors
#' @return length(alpha) vector of IPW estimates
#' @export

bscore <- function(predictors, B, G, theta, data){
  N <- length(unique(data[, G]))
  
  Bscores <- matrix(nrow = N, ncol = length(theta))  
  for(ii in 1:N){
    grp <- analysis_c[analysis_c$group == ii, ]
    X <- as.matrix(cbind(1, grp[, predictors]))
    Bscores[ii, ] <- cmp_scores_ll(theta = theta, B = grp[ , B], X = X)
  }
  
  return(Bscores)
}

#' V Matrix
#'  
#' @param predictors
#' @return V matrix
#' @export

V_matrix <- function(Bscores, ipw_obj){
  k <- ipw_obj$summary$nalphas
  N <- ipw_obj$summary$ngroups
  l <- ipw_obj$summary$ntreatments
  p <- ipw_obj$summary$npredictors
  alphas <- ipw_obj$summary$alphas
  
  hold_Vi <- array(dim = c(p+2, p+2, N, k, l))
  hold_V  <- array(dim = c(p+2, p+2, k, l),
                   dimnames = list(1:(p+2), 1:(p+2), 
                                   alphas, ipw_obj$summary$treatments))
  ## OUTCOMES ##
  for(ll in 1:l){
    for(kk in 1:k){
      psi_i_a <- cbind(Bscores, ipw_obj$ipw$outcome$group_resid[, kk, ll])
      
      for(ii in 1:N){
        hold_Vi[ , , ii, kk, ll] <- psi_i_a[ii, ]  %*% t(psi_i_a[ii, ])  
      }
      hold_V[ , , kk, ll] <- apply(hold_Vi, 1:2, sum, na.rm = T)/N  
    }
  }
  
  ## CONTRASTS ##
  for(ll in 1:l){
    for(kk in 1:k){
      for(pp in 1:p){
        psi_i_eff <- cbind(Bscores, ipw_obj$ipw$contrasts$group_resid[, kk, ll])
        
        for(ii in 1:N){
          hold_Vi[ , , ii, kk, ll] <- psi_i_a[ii, ]  %*% t(psi_i_a[ii, ])  
        }
      }
      
      

      hold_V[ , , kk, ll] <- apply(hold_Vi, 1:2, sum, na.rm = T)/N  
    }
  }

  return(hold_V)
}


#' Partial U
#'  
#' @param predictors
#' @return length(alpha) vector of IPW estimates
#' @export
partialU <- function(y, G, A, B, data, weights, weight_dervs,
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

