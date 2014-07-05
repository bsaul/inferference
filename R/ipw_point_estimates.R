#-----------------------------------------------------------------------------#
#' Calculate IPW point estimates
#'  
#' @param y vector of unweighted group means (i.e. output from \code{\link{group_means}}
#' @param a treatment level (0,1) for which to compute weights. Defaults to NULL which returns marginal.
#' @param weights weight matrix/array to use
#' @param rescale.factor factor by which to rescale values
#' @param na.rm exclude groups missing weights? defaults to FALSE.
#' @return length(alpha) vector of IPW estimates
#' @export
#-----------------------------------------------------------------------------#

ipw_point_estimates <- function(y, G, A, B, data, weights,
                                 rescale.factor, na.rm = FALSE){
  
  ## DEFINE OBJECTS NEEDED FOR FUNCTION ##
  out <- list()
  
  groups <- as.numeric(dimnames(weights)[[1]])
  alphas <- as.numeric(dimnames(weights)[[length(dim(weights))]])
  trt_lvls <- sort(unique(data[, A]))
  
  N <- length(groups)
  k <- length(alphas)
  l <- length(trt_lvls)
  p <- ifelse(is.matrix(weights), 1, dim(weights)[2])
  
  ## Generalize function to work on arrays. Add a dimension to matrices. ##
  if(is.matrix(weights)){
    weights <- array(c(weights, 1), dim=c(N, 1, k))
    predictors <- NULL
  } else {
    predictors <- dimnames(weights)[[2]]
  }
  
  ## replace any missing weights with 0 ##
  if(na.rm == FALSE) {
    weights[is.na(weights)] <- 0
  }
  
  ## SUMMARY ## 
  out$summary <- list(ngroups = N, 
                      nalphas = k,
                      ntreatments = l,
                      alphas = alphas,
                      treatments = trt_lvls)  
  
  ## CALCULATE MARGINAL ESTIMATES ####
  ybar <- group_means(Y = y, A = A, G = G, a = NA, data = data)
  
  grp_est <- apply(weights, 2:3, function(x) x * ybar) * rescale.factor
  dimnames(grp_est) <- list(groups, predictors, alphas)
  
  oa_est <- apply(grp_est, 2:3, mean, na.rm = TRUE)

  out$marginal_outcomes$groups <- drop(grp_est)
  out$marginal_outcomes$overall <- drop(oa_est)
  
  hold_diff <- array(dim = c(N, p, k),
                     dimnames = list(groups, predictors, alphas))
  
  for(pp in 1:p){
    hold_diff[ , pp, ] <- grp_est[ , pp, ] - oa_est[pp, ]
  }
  
  out$marginal_outcomes$group_resid <- drop(hold_diff)
  
  ## CALCULATE OVERALL EFFECTS ####
  
  hold_oe_overall <- array(dim = c(k, k, p),
                           dimnames = list(alphas, alphas, predictors))
  hold_oe_grp <- hold_oe_diff <- array(dim = c(k, k, N, p), 
                       dimnames = list(alphas, alphas, groups, predictors))
  
  for(pp in 1:p){
    hold_oe_overall[ , , pp] <- outer(oa_est[pp, ], oa_est[pp, ], '-')
    
    for(ii in 1:N){
      hold_oe_grp[ , , ii, pp] <- outer(grp_est[ii, pp, ], grp_est[ii, pp, ], '-')
      hold_oe_diff[ , , ii, pp] <- hold_oe_grp[ , , ii, pp] - hold_oe_overall[ , , pp] 
    }
  }
  
  out$marginal_outcomes$contrasts$overall <- drop(hold_oe_overall)
  out$marginal_outcomes$contrasts$groups <- drop(hold_oe_grp)
  out$marginal_outcomes$contrasts$group_resid <- drop(hold_oe_diff)
  
  ## CALCULATE OUTCOME ESTIMATES PER TREATMENT LEVEL####
  
  hold_grp <- hold_grp_diff <- array(dim = c(N, p, k, l),
                                     dimnames = list(groups, predictors, 
                                                     alphas, trt_lvls))
  hold_oal <- array(dim = c(p, k, l),
                    dimnames = list(predictors, alphas, trt_lvls))
  
  for(ll in 1:l){    
    a <- trt_lvls[ll]
    
    # Compute means per group
    ybar_trt <- group_means(Y = y, A = A, G = G, a = a, data = data)
    
    # Modify weights per treatment level
    weights_trt <- array(dim= c(N, p, k),
                         dimnames = list(1:N, predictors, alphas))
    for(pp in 1:p){
      weights_trt[ , pp, ] <- t(t(weights[ , pp, ])/((alphas^a*(1-alphas)^(1-a))))
    }
    
    #weights_trt <- t(t(weights)/((alphas^a*(1-alphas)^(1-a))))
    
    # Compute estimates
    grp_est <- apply(weights_trt, 2:3, function(x) x * ybar_trt) * rescale.factor
    oal_est <- apply(grp_est, 2:3, mean, na.rm = TRUE)
    
    hold_grp[ , , , ll] <- grp_est
    hold_oal[ , , ll]   <- oal_est
    
    for(kk in 1:k){
      for(ii in 1:N){
        hold_grp_diff[ii, , kk, ll] <- hold_grp[ii,  , kk, ll] - hold_oal[ , kk, ll]
      }
    }
  }
  
  out$outcomes <- list(groups = drop(hold_grp), 
                       overall = drop(hold_oal), 
                       group_resid = drop(hold_grp_diff))
  
  ## CALCULATE EFFECT CONTRASTS ####
  
  hold_contrasts_oal <- array(dim = c(k, l, k, l, p),
                              dimnames = list(alphas, trt_lvls, alphas, 
                                              trt_lvls, predictors))
  
  hold_contrast_grp <- array(dim = c(k, l, k, l, N, p),
                             dimnames = list(alphas, trt_lvls, alphas, trt_lvls, 
                                             groups, predictors))
  
  hold_resid <- array(dim = c(k, l, k, l, N, k, l, k, l, p),
                      dimnames = list(alphas, trt_lvls, alphas, trt_lvls, 
                                      groups, alphas, trt_lvls, alphas, trt_lvls, 
                                      predictors))
  
  for(pp in 1:p){
    hold_contrasts_oal[ , , , , pp] <- outer(hold_oal[pp, , ], hold_oal[pp, , ], '-')
    
    for(ii in 1:N){
      hold_contrast_grp[ , , , , ii, pp] <- outer(hold_grp[ii, pp, , ], hold_grp[ii, pp, , ], '-')
    }
    
    hold_resid[ , , , , , , , , , pp] <- outer(hold_contrast_grp[ , , , , , pp],
                                               hold_contrasts_oal[ , , , , pp], '-')
  }
  
  out$outcomes$contrasts <- list(overall = drop(hold_contrasts_oal),
                                 groups  = drop(hold_contrast_grp),
                                 group_resid = drop(hold_resid))
  
  ## DONE ####
  return(out)
}