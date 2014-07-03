#-----------------------------------------------------------------------------#
#' Calculate IPW point estimates
#'  
#' @param y vector of unweighted group means (i.e. output from \code{\link{group_means}}
#' @param a treatment level (0,1) for which to compute weights. Defaults to NULL which returns marginal.
#' @param weights weight matrix to use
#' @param rescale.factor factor by which to rescale values
#' @param na.rm exclude groups missing weights? defaults to FALSE.
#' @return length(alpha) vector of IPW estimates
#' @export
#-----------------------------------------------------------------------------#

ipw_point_estimates <- function(y, G, A, B, data, weights,
                                rescale.factor, na.rm = FALSE){
  out <- list()
  
  groups <- as.numeric(dimnames(weights)[[1]])
  alphas <- as.numeric(dimnames(weights)[[2]])
  trt_lvls <- sort(unique(data[, A]))
  
  N <- length(groups)
  k <- length(alphas)
  l <- length(trt_lvls)
  
  # Summary information
  out$summary <- list(ngroups = N, 
                      nalphas = k,
                      ntreatments = l,
                      alphas = alphas,
                      treatments = trt_lvls)  

  # replace any missing weights with 0
  if(na.rm == FALSE) {
    weights[is.na(weights)] <- 0
  }
    
  ## CALCULATE MARGINAL ESTIMATES ####
  ybar <- group_means(Y = y, A = A, G = G, a = NA, data = data)
  grp_est <- apply(weights, 2, function(x) x * ybar) * rescale.factor

  out$marginal_outcomes$groups <- grp_est
  out$marginal_outcomes$overall <- ymarg <- apply(grp_est, 2, mean, na.rm = TRUE)
  
  ## CALCULATE OVERALL EFFECTS ####
  out$ipw$marginal_outcomes$contrasts$overall <- outer(ymarg, ymarg, '-')
  hold_oe_grp <- array(dim = c(k, k, N), dimnames = list(alphas, alphas, groups))
  
  for(ii in 1:N){
    hold_oe_grp[ , , ii] <- outer(grp_est[ii, ], grp_est[ii, ], '-')
  }
  
  out$marginal_outcomes$contrasts$overall <- hold_oe_grp

  ## CALCULATE OUTCOME ESTIMATES PER TREATMENT LEVEL####

  hold_grp <- hold_grp_diff <- array(dim = c(N, k, l),
                    dimnames = list(groups, alphas, trt_lvls))
  hold_oal <- array(dim = c(k, l),
                    dimnames = list(alphas, trt_lvls))
  
  for(ll in 1:l){    
    a <- trt_lvls[ll]

    # Compute means per group
    ybar_trt <- group_means(Y = y, A = A, G = G, a = a, data = data)
    
    # Modify weights per treatment level
    weights_trt <- t(t(weights)/((alphas^a*(1-alphas)^(1-a))))
    
    # Compute estimates
    grp_est <- apply(weights_trt, 2, function(x) x * ybar_trt) * rescale.factor
    oal_est <- apply(grp_est, 2, mean, na.rm = TRUE)
    
    hold_grp[ , , ll] <- grp_est
    hold_oal[ , ll]   <- oal_est
    
    for(kk in 1:k){
      for(ii in 1:N){
        hold_grp_diff[ii, kk, ll] <- hold_grp[ii, kk, ll] - hold_oal[kk, ll]
      }
    }
  }

  out$outcomes <- list(groups = hold_grp, 
                       overall = hold_oal, 
                       group_resid = hold_grp_diff)
  
  ## CALCULATE EFFECT CONTRASTS ####
  hold_contrasts_oal <- outer(hold_oal, hold_oal, '-')
  
  hold_contrast_grp <- array(dim = c(k, l, k, l, N),
                             dimnames = list(alphas, trt_lvls,
                                             alphas, trt_lvls, groups))

  for(ii in 1:N){
    hold_contrast_grp[ , , , , ii] <- outer(hold_grp[ii, , ], hold_grp[ii, , ], '-')
  }
  
  out$outcomes$contrasts <- list(overall = hold_contrasts_oal,
                                 groups  = hold_contrast_grp,
                                 group_resid = outer(hold_contrast_grp, 
                                                hold_contrasts_oal, '-'))
  ## DONE ####
  return(out)
}