#-----------------------------------------------------------------------------#
#' Calculate IPW point estimates
#'  
#' @param y vector of unweighted group means (i.e. output from \code{\link{group_means}}
#' @param a treatment level (0,1) for which to compute weights. Defaults to NULL which returns marginal.
#' @param weights weight matrix/array to use
#' @param rescale.factor factor by which to rescale values. Defaults to 1.
#' @param na.rm exclude groups missing weights? defaults to FALSE.
#' @return length(alpha) vector of IPW estimates
#' @export
#-----------------------------------------------------------------------------#

ipw_point_estimates <- function(y, 
                                G, 
                                A, 
                                data, 
                                weights,
                                rescale.factor = 1, 
                                set.NA.to.0 = TRUE){
  
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
  if(set.NA.to.0 == TRUE) {
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

  ## DONE ####
  return(out)
}