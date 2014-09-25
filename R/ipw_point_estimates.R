#-----------------------------------------------------------------------------#
#' Calculate IPW point estimates
#'  
#' @param outcome quoted name of outcome variable in \code{data}
#' @param groups quoted name of group variable in \code{data}
#' @param treatment quoted name of treatment variable in \code{data}
#' @param weights weight matrix/array to use from either \code{\link{wght_matrix}}
#' or \code{\link{wght_deriv_array}}
#' @return list containing point estimates for marginal outcomes and estimates
#' per treatment level
#' @export
#-----------------------------------------------------------------------------#

ipw_point_estimates <- function(outcome, 
                                groups, 
                                treatment, 
                                data, 
                                weights){
  
  ## DEFINE OBJECTS NEEDED FOR FUNCTION ##
  out <- list()
  
  groups <- dimnames(weights)[[1]]
  alphas <- as.numeric(dimnames(weights)[[length(dim(weights))]])
  trt_lvls <- sort(unique(data[, treatment]))
  
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
  
  ## CALCULATE MARGINAL ESTIMATES ####
  ybar <- group_means(Y = outcome, A = treatment, G = groups, a = NA, data = data)
  
  grp_est <- apply(weights, 2:3, function(x) x * ybar) 
  dimnames(grp_est) <- list(groups, predictors, alphas)
  
  oa_est <- apply(grp_est, 2:3, sum, na.rm = TRUE)/N

  out$marginal_outcomes$groups <- drop(grp_est)
  out$marginal_outcomes$overall <- drop(oa_est)
  
  ## CALCULATE OUTCOME ESTIMATES PER TREATMENT LEVEL####
  
  hold_grp <- array(dim = c(N, p, k, l), dimnames = list(groups, predictors, 
                                                     alphas, trt_lvls))
  hold_oal <- array(dim = c(p, k, l),
                    dimnames = list(predictors, alphas, trt_lvls))
  
  for(ll in 1:l){    
    a <- trt_lvls[ll]
    
    # Compute means per group
    ybar_trt <- group_means(Y = outcome, A = treatment, G = groups, a = a, data = data)
    
    # Modify weights per treatment level
    weights_trt <- array(dim= c(N, p, k))
    
    for(pp in 1:p){
      weights_trt[ , pp, ] <- t(t(weights[ , pp, ])/((alphas^a * (1-alphas)^(1-a))))
    }
    
    # Compute estimates
    grp_est <- apply(weights_trt, 2:3, function(x) x * ybar_trt) 
    oal_est <- apply(grp_est, 2:3, sum, na.rm = TRUE)/N
    
    hold_grp[ , , , ll] <- grp_est
    hold_oal[ , , ll]   <- oal_est
  }
  
  out$outcomes <- list(groups = drop(hold_grp), 
                       overall = drop(hold_oal))

  ## DONE ####
  return(out)
}