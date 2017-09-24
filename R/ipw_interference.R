#-----------------------------------------------------------------------------#
# IPW Interference estimation 
#
# Prepares the object necessary to compute IPW effect estimates with 
# \code{\link{ipw_effect_calc}}.
# 
# @inheritParams interference
# @param Y outcome vector
# @param X covariate matrix
# @param A treatmeent vector
# @param B 'participation' vector. Defaults to A in the case there is no
# participation variable.
# @param G group assignment vector
# @param parameters a list of fixed_effects and random_effects
# @param variance_estimation currently supports 'robust' or 'naive'
# @param ... additional arguments passed to other functions such as 
# \code{\link{glmer}}, \code{\link{grad}}, and \code{integrand} or \code{likelihood}.
# @return Returns a list of overall and group-level IPW point estimates 
# (the output of \code{\link{ipw_point_estimates}}), overall and group-level IPW 
# point estimates (using the weight derivatives), scores (the output of 
# \code{\link{score_matrix}}), the computed weight matrix, and the computed 
# weight derivative array.
# @export
#-----------------------------------------------------------------------------#

ipw_interference <- function(propensity_integrand,
                             loglihood_integrand = propensity_integrand,
                             allocations,
                             Y, X, A, B = A, G, 
                             parameters,
                             variance_estimation,
                             runSilent   = TRUE, 
                             integrate_allocation,
                             ...)
{
  dots <- list(...)
  
  ## Warnings ##
  
  #### Arguments Necessary for Causal Estimation Functions ####
  integrand_args <- get_args(FUN = propensity_integrand, args_list = dots)
  point_est_args <- get_args(FUN = ipw_point_estimates, args_list = dots)
  loglihood_args <- get_args(FUN = loglihood_integrand, args_list = dots)
  grad_args      <- get_args(FUN = numDeriv::grad, args_list = dots)
  integrate_args <- get_args(FUN = stats::integrate, args_list = dots)
  
  weight_args <- append(append(integrand_args, integrate_args),
                        list(integrand   = propensity_integrand, 
                             allocations = allocations, 
                             X = X, A = A, G = G,
                             parameters = parameters,
                             runSilent  = runSilent, #BB 2015-06-23
                             integrate_allocation = integrate_allocation
                             ))
  #### Prepare output ####
  out <- list()  

  ## Compute Weights ##
  weights <- do.call(wght_matrix, args = weight_args)
  
  if(variance_estimation == 'robust'){
    weightd <- do.call(wght_deriv_array, args = append(weight_args, grad_args)) 
    out$weightd <- weightd
  }
  
  
  #### COMPUTE ESTIMATES AND OUTPUT ####
  estimate_args <- append(point_est_args, list(Y = Y, G = G, A = A))
  point_args    <- append(estimate_args, list(weights = weights))


  #### Calculate output ####
  out$point_estimates <- do.call(ipw_point_estimates, args = point_args)
  
  if(variance_estimation == 'robust'){
    U_args     <- append(estimate_args, list(weights = weightd))
    sargs      <- append(append(loglihood_args, grad_args), integrate_args)
    score_args <- append(sargs, list(integrand = loglihood_integrand,
                                     X = X, G = G, 
                                     A = B, # Use B for treatment in scores
                                     parameters = parameters,
                                     runSilent  = runSilent #BB 2015-06-23
                                     ))
    
    # set randomization scheme to 1 for scores for logit_integrand
    score_args$randomization <- 1
    
    out$Upart           <- do.call(ipw_point_estimates, args = U_args)
    out$scores          <- do.call(score_matrix, args = score_args)
  } 

  out$weights <- weights
 # out$variance_estimation <- variance_estimation #for use in ipw_effect_calc()
  
  return(out)
}