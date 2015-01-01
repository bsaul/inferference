#-----------------------------------------------------------------------------#
#' Estimate Causal Effects in presence of interference  
#'
#' 
#' @param estimation_type currently only supports `ipw`.
#' @param allocations the allocation (coverage) schemes for which to estimate effects. 
#' Must be in (0, 1].
#' @param data the analysis data.frame
#' @param groups quoted name of group variable in \code{data}
#' @param outcome quoted name of outcome variable in \code{data}
#' @param treatment quoted name of treatment variable in \code{data}
#' @param propensityB Optional. If the treatment variable is different for 
#' \code{integrand} and \code{likelihood}, include the quote name of the variable
#' to be passed to \code{likelihood}. Defaults to \code{treatment}. 
#' @param propensity_formula formula passed to model_method for estimation of propensity score 
#' parameters.
#' @param rescale.factor factor by which to rescale values. Defaults to 1.
#' @param conf.level Confidence level for confidence intervals. Defaults to 0.95.
#' @param ... additional arguments passed to other functions such as 
#' \code{\link{glmer}}, \code{\link{grad}}, and \code{integrand} or \code{likelihood}.
#' @return Returns a list of overall and group-level IPW point estimates 
#' (the output of \code{\link{ipw_point_estimates}}), overall and group-level IPW 
#' point estimates (using the weight derivatives), scores (the output of 
#' \code{\link{score_matrix_calc}}), the computed weight matrix, the computed 
#' weight derivative array, and a summary.
#' @export
#-----------------------------------------------------------------------------#

interference <- function(estimation_type = 'ipw',
                         allocations,
                         data,
                         groups,
                         outcome,
                         treatment, 
                         propensityB = treatment,
                         propensity_formula,
                         conf.level = 0.95,
                         rescale.factor = 1,   
                         ...)
{
  dots <- list(...)
  
  # Will have other estimation types in the future
  if(estimation_type == 'ipw')
  {
    ipw_args <- append(dots,
                       list(data = data, 
                            allocations = allocations,
                            groups = groups,
                            outcome = outcome,
                            treatment = treatment,
                            propensityB = propensityB,
                            propensity_formula = propensity_formula))
    out <- do.call(ipw_interference, args = ipw_args)
  }

  grid <- effect_grid(allocations = out$summary$allocations, 
                      treatments  = out$summary$treatments)
  
  out$estimates <- cbind(grid, t(ipw_effect_calc(out, 
                                                 alpha1      = grid$alpha1,
                                                 trt.lvl1    = grid$trt1,
                                                 alpha2      = grid$alpha2,
                                                 trt.lvl2    = grid$trt2,
                                                 marginal    = grid$marginal,
                                                 effect_type = grid$effect_type,
                                                 rescale.factor = rescale.factor,
                                                 conf.level = conf.level,
                                                 print = FALSE)))
  out$summary$conf.level <- conf.level 
  
  class(out) <- "interference"
  return(out)
}