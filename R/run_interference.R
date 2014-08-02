#-----------------------------------------------------------------------------#
#' Run Interference  
#'
#' Prepares the object necessary to compute IPW effect estimates with 
#' \code{\link{calc_effect}}, \code{\link{direct_effect}}, \code{\link{indirect_effect}},
#' \code{\link{total_effect}}, or \code{\link{overall_effect}}.
#' 
#' @param f.ab f.ab function to pass to the argument 'f' of \code{\link{integrate}}.
#' @param alphas the allocation schemes for which to estimate effects. Must be in (0, 1].
#' @param data the analysis data.frame
#' @param groups quoted name of group variable in \code{data}
#' @param outcome quoted name of outcome variable in \code{data}
#' @param treatment quoted name of treatment variable in \code{data}
#' @param predictors character vector of predictor variables in \code{data}
#' @param type see \code{\link{wght_calc}}
#' @param propensityB
#' @param family
#' @param additional arguments passed to other functions. See \code{\link{glmer}}
#' @return Returns a list of ..
#' @export
#-----------------------------------------------------------------------------#

run_interference <- function(f.ab,
                             alphas,
                             data,
                             groups,
                             outcome,
                             treatment, 
                             predictors,
                             type = 'c',
                             rescale.factor = 1,
                             propensityB = treatment,
                             family = binomial,
                             ...){
  #### FIT GLMER MODEL ####
  form <- paste(propensityB, '~',
                paste(predictors, collapse=' + '), 
                '+ (1|', groups, ')')
  
  glmer_args <- append(list(formula = form, data = data, family = family),
                       get_args(FUN = glmer, ...))
  
  fit <- do.call(glmer, args = glmer_args)
  
  theta_fit <- c(fixef(fit), random.var = VarCorr(fit)[groups][[1]])

  #### COMPUTE NECESSARY PIECES FOR ESTIMATION ####
  f.ab <- match.fun(f.ab)
  weight_args <- append(list(f.ab = f.ab, 
                             alphas = alphas, 
                             data = data,
                             groups = groups, 
                             predictors = predictors, 
                             treatment = treatment,
                             theta = theta_fit, 
                             type = type),
                        get_args(FUN = f.ab, ...))
  
  weights <- do.call(wght_matrix, args = weight_args)
  weightd <- do.call(wght_deriv_array, args = weight_args)
  
  #### COMPUTE ESTIMATES AND OUTPUT ####
  estimate_args <- append(list(y = outcome, 
                               G = groups, 
                               A = treatment, 
                               data = data,
                               rescale.factor = rescale.factor),
                          get_args(FUN = ipw_point_estimates, ...))
  args1 <- append(estimate_args, list(weights = weights))
  args2 <- append(estimate_args, list(weights = weightd))
  
  bscore_args <- append(list(predictors = predictors,
                             B = propensityB, 
                             G = groups,
                             theta = theta_fit, 
                             data = data),
                        get_args(FUN = bscore_calc, ...))
 
  out <- list()
  out$point_estimates <- do.call(ipw_point_estimates, args = args1)
  out$Upart  <- do.call(ipw_point_estimates, args = args2)
  out$bscores <- do.call(bscore_calc, args = bscore_args)

  return(out)
}