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
#' @return Returns a list of ..
#' @export
#-----------------------------------------------------------------------------#

run_interference <- function(f.ab,
                             alphas,
                             data,
                             groups,
                             outcome,
                             treatment, 
                             B,
                             predictors,
                             type = 'c',
                             rescale.factor = 1,
                             ...)
{
  ## FIT GLMER MODEL ##
  form <- paste(B, '~',
                paste(predictors, collapse=' + '), 
                '+ (1|', groups, ')')
  fit <- glmer(form, data = data, family = binomial)
  
  theta_fit <- c(fixef(fit), random.var = VarCorr(fit)[groups][[1]])

  ## NECESSARY PIECES FOR ESTIMATION ##

  WT_fit <- wght_matrix(f.ab = f.ab, alphas = alphas, data = data,
                        groups = groups, predictors = predictors, A = treatment,
                        theta = theta_fit, type = type)

  WTa_fit <- wght_deriv_array(f.ab = f.ab, alphas = alphas, data = data,
                              groups = groups, predictors = predictors, 
                              A = treatment,
                              theta = theta_fit, type = type)

  bscores_fit <- bscore_calc(predictors = predictors, B = B, G = groups,
                             theta = theta_fit, data = data)

  ## GET ESTIMATES ##

  out <- ipw_estimates(y = outcome, 
                       G = groups, 
                       A = treatment, 
                       B = B, 
                       data = data, 
                       weights = WT_fit, 
                       weight_dervs = WTa_fit, 
                       bscores = bscores_fit,
                       predictors = predictors, 
                       rescale.factor = rescale.factor)
  return(out)
}