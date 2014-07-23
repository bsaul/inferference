#-----------------------------------------------------------------------------#
#' Run Interference  
#'
#' @param predictors character vector of predictors from model
#' @param B character naming B variable in data
#' @param G character vector (length == 1) naming group variable in data
#' @param theta
#' @param data dataframe
#' @return N X length(theta) matrix of scores
#' @export
#-----------------------------------------------------------------------------#

run_interference <- function(f.ab,
                             alphas,
                             data,
                             groups,
                             outcome,
                             A, 
                             B,
                             predictors,
                             type = 'b',
                             rescale.factor = 1)
{
  ## FIT GLMER MODEL ##
  form <- paste(B, '~',
                paste(predictors, collapse=' + '), 
                '+ (1|', groups, ')')
  fit <- glmer(form, data = sim1, family = binomial)
  
  theta_fit <- c(fixef(fit), random.var = VarCorr(fit)[groups][[1]])

  ## NECESSARY PIECES FOR ESTIMATION ##

  WT_fit <- wght_matrix(f.ab = f.ab, alphas = alphas, data = data,
                        groups = groups, predictors = predictors, A = A,
                        theta = theta_fit, type = type)

  WTa_fit <- wght_deriv_array(f.ab = f.ab, alphas = alphas, data = data,
                              groups = groups, predictors = predictors, A = A,
                              theta = theta_fit, type = type)

  bscores_fit <- bscore_calc(predictors = predictors, B = B, G = groups,
                             theta = theta_fit, data = data)

  ## GET ESTIMATES ##

  out <- ipw_estimates(y = outcome, G = groups, A = A, B = B, 
                       data = data, 
                       weights = WT_fit, 
                       weight_dervs = WTa_fit, 
                       bscores = bscores_fit,
                       predictors = predictors, 
                       rescale.factor = rescale.factor)
  return(out)
}