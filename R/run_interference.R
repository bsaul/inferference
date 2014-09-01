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
#' @param propensityB. Defaults to \code{treatment}
#' @param family passed to \code{\link{glmer}}. See \code{\link{family}}
#' @param known_params if the parameter vector is known (e.g. for simulated data), 
#' this argument can be used to pass the known parameters to \code{\link{wght_calc}} 
#' and related functions. If this argument is given, fixed effect parameter estimation
#' is skipped. Defaults to NULL.
#' @param additional arguments passed to other functions. See \code{\link{glmer}}
#' @return Returns a list of ..
#' @export
#-----------------------------------------------------------------------------#

run_interference <- function(f.ab = logit_integrand,
                             likelihood = logit_integrand,
                             alphas,
                             data,
                             groups,
                             outcome,
                             treatment, 
                             predictors,
                             include.alpha = FALSE,
                             propensityB = treatment,
                             family = binomial(link = 'logit'),
                             known_params = NULL,
                             set.NA.to.0 = TRUE,
                             ...){
  ## Necessary bits ##
  dots <- list(...)
  f.ab <- match.fun(f.ab)
  likelihood <- match.fun(likelihood)
  
  ## Reorder data frame by groups ##
  data <- data[order(data[ , groups]), ]
  
  #### Unless parameter values are provided, use glmer() for estimates ####
  form <- paste(propensityB, '~',
                paste(predictors, collapse=' + '), 
                '+ (1|', groups, ')')
  
  if (is.null(known_params)){
    glmer_args <- append(list(formula = form, data = data, family = family),
                         get_args(FUN = glmer, args_list = dots))
    
    fit <- do.call(glmer, args = glmer_args)
    
    theta_fit <- c(fixef(fit), random.eff = log(sqrt(VarCorr(fit)[groups][[1]])))
  } else {
    theta_fit <- known_params
  }

  #### COMPUTE NECESSARY PIECES FOR ESTIMATION ####
  f.ab <- match.fun(f.ab)
  weight_args <- append(list(f.ab = f.ab, 
                             alphas = alphas, 
                             data = data,
                             groups = groups, 
                             predictors = predictors, 
                             treatment = treatment,
                             params = theta_fit, 
                             include.alpha = include.alpha),
                        get_args(FUN = f.ab, args_list = dots))

  weights <- do.call(wght_matrix, args = weight_args)
  weightd <- do.call(wght_deriv_array, args = weight_args)
  
  ## replace any missing weights with 0 ##
  if(set.NA.to.0 == TRUE) {
    weights[is.na(weights)] <- 0
    weightd[is.na(weightd)] <- 0
  }
  
  #### COMPUTE ESTIMATES AND OUTPUT ####
  estimate_args <- append(list(y = outcome, 
                               G = groups, 
                               A = treatment, 
                               data = data),
                          get_args(FUN = ipw_point_estimates, args_list = dots))
  args1 <- append(estimate_args, list(weights = weights))
  args2 <- append(estimate_args, list(weights = weightd))
  
  score_args <- append(list(integrand = likelihood,
                            predictors = predictors,
                            treatment = propensityB, 
                            groups = groups,
                            params = theta_fit, 
                            data = data),
                        get_args(FUN = likelihood, args_list = dots))
  score_args$r <- 1 # set randomization scheme to 1 for scores

  #### Prepare output ####
  out <- list()  

  #### Calculate output ####
  out$point_estimates <- do.call(ipw_point_estimates, args = args1)
  out$Upart   <- do.call(ipw_point_estimates, args = args2)
  out$scores  <- do.call(score_matrix_calc, args = score_args)
  ## replace any Bscores with 0 ##
  if(set.NA.to.0 == TRUE) {
    out$scores[is.na(out$scores)] <- 0
  }
  out$weights <- weights
  out$weightd <- weightd
  
  ## Summary ##
  trt_lvls <- sort(unique(data[, treatment]))
  N <- length(unique(data[ , groups]))
  k <- length(alphas)
  l <- length(trt_lvls)
  weights_na <- apply(weights, 2, function(x) sum(is.na(x)))
  scores_na <- apply(out$scores, 2, function(x) sum(is.na(x)))
  
  out$summary <- list(formula     = ifelse(is.null(known_params), form, NA),
                      ngroups     = N, 
                      nalphas     = k,
                      alphas      = alphas,
                      ntreatments = l,
                      treatments  = trt_lvls,
                      predictors  = predictors,
                      weights_na_count  = weights_na,
                      scores_na_count   = scores_na,
                      parameters  = theta_fit,
                      oracle      = !is.null(known_params))  
  
  print('Run_interference complete')
  return(out)
}