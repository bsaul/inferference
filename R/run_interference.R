#-----------------------------------------------------------------------------#
#' Run Interference  
#'
#' Prepares the object necessary to compute IPW effect estimates with 
#' \code{\link{calc_effect}}, \code{\link{direct_effect}}, \code{\link{indirect_effect}},
#' \code{\link{total_effect}}, or \code{\link{overall_effect}}.
#' 
#' @param integrand function used as the integrand in computing the IPW weights. 
#' Defaults to \code{\link{logit_integrand}}.
#' @param likelihood function used as the integrand in computing the scores in 
#' \code{\link{score_calc}}. Defaults to \code{\link{logit_integrand}}.
#' @param allocations the allocation (coverage) schemes for which to estimate effects. 
#' Must be in (0, 1].
#' @param data the analysis data.frame
#' @param groups quoted name of group variable in \code{data}
#' @param outcome quoted name of outcome variable in \code{data}
#' @param treatment quoted name of treatment variable in \code{data}
#' @param predictors character vector of predictor variables in \code{data}
#' @param propensityB Optional. If the treatment variable is different for 
#' \code{integrand} and \code{likelihood}, include the quote name of the variable
#' to be passed to \code{likelihood}. Defaults to \code{treatment}.
#' @param family passed to \code{\link{glmer}}. See \code{\link{family}}.
#' @param known_params if the parameter vector is known (e.g. for simulated data), 
#' this argument can be used to pass the known parameters to \code{\link{wght_calc}} 
#' and related functions. If this argument is given, fixed effect parameter estimation
#' is skipped. Defaults to NULL.
#' @param additional arguments passed to other functions such as 
#' \code{\link{glmer}}, \code{\link{grad}}, and \code{\integrand} or \code{likelihood}.
#' @return Returns a list of overall and group-level IPW point estimates 
#' (the output of \code{\link{ipw_point_estimates}}), overall and group-level IPW 
#' point estimates (using the weight derivatives), scores (the output of 
#' \code{\link{score_matrix_calc}}), the computed weight matrix, the computed 
#' weight derivative array, and a summary.
#' @export
#-----------------------------------------------------------------------------#

run_interference <- function(integrand = logit_integrand,
                             likelihood = logit_integrand,
                             allocations,
                             data,
                             groups,
                             outcome,
                             treatment, 
                             predictors,
                             propensityB = treatment,
                             family = binomial(link = 'logit'),
                             known_params = NULL,
                             set.NA.to.0 = TRUE,
                             ...){
  ## Necessary bits ##
  dots <- list(...)
  integrand  <- match.fun(integrand)
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
    
    theta_fit <- c(fixef(fit), random.eff = sqrt(VarCorr(fit)[groups][[1]]))
  } else {
    theta_fit <- known_params
  }

  #### COMPUTE NECESSARY PIECES FOR ESTIMATION ####
  integrand <- match.fun(integrand)
  weight_args <- append(list(integrand = integrand, 
                             allocations = allocations, 
                             data = data,
                             groups = groups, 
                             predictors = predictors, 
                             treatment = treatment,
                             params = theta_fit),
                        get_args(FUN = integrand, args_list = dots))

  weights <- do.call(wght_matrix, args = weight_args)
  weightd <- do.call(wght_deriv_array, args = weight_args)
  
  ## replace any missing weights with 0 ##
  if(set.NA.to.0 == TRUE) {
    weights[is.na(weights)] <- 0
    weightd[is.na(weightd)] <- 0
  }
  
  #### COMPUTE ESTIMATES AND OUTPUT ####
  estimate_args <- append(list(outcome   = outcome, 
                               groups    = groups, 
                               treatment = treatment, 
                               data      = data),
                          get_args(FUN = ipw_point_estimates, args_list = dots))
  args1 <- append(estimate_args, list(weights = weights))
  args2 <- append(estimate_args, list(weights = weightd))
  
  sargs      <- append(get_args(FUN = likelihood, args_list = dots),
                       get_args(FUN = grad, args_list = dots))
  score_args <- append(list(integrand = likelihood,
                            predictors = predictors,
                            treatment = propensityB, 
                            groups = groups,
                            params = theta_fit, 
                            data = data),
                        sargs)
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
  k <- length(allocations)
  l <- length(trt_lvls)
  weights_na <- apply(weights, 2, function(x) sum(is.na(x)))
  scores_na <- apply(out$scores, 2, function(x) sum(is.na(x)))
  
  out$summary <- list(formula      = ifelse(is.null(known_params), form, NA),
                      ngroups      = N, 
                      nallocations = k,
                      allocations  = allocations,
                      ntreatments  = l,
                      treatments   = trt_lvls,
                      predictors   = predictors,
                      weights_na_count  = weights_na,
                      scores_na_count   = scores_na,
                      parameters  = theta_fit,
                      oracle      = !is.null(known_params))  
  
  print('Run_interference complete')
  return(out)
}