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
#' @param set.NA.to.0 if TRUE, sets any weights that returned an NA value to 0. 
#' Defaults to TRUE.
#' @param additional arguments passed to other functions such as 
#' \code{\link{glmer}}, \code{\link{grad}}, and \code{integrand} or \code{likelihood}.
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
                             propensityB = treatment,
                             propensity_formula,
                             estimation_method = 'glmer',
                             estimation_options = list(family = binomial(link = 'logit')),
                             set.NA.to.0 = TRUE,
                             ...){
  ## Necessary bits ##
  dots <- list(...)
  integrand  <- match.fun(integrand)
  likelihood <- match.fun(likelihood)
#  est_method <- match.fun(estimation_method)
  oracle     <- estimation_method == 'oracle'
  
  ## Reorder data frame by groups ##
  data <- data[order(data[ , groups]), ]
  
  #### Compute Parameter Estimates ####
  
  estimation_args <- append(list(formula = propensity_formula, data = data), 
                            estimation_options)
    
  if(estimation_method == "glmer"){
    propensity_model <- do.call(glmer, args = estimation_args)
    fixed.effects <- getME(propensity_model, 'fixef')
    random.effect <- getME(propensity_model, 'theta')[1]
    XXp <- getME(propensity_model, "X")
  } else if(estimation_method == "glm"){
    propensity_model <- do.call(glm, args = estimation_args)
    fixed.effects <- coef(propensity_model)
    random.effect <- NULL
    XXp <- model.matrix(propensity_model)
  } else if(estimation_method == "oracle"){
    fixed.effects <- estimation_options[[1]]
    random.effect <- estimation_options[[2]]
    XXp <- model.matrix(propensity_formula, data)
  }
  
  ## Propensity Model Pieces ##
  AAp <- data[, treatment]
  GGp <- data[, groups]
  BBp <- data[, propensityB]

  #### Arguments Necessary for Causal Estimation Functions ####
  integrand_args <- get_args(FUN = integrand, args_list = dots)
  point_est_args <- get_args(FUN = ipw_point_estimates, args_list = dots)
  loglihood_args <- get_args(FUN = likelihood, args_list = dots)
  grad_args      <- get_args(FUN = grad, args_list = dots)
  integrate_args <- get_args(FUN = integrate, args_list = dots)
  
  weight_args <- append(append(integrand_args, integrate_args),
                        list(integrand = integrand, 
                             allocations = allocations, 
                             X = XXp, A = AAp, G = GGp,
                             fixed.effects = fixed.effects,
                             random.effect = random.effect))
  
  ## Compute Weights ##
  weights <- do.call(wght_matrix, args = weight_args)
  weightd <- do.call(wght_deriv_array, args = append(weight_args, grad_args))
  
  ## replace any missing weights with 0 ##
  if(set.NA.to.0 == TRUE) {
    weights[is.na(weights)] <- 0
    weightd[is.na(weightd)] <- 0
  }
  
  #### COMPUTE ESTIMATES AND OUTPUT ####
  estimate_args <- append(point_est_args, list(outcome   = outcome, 
                                               groups    = groups, 
                                               treatment = treatment, 
                                               data      = data))
  point_args <- append(estimate_args, list(weights = weights))
  U_args     <- append(estimate_args, list(weights = weightd))
  sargs      <- append(append(loglihood_args, grad_args), integrate_args)
  score_args <- append(sargs, list(integrand = likelihood,
                                   X = XXp, G = GGp, A = BBp,
                                   fixed.effects = fixed.effects,
                                   random.effect = random.effect))
  
  # set randomization scheme to 1 for scores for logit_integrand
  score_args$randomization <- 1
  # If integrate.allocation is used in the integrand (as in logit_integrand)
  # set this argument to FALSE
  if('integrate.allocation' %in% names(formals(integrand))){
    score_args$integrate.allocation <- FALSE
  }

  #### Prepare output ####
  out <- list()  

  #### Calculate output ####
  out$point_estimates <- do.call(ipw_point_estimates, args = point_args)
  out$Upart           <- do.call(ipw_point_estimates, args = U_args)
  out$scores          <- do.call(score_matrix_calc,   args = score_args)
  
  ## replace any Bscores with 0 ##
  if(set.NA.to.0 == TRUE) {
    out$scores[is.na(out$scores)] <- 0
  }
  out$weights <- weights
  out$weightd <- weightd
  
  #### Summary ####
  trt_lvls <- sort(unique(data[, treatment]))
  N <- length(unique(data[ , groups]))
  k <- length(allocations)
  l <- length(trt_lvls)
  weights_na <- apply(weights, 2, function(x) sum(is.na(x)))
  scores_na <- apply(out$scores, 2, function(x) sum(is.na(x)))

  if(!oracle){
    out$models <- list(propensity_model = propensity_model)
  } else {
    out$oracle_parameters <- estimation_options
  }
  
  out$summary <- list(oracle      = oracle,
                      ngroups      = N, 
                      nallocations = k,
                      allocations  = allocations,
                      ntreatments  = l,
                      treatments   = trt_lvls,
                      predictors   = dimnames(XXp)[[2]],
                      weights_na_count  = weights_na,
                      scores_na_count   = scores_na)  
  
  print('Run_interference complete')
  return(out)
}