#-----------------------------------------------------------------------------#
#' Estimate Causal Effects in presence of interference  
#'
#' @param formula
#' @param propensity_integrand function used as the integrand in computing the IPW weights. 
#' Defaults to \code{\link{logit_integrand}}.
#' @param loglihood_integrand function used as the integrand in computing the scores in 
#' \code{\link{score_calc}}. Defaults to \code{\link{logit_integrand}}.
#' @param allocations the allocation (coverage) schemes for which to estimate effects. 
#' Must be in [0, 1].
#' @param data the analysis data.frame
#' @param causal_estimation_method currently only supports `ipw`.
#' @param causal_estimation_options list of options passed to the corresponding 
#' estimation function. `ipw` is \code{\link{ipw_interference}}
#' @param rescale.factor factor by which to rescale values. Defaults to 1.
#' @param model_method One of three options: `glmer`, `glm`, or `oracle`
#' @param model_options List of options passed to \code{model_method}
#' @param set_NA_to_0 if TRUE, sets any weights that returned an NA value to 0. 
#' Defaults to TRUE
#' @param conf.level Confidence level for confidence intervals. Defaults to 0.95.
#' @param ... additional arguments passed to other functions such as 
#' \code{\link{glmer}}, \code{\link{grad}}, \code{\link{integrate}}, and \code{integrand} or \code{likelihood}.
#' @return Returns a list of overall and group-level IPW point estimates 
#' (the output of \code{\link{ipw_point_estimates}}), overall and group-level IPW 
#' point estimates (using the weight derivatives), scores (the output of 
#' \code{\link{score_matrix_calc}}), the computed weight matrix, the computed 
#' weight derivative array, and a summary.
#' @export
#-----------------------------------------------------------------------------#

interference <- function(formula,
                         propensity_integrand,
                         loglihood_integrand = propensity_integrand,
                         allocations,
                         data,
                         model_method = "glmer",
                         model_options = list(family = binomial(link = 'logit')),
                         causal_estimation_method = 'ipw',
                         causal_estimation_options = list(set_NA_to_0 = TRUE, 
                                                          variance_estimation = 'robust'),
                         conf.level = 0.95,
                         rescale.factor = 1,   
                         ...)
{
  ## Necessary bits ##
  dots <- list(...)
  integrandFUN    <- match.fun(propensity_integrand)
  likelihoodFUN   <- match.fun(loglihood_integrand)
  oracle          <- model_method == 'oracle'
  cformula        <- Formula::Formula(formula)
  len_lhs         <- length(cformula)[1]
  len_rhs         <- length(cformula)[2]

  
  ## For the sake of consistency downstream, reorder data frame by groups ##
  group_var <- attr(terms(cformula, lhs = 0, rhs = len_rhs), 'term.labels')
  data <- data[order(data[ , group_var]), ]
  
  ## Parse out the formula into necessary pieces ##
  Y <- Formula::model.part(cformula, data = data, lhs = 1, drop = TRUE)
  A <- Formula::model.part(cformula, data = data, lhs = 2, drop = TRUE)
  G <- Formula::model.part(cformula, data = data, rhs = len_rhs, drop = TRUE)
  
  # Used when there is 'participation' variable
  if(len_lhs > 2){
    B <- Formula::model.part(cformula, data = data, lhs = len_lhs, drop = TRUE)
  } else {
    B <- A
  }
  
  propensity_formula <- formula(terms(cformula, lhs = len_lhs, rhs = -2))
  random.count <- length(lme4::findbars(propensity_formula))

  trt_lvls     <- sort(unique(A))
  N            <- length(unique(G))
  k            <- length(allocations)
  l            <- length(trt_lvls)

  ## Warnings ##
  if(model_method == 'glm' & random.count > 0 ){
    stop('propensity_formula appears to include a random effect when "glm" was chosen \n 
         for parameter estimation. Set model_method to "glmer" to include a random effect')
  }
  
  if(propensity_integrand == "logit_integrand" & random.count > 1){
    stop('Logit integrand is designed to handle only 1 random effect.')
  }
  
  if(min(allocations) < 0 | max(allocations) > 1){
    stop('Allocations must be between 0 and 1 (inclusive)')
  }
  
  if(length(allocations) < 2){
    warning('At least 2 allocations must be specified in order to estimate indirect effects')
  }
  
  #### Compute Parameter Estimates ####

  estimation_args <- append(list(formula = propensity_formula, data = data), 
                            model_options)
  
  if(model_method == "glmer"){
    propensity_model <- do.call(lme4::glmer, args = estimation_args)
    fixed.effects  <- lme4::getME(propensity_model, 'fixef')
    random.effects <- lme4::getME(propensity_model, 'theta')
    X <- lme4::getME(propensity_model, "X")
  } else if(model_method == "glm"){
    propensity_model <- do.call("glm", args = estimation_args)
    fixed.effects  <- coef(propensity_model)
    random.effects <- NULL
    X <- model.matrix(propensity_model)
  } else if(model_method == "oracle"){
    fixed.effects  <- model_options[[1]]
    random.effects <- model_options[[2]]
    X <- model.matrix(propensity_formula, data)
    
    if(length(fixed.effects) != ncol(X)){
      stop('oracle fixed effects vector must have length of # of columns of X')
    }
  }
  
  #### Compute Effect Estimates ####
  out <- list()
  grid <- effect_grid(allocations = allocations, treatments  = trt_lvls)
  
  # Will have other causal estimation types in the future
  if('ipw' %in% causal_estimation_method)
  {
    ipw_args <- append(append(dots, causal_estimation_options),
                       list(propensity_integrand = integrandFUN, 
                            loglihood_integrand  = likelihoodFUN,
                            allocations          = allocations,
                            fixed.effects        = fixed.effects, 
                            random.effects       = random.effects,
                            Y = Y, X = X, A = A, B = B, G = G))
  
    ipw <- do.call(ipw_interference, args = ipw_args)
    out <- append(out, ipw)
    
    print('Computing effect estimates...')
    
    out$estimates <- cbind(grid, t(ipw_effect_calc(ipw, 
                                                   alpha1      = grid$alpha1,
                                                   trt.lvl1    = grid$trt1,
                                                   alpha2      = grid$alpha2,
                                                   trt.lvl2    = grid$trt2,
                                                   marginal    = grid$marginal,
                                                   effect_type = grid$effect_type,
                                                   rescale.factor = rescale.factor,
                                                   conf.level = conf.level,
                                                   print = FALSE)))
  }

  
  #### Summary ####
  
  # count missing weights by allocation
  weights_na <- apply(out$weights, 2, function(x) sum(is.na(x)))
  
  if(!oracle){
    out$models <- list(propensity_model = propensity_model)
  } else {
    out$models <- list(propensity_model = model_options) 
  }
  
  out$summary <- list(causal_model = deparse(formula),
                      oracle       = oracle,
                      conf.level   = conf.level,
                      ngroups      = N, 
                      nallocations = k,
                      npredictors  = length(fixed.effects),
                      ntreatments  = l,
                      allocations  = allocations,
                      treatments   = trt_lvls,
                      weights_na_count  = weights_na)
  
  class(out) <- "interference"
  
  print('Interference complete')
  return(out)
}