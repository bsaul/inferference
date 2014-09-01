#' Log Likelihood 
#' 
#' Used by \code{\link{score_matrix_calc}} to compute the log likelihood.
#' 
#' @param x
#' @param pos
#' @param integrand Defaults to logit_integrand
#' @return value of log likelihood
#' @export

log_likelihood <- function(x, pos, integrand = logit_integrand, ...){
  integrand <- match.fun(integrand)
  dots <- list(...)
  args <- append(get_args(integrand, dots), 
                 list(f = integrand, lower = -Inf, upper = Inf, x = x, pos = pos))
  
  attempt <- try(do.call("integrate", args = args))
  val <- ifelse(is(attempt, 'try-error'), NA, attempt$value)

  return(log(val))
}

#' Compute scores for a single group
#' 
#' Used by \code{\link{score_matrix_calc}} to log likelihood derivatives for
#' a single group.
#' 
#' @param integrand function to used for the integrand. 
#' Defaults to \code{\link{logit_integrand}}
#' @param hide.errors Hide errors printed from \code{\link{grad}}. 
#' Defaults to true
#' @param theta See \code{\link{logit_integrand}}.
#' @param ... additional arguments pass to the integrand function
#' @return length(theta) vector of scores
#' @export

score_calc <- function(integrand = logit_integrand,
                       hide.errors = TRUE,
                       params,
                       ...){
  ## Necessary bits ##
  integrand <- match.fun(integrand)
  dots <- list(...)
  
  ## Compute the derivative of the log likelihood for each parameter ##
  scores <- sapply(1:length(params), function(i){
    args <- append(get_args(integrand, dots),
                   list(func = log_likelihood, params = params, x = params[i], pos = i,
                        method = 'simple'))
    attempt <- try(do.call('grad', args = args), silent = hide.errors)
    return(ifelse(is(attempt, 'try-error'), NA, attempt))
  })
  
  return(scores)
}