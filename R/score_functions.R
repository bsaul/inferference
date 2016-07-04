#-----------------------------------------------------------------------------#
#' Log Likelihood
#' 
#' Used by \code{\link{score_matrix}} to compute the log likelihood.
#' 
#' @param x used by \code{\link{grad}} to take the derivative of the
#' \code{\link{integrate}}(\code{integrand}) with respect to each value of the
#' \code{param} argument in \code{integrand}
#' @param pos The position of theta for which to take the derivative.
#' @param integrand Defaults to logit_integrand
#' @param ... additional arguments passed to \code{integrand} function.
#' @return value of log likelihood
#' @export
#-----------------------------------------------------------------------------#

log_likelihood <- function(parameters,
                           integrand,
                           ...)
{
  ## Necessary pieces ##
  integrand <- match.fun(integrand)
  dots      <- list(...)
  dot.names <- names(dots)
  
  ## Integrate() arguments ##
  if(!'lower' %in% dot.names){
    dots$lower <- -Inf
  }
  
  if(!'upper' %in% dot.names){
    dots$upper <- Inf
  }
  
  int.args <- append(get_args(integrate, dots),
                     get_args(integrand, dots))
  args <- append(int.args, list(f = integrand, parameters = parameters))
  
  ## Calculuation ##
  attempt <- try(do.call(integrate, args = args))
  val <- if(is(attempt, 'try-error')) NA else attempt$value

  return(log(val))
}

#-----------------------------------------------------------------------------#
#' Compute scores for a single group
#' 
#' Used by \code{\link{score_matrix}} to log likelihood derivatives for
#' a single group.
#' 
#' @param integrand function to used for the integrand.
#' Defaults to \code{\link{logit_integrand}}.
#' @param hide.errors Hide errors printed from \code{\link{grad}}.
#' Defaults to true.
#' @param fixed.effects vector of fixed effect parameters.
#' @param random.effects OPTIONAL vector random effect parameters.
#' @param ... additional arguments pass to the integrand function.
#' @return length(theta) vector of scores
#' @export
#-----------------------------------------------------------------------------#


score_calc <- function(parameters,
                       integrand,
                       hide.errors = TRUE,
                       ...)
{
  ## Necessary bits ##
  integrand <- match.fun(integrand)
  dots <- list(...)
  
  ## Function arguments ##
  int.args <- append(get_args(integrand, dots),
                     get_args(integrate, dots))
  fargs    <- append(int.args, get_args(numDeriv::grad, dots))
  
  args     <- append(fargs,
                 list(func = log_likelihood,
                      x    = parameters,
                      integrand = integrand))
  ## Compute the derivative of the log likelihood for each parameter ##
  do.call(numDeriv::grad, args = args)
}  
