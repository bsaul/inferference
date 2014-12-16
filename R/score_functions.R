#' Log Likelihood 
#' 
#' Used by \code{\link{score_matrix_calc}} to compute the log likelihood.
#' 
#' @param x used by \code{\link{grad}} to take the derivative of the 
#' \code{\link{integrate}}(\code{integrand}) with respect to each value of the 
#' \code{param} argument in \code{integrand} 
#' @param pos 
#' @param integrand Defaults to logit_integrand
#' @return value of log likelihood
#' @export

log_likelihood <- function(x, 
                           pos, 
                           integrand = logit_integrand, 
                           ...)
{
  ## Necessary pieces ##
  integrand <- match.fun(integrand)
  dots      <- list(...)
  dot.names <- names(dots)
  
  # If include.allocation is used in the integrand (as in logit_integrand)
  # set this argument to FALSE
  if('include.allocation' %in% names(formals(integrand))){
    dots$include.allocation <- FALSE
  }
  
  ## Integrate() arguments ##
  if(!'lower' %in% dot.names){
    dots$lower <- -Inf
  }
  
  if(!'upper' %in% dot.names){
    dots$upper <- Inf
  }
  
  int.args <- append(get_args(integrate, dots),
                     get_args(integrand, dots))
  args <- append(int.args, list(f = integrand, x = x, pos = pos))
  
  ## Calculuation ##
  attempt <- try(do.call(integrate, args = args))
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
#' @param params See \code{\link{logit_integrand}}.
#' @param ... additional arguments pass to the integrand function
#' @return length(theta) vector of scores
#' @export

score_calc <- function(integrand = logit_integrand,
                       hide.errors = TRUE,
                       fixed.effects,
                       random.effect,
                       ...)
{
  ## Necessary bits ##
  params <- c(fixed.effects, random.effect)
  integrand <- match.fun(integrand)
  dots <- list(...)
  
  ## Function arguments ##
  int.args <- append(get_args(integrand, dots),
                     get_args(integrate, dots))
  fargs    <- append(int.args, get_args(grad, dots))
  
  ## Compute the derivative of the log likelihood for each parameter ##
  scores <- sapply(1:length(params), function(i){
    args <- append(fargs,
                   list(func = log_likelihood, 
                        fixed.effects = fixed.effects,
                        random.effect = random.effect,
                        x = params[i], 
                        pos = i))
    
    attempt <- try(do.call(grad, args = args), silent = hide.errors)
    return(ifelse(is(attempt, 'try-error'), NA, attempt))
  })
  
  return(scores)
}