#-----------------------------------------------------------------------------#
#' Compute IPW weight
#' 
#' Calculates the IPW for a single group. Used by \code{\link{wght_matrix}} to
#' create a matrix of weights for each group and allocation scheme.
#' 
#' If \code{integrate.allocation} is an argument in the integrand function and
#' \code{integrate.allocation == TRUE}, then the weight is calcuated as:
#' 
#' \deqn{\frac{1}{Pr(A|X)}}{1 / integrate(integrand)}
#' 
#' Otherwise, the weight is computed by:
#' \deqn{\frac{\prod_{j=1}^n \alpha^A_j (1 - \alpha)^(1-
#' A_j)}{Pr(A|X)}}{prod(allocation^A(1 - allocation)^A) / integrate(integrand)}
#' 
#' @param parameters vector of parameter values
#' @param integrand function to pass to the argument 'f' of \code{\link{integrate}}.
#' @param allocation the allocation ratio for which to compute the weight
#' @param ... other arguments passed to integrand.
#' @return scalar result of the integral
#' @export
#'
#-----------------------------------------------------------------------------#

wght_calc <- function(parameters,
                      integrand, 
                      allocation,
                      ...)
{  
  ## Necessary pieces ##
  integrand         <- match.fun(integrand)
  integrand.formals <- names(formals(integrand))
  dots              <- list(...)
  dot.names         <- names(dots)
  A                 <- dots[[match.arg('A', dot.names)]]
  
  ## Warnings ##
  if(!'A' %in% dot.names){
    stop("The argument 'A' (treatment assignment) must be specified")
  }
  
  ## Integrate() arguments ##
  if(!'lower' %in% dot.names){
    dots$lower <- -Inf
  }
  
  if(!'upper' %in% dot.names){
    dots$upper <- Inf
  }
  
  int.args <- append(get_args(stats::integrate, dots),
                     list(f = integrand, parameters = parameters))
  
  args <- append(get_args(integrand, dots), int.args)
  
  # Allocation is optional in user-defined integrands. Include this argument 
  # when necessary. Note that allocation will either be used in this function
  # or passed to the integrand function.
  if("allocation" %in% integrand.formals){
    args$allocation <- allocation
  }
  
  ## Compute the integral ##
  # if any of the products within the integrand return Inf, then return NA
  # else return the result of integration
  
  f <- try(do.call("stats::integrate", args = args), silent = TRUE)
  PrA <- if(methods::is(f, 'try-error')) NA else f$value

  ## Compute the weight ##
  weight <- 1/PrA

  return(weight)
}
