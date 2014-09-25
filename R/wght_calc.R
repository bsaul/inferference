#' Compute IPW weight
#' 
#' Calculates the IPW for a single group. Used by \code{\link{wght_matrix}} to 
#' create a matrix of weights for each group and allocation scheme.
#' 
#' If \code{include.allocation} is an argument in the integrand function and 
#' \code{include.allocation == TRUE}, then the weight is calcuated as:
#' 
#' \deqn{\frac{1}{Pr(A|X)}}{1 / integrate(integrand)}
#' 
#' Otherwise, the weight is computed by:
#' \deqn{\frac{\prod_{j=1}^n \alpha^A_j (1 - \alpha)^(1- 
#' A_j)}{Pr(A|X)}}{prod(allocation^A(1 - allocation)^A) / integrate(integrand)}
#' 
#' @param integrand function to pass to the argument 'f' of \code{\link{integrate}}.
#' @param allocation the allocation ratio for which to compute the weight
#' @param x necessary argument for \code{\link{grad}}. Defaults to NULL, so if 
#' not evaluting a derivative with \code{\link{wght_deriv_calc}}, this can be
#' ignored.
#' @param pos necessary taking a derivative when using to \code{\link{PrAX_integrand}} 
#' Defaults to NULL, so if not evaluting a derivative with \code{\link{wght_deriv_calc}}, 
#' this can be ignored.
#' @param ... other arguments passed to integrand. 
#' @return scalar result of the integral
#' @export

wght_calc <- function(integrand, 
                      allocation,
                      x = NULL, 
                      pos = NULL, 
                      ...){  

  # Necessary pieces
  integrand         <- match.fun(integrand)
  integrand.formals <- names(formals(integrand))
  dots              <- list(...)
  dot.names         <- names(dots)
  A                 <- dots[[match.arg('A', dot.names)]]
  pp                <- prod(allocation^A * (1-allocation)^(1-A))
  
  # Warnings #
  if(!'A' %in% dot.names){
    stop("The argument 'A' (treatment assignment) must be specified")
  }
  
  # Integration arguments
  args <- append(get_args(integrand, dots), 
                 list(f = integrand, lower = -Inf, upper = Inf,
                      x = x, pos = pos))
  
  # Allocation is optional in user-defined integrands
  # include this arguments when necessary: 
  if("allocation" %in% integrand.formals){
    args$allocation <- allocation
  }
  
  ## Compute the integral
  # if any of the products within the integrand return Inf, then return NA
  # else return the result of integration
  
  f <- try(do.call("integrate", args = args), silent = TRUE)
  PrA <- ifelse(is(f, 'try-error'), NA, f$value)
  

  # Compute the weight
  # TODO: the if-else logic is clunky. Basically, there are 3 options:
  # 1) include.allocation is not included in the formals of a user-defined 
  # integrand 2) include.alloction is set to FALSE # 3) include.allocation 
  # is set to TRUE.

  if(!'include.allocation' %in% dot.names){
    weight <- pp/PrA
  } else {
    if(args$include.allocation == TRUE){
      # In this case the pp term is included in PrA (if using logit_integrand)
      weight <- 1/PrA
      } else {
      weight <- pp/PrA
    }
  }
  return(weight)
}
