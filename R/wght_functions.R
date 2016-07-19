#-----------------------------------------------------------------------------#
#' Compute IPW weight
#' 
#' Calculates the IPW for a single group. Used by \code{\link{wght_matrix}} to
#' create a matrix of weights for each group and allocation scheme.
#' 
#' If \code{allocation} is an argument in the integrand function and
#' \code{integrate_allocation == TRUE}, then the weight is calcuated as:
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
#' @param integrate_allocation Indicator of whether the integrand function uses 
#' the allocation parameter. Defaults to TRUE.
#' @param ... other arguments passed to integrand.
#' @return scalar result of the integral
#' @export
#' @importFrom methods is
#'
#-----------------------------------------------------------------------------#

wght_calc <- function(parameters,
                      integrand, 
                      allocation,
                      integrate_allocation = TRUE,
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
  
  f <- try(do.call(stats::integrate, args = args), silent = TRUE)
  PrA <- if(is(f, 'try-error')) NA else f$value
  
  ## Compute the weight ##
  # weight <- 1/PrA
  
  if(integrate_allocation == TRUE){
    weight <- 1/PrA
  } else {
    ppp    <- prod(allocation^A * (1-allocation)^(1-A))
    weight <- ppp/PrA
  } 
  
  weight
}

#-----------------------------------------------------------------------------#
#' Compute the derivative(s) of a weight
#' 
#' Takes the derivative of the \code{\link{wght_calc}} function with respect to each
#' parameter in \code{params}.
#' 
#' @inheritParams wght_calc
#' @return vector of derivatives with respect to element of params
#' @export
#-----------------------------------------------------------------------------#

wght_deriv_calc <- function(integrand,
                            parameters,
                            allocation,
                            integrate_allocation = TRUE,
                            ...)
{  
  ## Necessary pieces ##
  integrand <- match.fun(integrand)
  dots <- list(...)
  
  ## Integrand and  arguments ##
  int.args <- append(get_args(integrand, dots),
                     get_args(stats::integrate, dots))
  
  args <- append(append(int.args, get_args(numDeriv::grad, dots)),
                 list(func       = wght_calc, 
                      integrand  = integrand, 
                      allocation = allocation,
                      x          = parameters))
  
  dervs <- do.call(numDeriv::grad, args = args)
  
  dervs
}

#-----------------------------------------------------------------------------#
# Create a matrix of group IP weights 
#' 
#' Creates a number of groups by number of allocation schemes matrix of group weights.
#' Allocation schemes are selected by the user.
#' 
#' Groups should be numbered 1, ..., N
#' 
#' @param allocations coverage levels in [0, 1]. Can be vector.
#' @param X covariate matrix
#' @param A vector of treatment assignments
#' @param G vector of group assignments
#' @param parameters vector of parameters passed to \code{integrand}
#' @param runSilent if FALSE, errors are printed to console. Defaults to TRUE.
#' @inheritParams wght_calc
#' @return a length(unique(group)) X length(alphas) matrix of group weights
#' @export
#
#-----------------------------------------------------------------------------#

wght_matrix <- function(integrand, 
                        allocations, 
                        X, A, G,
                        parameters,
                        runSilent = TRUE, 
                        integrate_allocation = TRUE,
                        ...)
{
  ## Gather necessary bits ##
  p  <- ncol(X)
  aa <- sort(allocations) # Make sure alphas are sorted
  gg <- sort(unique(G))
  
  ## Compute weight for each group and allocation level ##
  if(!runSilent) print('Calculating matrix of IP weights...') 
  
  w.list <- lapply(aa, function(allocation){
    w <- by(cbind(X, A), INDICES = G, simplify = FALSE, 
            FUN = function(x) {
              wght_calc(parameters = parameters,
                        integrand  = integrand, 
                        allocation = allocation, 
                        integrate_allocation = integrate_allocation,
                        A = x[, p+1], X = x[, 1:p], ...)})
    as.numeric(w)
  }) 
  
  ## Reshape list into matrix ##
  w.matrix <- matrix(unlist(w.list, use.names = FALSE), 
                     ncol = length(allocations), 
                     byrow = FALSE,
                     dimnames = list(gg, aa))
  
  return(w.matrix)
}

#-----------------------------------------------------------------------------#
#' Create an array of group weight derivatives
#' 
#' Uses \code{\link{wght_deriv_calc}} to compute the weight derivatives for each
#' group per coverage level
#' 
#' @inheritParams wght_matrix
#' @inheritParams wght_calc
#' @return a length(unique(group)) X length(params) X length(alphas) array of
#' group weight derivatives
#' @export
#-----------------------------------------------------------------------------#

wght_deriv_array <- function(parameters, 
                             integrand, 
                             allocations, 
                             X, A, G,
                             runSilent = TRUE, 
                             integrate_allocation = TRUE,
                             ...)
{
  ## Gather necessary bits ##
  integrand <- match.fun(integrand)
  XX <- cbind(X, A)
  p  <- ncol(X) # number of predictors
  pp <- length(parameters)
  aa <- sort(allocations) # Make sure alphas are sorted
  gg <- sort(unique(G))
  k  <- length(allocations) 
  N  <- length(unique(G))
  dots <- list(...)
  
  ## Warnings ##
  
  ## Compute weight (derivative) for each group, parameter, and alpha level ##
  if(!runSilent) print('Calculating array of IP weight derivatives...')
  
  w.list <- lapply(aa, function(allocation){
    w <- by(XX, INDICES = G, simplify = TRUE, 
            FUN = function(x) {
              wght_deriv_calc(parameters = parameters,
                              integrand  = integrand, 
                              allocation = allocation, 
                              integrate_allocation = FALSE,
                              A = x[, p+1], 
                              X = x[, 1:p], 
                              ...)})
    w2 <- matrix(unlist(w, use.names = FALSE), ncol = pp, byrow = TRUE)
    return(w2)}) 
  
  ## Reshape list into array ##
  out <- array(unlist(w.list, use.names = FALSE), 
               dim = c(N, pp, k),
               dimnames = list(gg, names(parameters), aa))
  
  return(out)
}