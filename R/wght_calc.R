#' Compute IPW weight
#' 
#' Calculates the IPW for a single group. Used by \code{\link{wght_matrix}} to 
#' create a matrix of weights for each group and allocation scheme.
#' 
#' Type c =
#' \deqn{\frac{\prod_{j=1}^n \alpha^A_j (1 - \alpha)^(1- 
#' A_j)}{Pr(A|X)}}{prod(alpha^A(1 - alpha)^A) / integrate(PrAX_integrand)}
#' Type b =
#' \deqn{\frac{1}{Pr(A|X)}}{1 / integrate(PrAX_integrand)}
#' 
#' Type b incorporates the numerator of type c into the \eqn{Pr(A|X)} integral, 
#' resulting in a slower computation but more accurate results (especially for large groups)
#' 
#' @param f.ab function to pass to the argument 'f' of \code{\link{integrate}}.
#' @param type see description
#' @param ... other arguments passed to f.ab. If type = 'c', then arguments A 
#' and alpha must be defined here
#' @return scalar result of the integral
#' @export

wght_calc <- function(f.ab, type = NULL, ...){  
  
  f.ab <- match.fun(f.ab)
  dots <- list(...)
  args <- append(get_args(f.ab, ...), list(f = f.ab, lower = -Inf, upper = Inf))
  
  if("type" %in% names(formals(f.ab))){
    args$type <- type
  }
  
  # if any of the products within the integrand return Inf, then return NA
  # else return the result of integration
  
  f <- try(do.call("integrate", args = args))
  PrA <- ifelse(is(f, 'try-error'), NA, f$value)
  
  if (type == 'c'){
    if((!'A' %in% names(dots)) | (!'alpha' %in% names(dots))){
      stop("If using type 'c', A and alpha arguments must both be specified")
    }
    A <- dots[[match.arg('A', names(dots))]]
    alpha <- dots[[match.arg('alpha', names(dots))]]
    
    pp <- prod(alpha^A * (1-alpha)^(1-A))
    weight <- pp/PrA
  }
  else if(type == 'b'){
    weight <- 1/PrA
  }
  return(weight)
}