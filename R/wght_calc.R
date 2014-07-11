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
#' 
#' @param type see description
#' @param A vector of treatments
#' @param alpha allocation strategy
#' @param ... other arguments passed to \code{\link{PrAX_integrand}}
#' @return scalar result of the integral
#' @export

wght_calc <- function(type, A, alpha, ...){
  
  # if any of the products within the integrand return Inf, then return NA
  # else return the result of integration

  f <- try(integrate(PrAX_integrand, -Inf, Inf, type = type, 
                     alpha = alpha, A = A, ...))
  PrA <- ifelse(is(f, 'try-error'), NA, f$value)
  
  if (type == 'c'){
    pp <- prod(alpha^A * (1-alpha)^(1-A))
    weight <- pp/PrA
  }
  else if(type == 'b'){
    weight <- 1/PrA
  }
  return(weight)
}