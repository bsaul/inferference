#' Compute a group weight
#' 
#' Type c =
#' \eqn{\frac{\pi}{Pr(A|X)}}
#' Type b =
#' \eqn{\frac{1}{Pr(A|X)}}
#' where the \eqn{\pi} term is brought into the dominator
#' 
#' @param f the function to integrad. Defaults to \code{\link{PrAX_integrand}}
#' @param llimit lower limit of the integral defaults to \code{-Inf}
#' @param ulimit upper limit of the integral defaults to \code{Inf}
#' @param type see description
#' @param alpha allocation strategy
#' @param A vector of treatments
#' @param ... other arguments passed to integrand
#' @return scalar which is the result of the integral
#' @export

wght_calc <- function(type, A, alpha, hide.errors = TRUE, ...){
  
  # if any of the products within the integrand return Inf, then return NA
  # else return the result of integration
  f <- try(integrate(PrAX_integrand, -Inf, Inf, type = type, 
                     alpha = alpha, A = A, ...),
           silent = hide.errors)
  PrA <- ifelse(is(f, 'try-error'), NA, f$value)
  
  ## TODO: I tried to write this function so that integrate() could have 
  # different integrands passed to it from the arguments. It worked for 
  # computing weights, but was fraught with 'unused argument' errors 
  # when taking derivatives. I couldn't figure out the problem on first blush
  # and moved on.
  
  if (type == 'c'){
    pp <- prod(alpha^A * (1-alpha)^(1-A))
    weight <- pp/PrA
  }
  else if(type == 'b'){
    weight <- 1/PrA
  }
  return(weight)
}