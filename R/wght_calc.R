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
#' resulting in a slower computation but more accurate results (especially 
#' for large groups)
#' 
#' @param f.ab function to pass to the argument 'f' of \code{\link{integrate}}.
#' @param type see description. Defaults to 'b'.
#' @param x necessary argument for \code{\link{grad}}. Defaults to NA, so if 
#' not evaluting a derivative with \code{\link{wght_deriv_calc}}, this can be
#' ignored.
#' @param pos necessary taking a derivative when using to \code{\link{PrAX_integrand}} 
#' Defaults to NA, so if not evaluting a derivative with \code{\link{wght_deriv_calc}}, 
#' this can be ignored.
#' @param ... other arguments passed to f.ab. If type = 'c', then arguments A 
#' and alpha must be defined here
#' @return scalar result of the integral
#' @export

wght_calc <- function(f.ab, 
                      type = 'b',
                      x = NA, 
                      pos = NA, 
                      alpha,
                      ...)
{  
  
  f.ab <- match.fun(f.ab)
  f.ab.names <- names(formals(f.ab))
  dots <- list(...)
  dot.names <- names(dots)
  args <- append(get_args(f.ab, ...), 
                 list(f = f.ab, lower = -Inf, upper = Inf))
  
  ## TODO: This could get cleaned up ##
  if("type" %in% f.ab.names){
    args$type <- type
  }
  if("x" %in% f.ab.names){
    args$x <- x
  }
  if("pos" %in% f.ab.names){
    args$pos <- pos
  }
  if("alpha" %in% f.ab.names){
    args$alpha <- alpha
  }
  # END TODO ##
    
  # if any of the products within the integrand return Inf, then return NA
  # else return the result of integration
  
  f <- try(do.call("integrate", args = args))
  PrA <- ifelse(is(f, 'try-error'), NA, f$value)
  
  if (type == 'c'){
    if(!'A' %in% dot.names){
      stop("If using type 'c', A (treatment assignment) arguments must 
           be specified")
    }
    A <- dots[[match.arg('A', dot.names)]]
    #alpha <- dots[[match.arg('alpha', dot.names)]]
    
    pp <- prod(alpha^A * (1-alpha)^(1-A))
    weight <- pp/PrA
  }
  else if(type == 'b'){
    weight <- 1/PrA
  }
  return(weight)
}
