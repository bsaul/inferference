#' Compute the derivative(s) of a weight 
#' 
#' @param theta
#' @return vector of derivatives with respect to element of theta
#' @export

wght_deriv_calc <- function(theta, A, X, type, alpha){
  dervs <- sapply(1:length(theta), function(i){
    f <- try(grad(wght_calc, x = theta[i], pos = i, theta = theta, 
                  type = type, alpha = alpha, A = A, X = X))
    return(ifelse(is(f, 'try-error'), NA, f))
  })
  return(dervs)
}