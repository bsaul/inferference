#' Calculate Pr(A|X)
#' 
#' \code{grod * dnorm(sd = sqrt(variance))}
#' 
#' @param b random effect necessary for \code{\link{integrate}}
#' @param X n by length(theta) - 1 matrix of covariates. Make sure the order of columns in X corresponds to theta
#' @param theta p + 1 vector of fixed effects plus the random effect variance. The variance estimate must be last.
#' @param A vector of observed treatments (0,1)
#' @param alpha The allocation strategy. This can be empty if type = 'c'.
#' @param r Randomization probability defaults to 2/3
#' @param type Either 'b' for including alpha in the product or 'c' does not include alpha.
#' @param x The argument passed to be \code{\link{grad}}. Only used if \code{pos} is not NA.
#' @param pos The position of theta for which to take the derivative. Defaults to NA.
#' 
#' @return value of the integrand
#' @export
#' 
PrAX_integrand <- function(b, x, pos = NA, X, A, theta, alpha, r = 2/3, type){
  n <- length(A)

  if(!is.na(pos)){
    theta[pos] <- x
  }

  hh <- 1
  for(jj in 1:n){
    pr.b <- r*plogis(X[jj, ] %*% theta[-length(theta)] + b)
    if(type == 'c'){
      hh <- hh * (pr.b)^A[jj] * (1-pr.b)^(1-A[jj])
    }
    else if(type == 'b'){
      hh <- hh * (pr.b/alpha)^A[jj] * ((1-pr.b)/(1 - alpha))^(1-A[jj])
    }
  }
  
  ans <- hh * dnorm(b, mean=0, sqrt(theta[length(theta)]))
  return(ans) 
}
