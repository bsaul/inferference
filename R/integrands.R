#' Integrand of Pr(A|X)
#' 
#' Computes the following function:
#' \deqn{\prod_{j=1}^{n} (r h_{j}(b))^{A_j}  (1 - r h_{j}(b))^{1 - A_j} 
#' f_b(b; \theta_b)}{ prod(r * plogis(X * fixef + b)^A * 
#' (1 - r * plogis(X * fixef+ b))^(1 - A)) * 
#' dnorm(sd = sqrt(ranef))} 
#' where \eqn{r} is the randomization scheme. \eqn{X} is the covariate(s) vectors. 
#' \eqn{fixef} is the vector of fixed effects. \eqn{b} is the random (group-level) effect.
#' \eqn{ranef} is the random effect variance.
#' 
#' @param b vector argument of values necessary for \code{\link{integrate}}
#' @param x The argument passed to be \code{\link{grad}}. Only used if \code{pos} is not NA.
#' @param pos The position of theta for which to take the derivative. Defaults to NA.
#' @param X n by length(theta) - 1 matrix of covariates. Make sure the order of columns in X corresponds to theta
#' @param theta p + 1 vector of fixed effects plus the random effect variance. The variance estimate must be the last element.
#' @param A vector of observed treatments (0,1)
#' @param alpha The allocation strategy. Required if include.alpha == TRUE. 
#' Defaults to NA.
#' @param r Randomization probability. Defaults to 2/3.
#' @param include.alpha Either TRUE for including alpha in the product or FALSE 
#' does not include alpha.
#' 
#' @return value of the integrand
#' @export
#' 

logit_integrand <- function(b, 
                           x, 
                           pos = NA, 
                           X, 
                           A, 
                           theta, 
                           alpha = NA, 
                           r = 1, 
                           include.alpha = FALSE){
  
  ## Warnings ##
  if(length(theta) - 1 != ncol(X)){
    stop('The number of fixed effect parameters is not equal to the number \n
         of columns in the covariate matrix')
  }
  
  if(!is.na(pos)){
    theta[pos] <- x
  }
  
  if(!is.matrix(X)){
    X <- as.matrix(X)
  }
  
  theta.fix <- theta[1:ncol(X)]
  
  pr.b <- r * (plogis(drop(outer(X %*% theta.fix, b, '+'))))  
  
  if(include.alpha == FALSE){
    #hh <- pr.b^A * (1 - pr.b)^(1 - A)
    hh <- dbinom(A, 1, pr.b)
  } else {
    hh <- (pr.b/alpha)^A * ((1-pr.b)/(1 - alpha))^(1-A)
  }
  
  hh_ <- apply(hh, 2, prod)
  
  return(hh_ * dnorm(b, mean=0, sqrt(theta[length(theta)])))
}
