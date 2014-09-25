#' Integrand of Pr(A|X)
#' 
#' Computes the following function:
#' \deqn{\prod_{j=1}^{n} (r h_{j}(b))^{A_j}  (1 - r h_{j}(b))^{1 - A_j} 
#' f_b(b; \theta_b)}{ prod(r * plogis(X * fixef + b)^A * 
#' (1 - r * plogis(X * fixef+ b))^(1 - A)) * 
#' dnorm(sd = sqrt(ranef))} 
#' where \eqn{r} is the randomization scheme. \eqn{X} is the covariate(s) vectors. 
#' \eqn{fixef} is the vector of fixed effects. \eqn{b} is the random (group-level) effect.
#' \eqn{ranef} is the random effect variance. Used by \code{\link{wght_calc}} and
#' \code{\link{score_calc}}.
#' 
#' @param b vector argument of values necessary for \code{\link{integrate}}
#' @param x Used by \code{\line{grad}} for taking the derivative with respect an element of
#' params. Only used if \code{pos} is not NULL.
#' @param pos The position of theta for which to take the derivative. Defaults to NULL.
#' @param X n by length(params) - 1 matrix of covariates. Make sure the order of columns in X corresponds to params.
#' @param params p + 1 vector of fixed effects plus the random effect variance. The variance estimate must be the last element.
#' @param A vector of observed treatments (0,1)
#' @param allocation The allocation strategy. Required if include.alpha == TRUE. 
#' Defaults to NA.
#' @param r Randomization probability. Defaults to 1.
#' @param include.allocation Either TRUE for including allocation in the product or FALSE 
#' does not include allocation. See \code{\link{wght_calc}} for more information.
#' 
#' @return value of the integrand
#' @export
#' 

logit_integrand <- function(b, 
                            X, 
                            A, 
                            params, 
                            x = NULL, 
                            pos = NULL, 
                            allocation = NULL, 
                            r = 1, 
                            include.allocation = FALSE){
  
  ## Warnings ##
  if(length(params) - 1 != ncol(X)){
    stop('The number of fixed effect parameters is not equal to the number \n
         of columns in the covariate matrix')
  }
  
  if(!is.null(pos)){
    params[pos] <- x
  }
  
  # X needs to be a matrix
  if(!is.matrix(X)){
    X <- as.matrix(X)
  }
  
  theta.fix <- params[1:ncol(X)]
  theta.ran <- params[length(params)]
  
  pr.b <- r * (plogis(drop(outer(X %*% theta.fix, b, '+'))))  
  
  if(include.allocation == FALSE){
    hh <- dbinom(A, 1, pr.b)
  } else {
    hh <- (pr.b/allocation)^A * ((1-pr.b)/(1 - allocation))^(1-A)
  }
  
  hha <- apply(hh, 2, prod)
  
  return(hha * dnorm(b, mean=0, theta.ran))
}
