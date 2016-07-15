#-----------------------------------------------------------------------------#
#' Default integrand for the group-level propensity score
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
#' @param b vector argument of values necessary for \code{\link{integrate}}.
#' @param X n by length(fixed effects) matrix of covariates.
#' @param parameters vector of fixed effect (and random effect if applicable). 
#' Random effect should be last element in vector.
#' @param A vector of observed treatments (0,1)
#' @param allocation The allocation strategy. Required if include.allocations == TRUE. 
#' Defaults to NA.
#' @param randomization Randomization probability. Defaults to 1.
#' 
#' @return value of the integrand
#' @export
#' 
#-----------------------------------------------------------------------------#

logit_integrand <- function(b, X, A, 
                            parameters,
                            allocation = A, 
                            randomization = 1)
{
  ## In the case of an intercept-only model, X needs to be converted to matrix
  # for the warning to work
  if(!is.matrix(X)){
    X <- as.matrix(X)
  }
  
  theta <- parameters 
  p <- ncol(X)
  
  ## Warnings ##
  # if(p != ncol(X)){
  #   stop('The number of fixed effect parameters is not equal to the number \n
  #        of columns in the covariate matrix')
  # }
  
  if(length(A) != nrow(X)){
    stop('Length of treatment vector is not equal to number of observations in
         X matrix')
  }
  
  # Check whether to ignore random effect
  ignore_re <- (length(theta) == p || theta[p + 1] <= 0)

  ## Calculations ## 
  if(ignore_re){
    pr.b <- randomization * (stats::plogis(X %*% theta[1:p]))
  } else {
    pr.b <- randomization * (stats::plogis(drop(outer(X %*% theta[1:p], b, '+'))))
  }
  
  hh <- (pr.b/allocation)^A * ((1-pr.b)/(1 - allocation))^(1-A)
  
  if(ignore_re){
    # in this way dnorm integrates to one when integrating from -Inf to Inf
    out <- exp(sum(log(hh))) * stats::dnorm(b, mean=0, sd = 1) 
  } else {
    hha <- apply(hh, 2, function(x) exp(sum(log(x))))
    out <- hha * stats::dnorm(b, mean=0, sd = theta[p + 1])
  }
  
  return(out)
}

