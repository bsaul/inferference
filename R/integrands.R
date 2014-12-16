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
                             fixed.effects,
                             random.effect = NULL,
                             x = NULL, 
                             pos = NULL, 
                             allocation = NULL, 
                             randomization = 1, 
                             integrate.allocation = FALSE)
{
  ## Warnings ##
  if(length(fixed.effects) != ncol(X)){
    stop('The number of fixed effect parameters is not equal to the number \n
         of columns in the covariate matrix')
  }
  
  if(length(A) != nrow(X)){
    stop('Length of treatment vector is not equal to number of observations')
  }
  
  ## For taking derivative w.r.t. a parameter ##
  params <- c(fixed.effects, random.effect)
  if(!is.null(pos)){
    params[pos] <- x
  }
  
  ## X needs to be a matrix ##
  if(!is.matrix(X)){
    X <- as.matrix(X)
  }
  
  if(is.null(random.effect) || random.effect <= 0){
    pr.b <- randomization * (plogis(X %*% fixed.effects))
  } else {
    pr.b <- randomization * (plogis(drop(outer(X %*% params[1:length(fixed.effects)], b, '+'))))
  }
  
  if(integrate.allocation == FALSE){
    hh <- dbinom(A, 1, pr.b)
  } else {
    hh <- (pr.b/allocation)^A * ((1-pr.b)/(1 - allocation))^(1-A)
  }
  
  if(is.null(random.effect) || random.effect <= 0){
    # in this way dnorm integrates to one when integrating from -Inf to Inf
    out <- prod(hh) * dnorm(b, mean=0, sd = 1) 
  } else {
    hha <- apply(hh, 2, prod)
    out <- hha * dnorm(b, mean=0, sd = params[length(fixed.effects) + 1])
  }
  
  return(out)
}

