#' Log Likelihood integrand
#' 
#' a function
#' 
#' @param b
#' @param B
#' @param X
#' @param x
#' @param pos
#' @return TBD
#' @export

lik_integrand <- function(b, theta, B, X, x, pos){  
  n <- length(B)
  theta[pos] <- x
  theta.fix <- theta[1:ncol(X)]
  
  p <- plogis(drop(outer(X %*% theta.fix, b, '+'))) 
  ans <- dbinom(B, 1, p)
  
  ans <- apply(ans, 2, prod)
  return(ans * dnorm(b, 0, theta[length(theta)]))
}

#' Log Likelihood 
#' 
#' TBD
#' 
#' @param x
#' @param pos
#' @param theta
#' @param B
#' @param X
#' @return value of log likelihood
#' @export

ll <- function(x, pos, theta, B, X){
  y <- log(integrate(lik_integrand, lower=-Inf,upper=Inf, theta = theta, 
                     B = B, X = X, x = x, pos = pos)$value)
  return(y)
}

#' Compute scores
#' 
#' TBD
#' 
#' @param theta
#' @param B
#' @param X 
#' @return Score
#' @export

cmp_scores_ll <- function(theta, B, X, hide.errors = TRUE){
  scores <- sapply(1:length(theta), function(i){
    f <- try(grad(ll, x = theta[i], pos = i, theta = theta, B = B, X = X),
             silent = hide.errors)
    return(ifelse(is(f, 'try-error'), NA, f))
  })
  return(scores)
}