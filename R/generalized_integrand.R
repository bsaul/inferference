## Generalized integrand for the group-level propensity score ##

generalized_integrand <- function(b, X, A, 
                            fixed.effects,
                            random.effects = NULL,
                            x = NULL, 
                            pos = NULL, 
                            allocation = NULL, 
                            randomization = 1, 
                            integrate.allocation = FALSE,
                            reff_dist = list(
                              dens = stats::dnorm,
                              mean=0),
                            inv_link_fun)
{
  p  <- length(fixed.effects)
  re <- random.effects[1]
  
  ## In the case of an intercept-only model, X needs to be converted to matrix
  # for the warning to work
  if(!is.matrix(X)){
    X <- as.matrix(X)
  }
  
  ## Warnings ##
  if(p != ncol(X)){
    stop('The number of fixed effect parameters is not equal to the number \n
         of columns in the covariate matrix')
  }
  
  if(length(A) != nrow(X)){
    stop('Length of treatment vector is not equal to number of observations')
  }
  
  if(!all(is.element(c("dens", "mean"), names(reff_dist)))){
    stop("reff_dist needs to provide a density and a mean parameter for generalized_integrand")
  }
  
  if(is.null(inv_link_fun)){
    stop("Must specify inv_link_fun in generalized_integrand")
  }
  
  ## For taking derivative w.r.t. a parameter ##
  params <- c(fixed.effects, re)
  if(!is.null(pos)){
    params[pos] <- x
  }
  
  ## Calculations ## 
  if(is.null(re) || re <= 0){
    pr.b <- randomization * (inv_link_fun(X %*% params[1:p]))
  } else {
    pr.b <- randomization * (inv_link_fun(drop(outer(X %*% params[1:p], b, '+'))))
  }
  if(integrate.allocation == FALSE){
    hh <- dbinom(A, 1, pr.b)
  } else {
    hh <- (pr.b/allocation)^A * ((1-pr.b)/(1 - allocation))^(1-A)
  }
  if(is.null(re) || re <= 0){
    # in this way reff_dens integrates to one when integrating from -Inf to Inf
    out <- prod(hh) * reff_dist$dens(b, reff_dist$mean=0, sd = 1) 
  } else {
    hha <- apply(hh, 2, prod)
    out <- hha * reff_dist$dens(b, reff_dist$mean=0, sd = params[p + 1])
  }
  
  return(out)
  }

