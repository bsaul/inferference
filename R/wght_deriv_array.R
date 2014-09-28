#' Create an array of group weight derivatives 
#' 
#' Uses \code{\link{wght_deriv_calc}} to compute the weight derivatives for each 
#' group per coverage level
#' 
#' @param integrand the function used in the weight calculation. Defaults to 
#' \code{\link{logit_integrand}} 
#' @param allocations coverage levels in (0, 1), possibly (probably) vector valued
#' @param data data frame
#' @param groups quoted string for name of variable in data containing group membership
#' @param predictors character vector of names of predictor variables in data
#' @param treatment character vector of name of treatment variable in data
#' @param params p + 1 vector of fixed effects plus the random effect variance. 
#' The variance estimate must be last.
#' @param ... additional arguments passed to integrand
#' @return a length(unique(group)) X length(params) X length(alphas) array of 
#' group weight derivatives
#' @export

wght_deriv_array <- function(integrand = logit_integrand, 
                             allocations, 
                             data, 
                             groups, 
                             predictors, 
                             treatment, 
                             params, 
                             ...){
  ## Gather necessary bits ##
  integrand <- match.fun(integrand)
  G  <- data[, groups]
  X  <- cbind(1, data[, predictors])
  A  <- data[, treatment]
  p  <- ncol(X) # number of predictors
  aa <- sort(allocations) # Make sure alphas are sorted
  gg <- sort(unique(G))
  k  <- length(allocations) 
  N  <- length(unique(G))

  ## Warnings ##
  if(length(params) != (p + 1)){
    stop("The length of params is not equal to the number of predictors + 2 ")
  }
  
  ## Compute weight (derivative) for each group, parameter, and alpha level ##
  w.list <- lapply(aa, function(allocation){
    w <- by(cbind(X, A), INDICES = G, simplify = TRUE, 
            FUN = function(x) {
              x <- as.matrix(x) # PrAX expects a matrix
              wght_deriv_calc(integrand = integrand, 
                              allocation = allocation, 
                              A = x[, p+1], 
                              X = x[, 1:p], 
                              params = params, ...)})
    w2 <- matrix(unlist(w, use.names = FALSE), ncol = p+1, byrow = TRUE)
    return(w2)}) 
  
  ## Reshape list into array ##
  out <- array(unlist(w.list, use.names = FALSE), 
               dim = c(N, p+1, k),
               dimnames = list(gg, names(params), aa))
  
  return(out)
}