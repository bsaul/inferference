#' Create an array of group weight derivatives 
#' 
#' Uses \code{\link{wght_deriv_calc}} to compute the weight derivatives for each 
#' group per coverage level
#'  
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

wght_deriv_array <- function(integrand, 
                             allocations, 
                             X, A, G,
                             fixed.effects,
                             random.effect,
                             ...)
{
  ## Gather necessary bits ##
  integrand <- match.fun(integrand)
  XX <- cbind(X, A)
  p  <- length(fixed.effects) # number of predictors
  aa <- sort(allocations) # Make sure alphas are sorted
  gg <- sort(unique(G))
  k  <- length(allocations) 
  N  <- length(unique(G))

  ## Warnings ##

  
  ## Compute weight (derivative) for each group, parameter, and alpha level ##
  w.list <- lapply(aa, function(allocation){
    w <- by(XX, INDICES = G, simplify = TRUE, 
            FUN = function(x) {
              wght_deriv_calc(integrand = integrand, 
                              allocation = allocation, 
                              A = x[, p+1], 
                              X = x[, 1:p], 
                              fixed.effects = fixed.effects, 
                              random.effect = random.effect,
                              ...)})
    w2 <- matrix(unlist(w, use.names = FALSE), ncol = p+1, byrow = TRUE)
    return(w2)}) 
  
  ## Reshape list into array ##
  out <- array(unlist(w.list, use.names = FALSE), 
               dim = c(N, p+1, k),
               dimnames = list(gg, names(c(fixed.effects, random.effect)), aa))
  
  return(out)
}