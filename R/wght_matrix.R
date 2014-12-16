#' Create a matrix of group weights 
#' 
#' Creates a number of groups by number of allocation schemes matrix of group weights.
#' Allocation schemes are selected by the user. 
#' 
#' Groups should be numbered 1, ..., N
#' 
#' \code{predictors} does not need to have a value corresponding to reference parameter. 
#' The function adds a column of 1s with \code{cbind()} as a first step.
#'  
#' @param integrand the function used in the weight calculation. Defaults to 
#' \code{\link{logit_integrand}}
#' @param allocations coverage levels in (0, 1), possibly (probably) vector valued
#' @param data data frame
#' @param groups quoted string for name of variable in data containing group membership
#' @param predictors character vector of names of predictor variables in data
#' @param A character vector of name of treatment variable in data
#' @param params p + 1 vector of fixed effects plus the random effect variance. 
#' The variance estimate must be last.
#' @param ... additional arguments passed to \code{integrand}
#' @return a length(unique(group)) X length(alphas) matrix of group weights 
#' @export

wght_matrix <- function(integrand = logit_integrand, 
                        allocations, 
                        X, A, G,
                        fixed.effects,
                        random.effect = NULL,
                        ...)
{
  ## Gather necessary bits ##
  XX <- cbind(X, A)
  p  <- length(fixed.effects)
  aa <- sort(allocations) # Make sure alphas are sorted
  gg <- sort(unique(G))
  
  ## Compute weight for each group and alpha level ##
  w.list <- lapply(aa, function(allocation){
    w <- by(cbind(X, A), INDICES = G, simplify = FALSE, 
            FUN = function(x) {
              wght_calc(integrand = integrand, 
                        allocation = allocation, 
                        A = x[, p+1], X = x[, 1:p], 
                        fixed.effects = fixed.effects, 
                        random.effect = random.effect, ...)})
    as.numeric(w)
  }) 
  
  ## Reshape list into matrix ##
  w.matrix <- matrix(unlist(w.list, use.names = FALSE), 
                     ncol = length(allocations), 
                     byrow = FALSE,
                     dimnames = list(gg, aa))
  
  return(w.matrix)
}
