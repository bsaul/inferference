#' Create a matrix of group weights 
#' 
#' Creates a number of groups by number of allocation scheme matrix of group weights.
#' Allocation schemes are selected by the user. 
#' 
#' Groups should be numbered 1, ..., N
#' 
#' \code{predictors} does not need to have a value corresponding to reference parameter. 
#' The function adds a column of 1s with \code{cbind()} as a first step.
#'  
#' @param f.ab the integrand used in the weight calculation. Defaults to 
#' \code{\link{logit_integrand}}
#' @param alphas coverage levels in (0, 1), possibly (probably) vector valued
#' @param data data frame
#' @param groups quoted string for name of variable in data containing group membership
#' @param predictors character vector of names of predictor variables in data
#' @param A character vector of name of treatment variable in data
#' @param params p + 1 vector of fixed effects plus the random effect variance. 
#' The variance estimate must be last.
#' @param type type of weight to compute. See \code{\link{wght_calc}}
#' @param ... additional arguments passed to \code{f.ab}
#' @return a length(unique(group)) X length(alphas) matrix of group weights 
#' @export

wght_matrix <- function(f.ab = logit_integrand, 
                        alphas, 
                        data, 
                        groups, 
                        predictors, 
                        treatment, 
                        params, 
                        include.alpha, 
                        ...){
  ## Gather necessary bits ##
  G  <- data[, groups]
  X  <- cbind(1, data[, predictors])
  A  <- data[, treatment]
  p  <- ncol(X) # number of predictors
  aa <- sort(alphas) # Make sure alphas are sorted
  gg <- unique(G)
  
  ## Compute weight for each group and alpha level ##
  w.list <- lapply(aa, function(alpha){
    w <- by(cbind(X, A), INDICES = G, simplify = FALSE, 
            FUN = function(x) {
              x <- as.matrix(x) # PrAX expects a matrix
              wght_calc(f.ab = f.ab, include.alpha = include.alpha, 
                        alpha = alpha, 
                        A = x[, p+1], X = x[, 1:p], 
                        params = params, ...)})
    as.numeric(w)
  }) 
  
  ## Reshape list into matrix ##
  w.matrix <- matrix(unlist(w.list, use.names = FALSE), 
                     ncol = length(alphas), 
                     byrow = FALSE,
                     dimnames = list(gg, aa))
  
  return(w.matrix)
}
