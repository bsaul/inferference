#' Create a matrix of group weights 
#' 
#' - creates a length(group) X l matrix of group weights 
#'  
#' @param alphas coverage levels in (0, 1), possibly (probably) vector valued
#' @param data data frame
#' @param groups quoted string for name of variable in data containing group membership
#' @param predictors character vector of names of predictor variables in data
#' @param A character vector of name of treatment variable in data
#' @param theta p + 1 vector of fixed effects plus the random effect variance. The variance estimate must be last.
#' @param type type of weight to compute. See \code{\link{wght_calc}}
#' @return a length(unique(group)) X length(alphas) matrix of group weights 
#' @export

wght_matrix <- function(alphas, data, groups, predictors, A, theta, type){
  
  # Make sure alphas are sorted
  alphas <- sort(alphas)
  
  G <- data[, groups]
  X <- cbind(1, data[, predictors])
  p <- ncol(X)
  A <- data[, A]
  
  w.list <- lapply(alphas, function(alpha){
    w <- by(cbind(X, A), INDICES = G, simplify = TRUE, 
            FUN = function(x) {
              x <- as.matrix(x) # PrAX expects a matrix
              wght_calc(type = type, alpha = alpha, 
                        A = x[, p+1], X = x[, 1:p], theta = theta)})
    as.numeric(w)
  }) 
  
  w.matrix <- matrix(unlist(w.list), ncol = length(alphas), byrow = FALSE,
                     dimnames = list(sort(unique(G)), alphas))
  
  return(w.matrix)
}
