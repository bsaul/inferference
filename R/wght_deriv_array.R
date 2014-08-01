#' Create an array of group weight derivative 
#' 
#' - creates a length(group) X length(theta) X length(alphas) array 
#'  
#' @param alphas coverage levels in (0, 1), possibly (probably) vector valued
#' @param data data frame
#' @param groups quoted string for name of variable in data containing group membership
#' @param predictors character vector of names of predictor variables in data
#' @param A character vector of name of treatment variable in data
#' @param theta p + 1 vector of fixed effects plus the random effect variance. The variance estimate must be last.
#' @param type type of weight to compute. See \code{\link{wght_calc}}
#' @param ... additional arguments passed to f.ab
#' @return a length(unique(group)) X length(alphas) matrix of group weights 
#' @export

wght_deriv_array <- function(f.ab, 
                             alphas, 
                             data, 
                             groups, 
                             predictors, 
                             treatment, 
                             theta, 
                             type, 
                             ...){
  
  # Make sure alphas are sorted
  alphas <- sort(alphas)
  
  G <- factor(data[, groups])
  X <- cbind(1, data[, predictors])
  A <- data[, treatment]
  p <- ncol(X)
  N <- length(unique(G))
  k <- length(alphas) 
  
  if(length(theta) != (p + 1)){
    stop("The length of theta is not equal to the number of predictors + 2 ")
  }
  
  w.list <- lapply(alphas, function(alpha){
    w <- by(cbind(X, A), INDICES = G, simplify = TRUE, 
            FUN = function(x) {
              x <- as.matrix(x) # PrAX expects a matrix
              wght_deriv_calc(f.ab = f.ab, type = type, alpha = alpha, 
                              A = x[, p+1], X = x[, 1:p], 
                              theta = theta, ...)})
    w2 <- matrix(unlist(w), ncol = p+1, byrow = TRUE,
                 dimnames= list(1:N, names(theta)))
    return(w2)}) 
  
  out <- array(unlist(w.list), dim = c(N, p+1, k),
               dimnames = list(1:N, names(theta), alphas))
  
  return(out)
}