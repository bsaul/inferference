#' Create an array of group weight derivative 
#' 
#' - creates a length(group) X length(theta) X length(alphas) array 
#'  
#' @param alphas coverage levels in (0, 1), possibly (probably) vector valued
#' @param data data frame
#' @param groups quoted string for name of variable in data containing group membership
#' @param predictors character vector of names of predictor variables in data
#' @param A character vector of name of treatment variable in data
#' @param params p + 1 vector of fixed effects plus the random effect variance. The variance estimate must be last.
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
  k  <- length(alphas) 
  N  <- length(unique(G))

  ## Warnings ##
  if(length(params) != (p + 1)){
    stop("The length of params is not equal to the number of predictors + 2 ")
  }
  
  ## Compute weight (derivative) for each group, parameter, and alpha level ##
  w.list <- lapply(aa, function(alpha){
    w <- by(cbind(X, A), INDICES = G, simplify = TRUE, 
            FUN = function(x) {
              x <- as.matrix(x) # PrAX expects a matrix
              wght_deriv_calc(f.ab = f.ab, include.alpha = include.alpha,
                              alpha = alpha, 
                              A = x[, p+1], X = x[, 1:p], 
                              params = params, ...)})
    w2 <- matrix(unlist(w, use.names = FALSE), ncol = p+1, byrow = TRUE)
    return(w2)}) 
  
  ## Reshape list into array ##
  out <- array(unlist(w.list, use.names = FALSE), 
               dim = c(N, p+1, k),
               dimnames = list(gg, names(params), aa))
  
  return(out)
}