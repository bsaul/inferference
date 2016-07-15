#-----------------------------------------------------------------------------#
#' Calculate matrix of log Likelihood derivatives
#' 
#' @param integrand function passed to \code{\link{log_likelihood}}. Defaults to
#' \code{\link{logit_integrand}}
#' @param X covariate matrix
#' @param A vector of treatment assignments
#' @param G vector of group assignments
#' @param parameters vector of parameters passed to \code{integrand}
#' @param runSilent If FALSE, prints errors to console. Defaults to TRUE.
#' @param ... additional arguments passed to \code{integrand} or \code{\link{grad}}.
#' For example, one can change the \code{method} argument in \code{grad}.
#' @return N X length(params) matrix of scores
#' @export
#-----------------------------------------------------------------------------#

score_matrix <- function(integrand,
                         X, A, G, 
                         parameters,
                         runSilent = TRUE, 
                         ...)
{
  ## Warnings ##
  # if(length(fixed.effects) != ncol(X)){
  #   stop("The length of params is not equal to the number of predictors")
  # }
  
  ## Necessary bits ##
  integrand <- match.fun(integrand)
  dots <- list(...)
  XX <- cbind(X, A)
  pp <- ncol(X)
  gg <- sort(unique(G))
  
  ## Compute score for each group and parameter ##
  int.args <- append(get_args(integrand, dots),
                     get_args(stats::integrate, dots))
  fargs <- append(int.args, get_args(numDeriv::grad, dots))
  
  if(!runSilent) print("Calculating matrix of scores...")
  
  s.list <- by(XX, INDICES = G, simplify = TRUE, 
               FUN = function(xx) {
                 args <- append(fargs, 
                              list(integrand  = integrand, 
                                   parameters = parameters,
                                   A = xx[ , (pp + 1)],
                                   X = xx[ , 1:pp]))
                 do.call(score_calc, args = args)
               })

  ## Reshape list into matrix ##
  out <- matrix(unlist(s.list, use.names = FALSE), 
                ncol = length(parameters), 
                byrow = TRUE,
                dimnames = list(gg, names(parameters)))
  
  out 
}