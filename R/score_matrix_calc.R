#-----------------------------------------------------------------------------#
#' Calculate matrix of log Likelihood derivatives  
#'
#' @param integrand function passed to \code{\link{log_likelihood}}. Defaults to
#' \code{\link{logit_integrand}}
#' @param predictors character vector of predictors from model
#' @param treatment character string naming treatment variable in data. As in Perez
#' 2014, this does not necessarily need to be the same 'treatment' used to compute 
#' the weights as in \code{\link{wght_matrix_calc}}. The \code{propensityB}
#' argument in \code{\link{run_interference}} may be used to pass a different 
#' indicator variable other than the actual treatment or exposure. In the case 
#' of the Perez paper.
#' @param data data frame
#' @param groups character string of the group variable in data 
#' @param params p + 1 vector of fixed effects plus the random effect variance. 
#' The variance estimate must be last.
#' @return N X length(params) matrix of scores
#' @export
#-----------------------------------------------------------------------------#

score_matrix_calc <- function(integrand = logit_integrand,
                              predictors, 
                              treatment, 
                              groups, 
                              params, 
                              data,
                              ...){
  ## Warnings ##
  if(length(params) != (length(predictors) + 2)){
    stop("The length of params is not equal to the number of predictors + 2 ")
  }
  
  ## Necessary bits ##
  integrand <- match.fun(integrand)
  dots <- list(...)
  G  <- data[, groups]
  XX <- cbind(1, data[, c(predictors, treatment)])
  pp <- ncol(XX) - 1
  gg <- sort(unique(G))
  
  ## Compute score for each group and parameter ##
  s.list <- by(XX, INDICES = G, simplify = TRUE, 
               FUN = function(xx) {
               xx <- as.matrix(xx)
               args <- append(list(integrand = integrand, params = params,
                            A = xx[ , (pp + 1)],
                            X = xx[ , 1:pp]), 
                            get_args(integrand, dots))
               return(do.call(score_calc, args = args))})
  
  ## Reshape list into matrix ##
  out <- matrix(unlist(s.list, use.names = FALSE), 
                ncol = pp + 1, 
                byrow = TRUE,
                dimnames = list(gg, names(params)))
  
  return(out)
}