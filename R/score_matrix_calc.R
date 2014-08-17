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
#' @param theta p + 1 vector of fixed effects plus the random effect variance. 
#' The variance estimate must be last.
#' @return N X length(theta) matrix of scores
#' @export
#-----------------------------------------------------------------------------#

score_matrix_calc <- function(integrand = logit_integrand,
                              predictors, 
                              treatment, 
                              groups, 
                              theta, 
                              data,
                              set.NA.to.0 = TRUE, 
                              ...){
  ## Warnings ##
  if(length(theta) != (length(predictors) + 2)){
    stop("The length of theta is not equal to the number of predictors + 2 ")
  }
  
  ## Necessary bits ##
  integrand <- match.fun(integrand)
  dots <- list(...)
  G  <- data[, groups]
  XX <- cbind(1, data[, c(predictors, treatment)])
  p  <- ncol(XX) - 1
  gg <- unique(G)
  
  ## Compute score for each group and parameter ##
  s.list <- by(XX, INDICES = G, simplify = TRUE, 
               FUN = function(xx) {
               xx <- as.matrix(xx)
               args <- append(list(integrand = integrand, theta = theta,
                            A = xx[ , (p + 1)],
                            X = xx[ , 1:p]), 
                            get_args(integrand, dots))
               return(do.call('score_calc', args = args))})
  
  ## Reshape list into matrix ##
  out <- matrix(unlist(s.list, use.names = FALSE), 
                ncol = p + 1, 
                byrow = TRUE,
                dimnames = list(gg, names(theta)))
  
  ## replace any Bscores with 0 ##
  if(set.NA.to.0 == TRUE) {
    out[is.na(out)] <- 0
  }
  
  return(out)
}