#' Compute the derivative(s) of a weight 
#' 
#' Takes the derivative of the \code{\link{wght_calc}} function with respect to each 
#' parameter in \code{params}.
#' 
#' @param integrand the function to passed to the argument 'f' of \code{\link{integrate}},
#' which is part of \code{\link{wght_calc}}.
#' @param params parameters with which to take derivatives with respect to
#' @param allocation the allocation ratio for which to compute the weights
#' @param hide.errors print \code{grad} error messages. Defaults to TRUE.
#' @param ... additional arguments passed to integrand
#' @return vector of derivatives with respect to element of params
#' @export

wght_deriv_calc <- function(integrand,
                            params, 
                            allocation,
                            hide.errors = TRUE,
                            ...)
{  
  integrand <- match.fun(integrand)
  dots <- list(...)
  
  fargs <- append(get_args(integrand, dots),
                  get_args(grad, dots))
  
  args <- append(fargs,
                 list(func   = wght_calc, 
                      integrand   = integrand, 
                      allocation  = allocation,
                      params  = params))
  
  dervs <- sapply(1:length(params), function(i){
    args$x <- params[i]; args$pos <- i
    f <- try(do.call('grad', args = args), silent = hide.errors)
    return(ifelse(is(f, 'try-error'), NA, f))
  })
  return(dervs)
}
