#' Compute the derivative(s) of a weight 
#' 
#' @param f.ab the function to passed to the argument 'f' of \code{\link{integrate}},
#' which is part of \code{\link{wght_calc}}.
#' @param params parameters with which to take derivatives with respect to
#' @param type see \code{\link{wght_calc}}. Defaults to 'b'.
#' @param method see \code{\link{grad}}.
#' @param method.args see \code{\link{grad}}
#' @param hide.errors print \code{grad} error messages. Defaults to TRUE.
#' @param ... additional arguments passed to f.ab
#' @return vector of derivatives with respect to element of params
#' @export

wght_deriv_calc <- function(f.ab,
                            params, 
                            alpha,
                            include.alpha,
                            method = "simple",
                            method.args = list(eps = 1e-4),
#                             method = "Richardson", 
#                             method.args = list(eps=1e-4, d=0.0001, 
#                                                 zero.tol=sqrt(.Machine$double.eps/7e-7), 
#                                                 r=4, v=2, show.details=FALSE),
                            hide.errors = TRUE,
                            ...)
{  
  f.ab <- match.fun(f.ab)
  args <- append(list(func   = wght_calc, 
                      f.ab   = f.ab, 
                      alpha  = alpha,
                      params  = params,
                      include.alpha   = include.alpha,
                      method = method, 
                      method.args = method.args), 
                 get_args(f.ab, ...)) # Get the necessary arguments for f.ab
  
#   if("params" %in% names(formals(f.ab))){
#     args$params <- params
#   }
  
#   print(names(args))
  dervs <- sapply(1:length(params), function(i){
    args$x <- params[i]; args$pos <- i
    f <- try(do.call('grad', args = args), silent = hide.errors)
    return(ifelse(is(f, 'try-error'), NA, f))
  })
  return(dervs)
}
