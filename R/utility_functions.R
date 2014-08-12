#' Get Arguments
#' 
#' Extracts the names of the arguments from a function, and creates a list 
#' of those arguments where they exist in ... . 
#' 
#' @param FUN function for which to find arguments
#' @param ... any arguments. Those necessary for FUN must be named as appropriate for FUN
#' @return list of arguments for FUN
#' @export
#' @examples
#' myargs <- get_args(lm, formula = Sepal.Length ~ Sepal.Width, data = iris )
#' summary(do.call('lm', myargs))

get_args <- function(FUN, args_list = NULL, ...){
  dots <- append(args_list, list(...))
  arg_names <- names(formals(FUN))
  
  args <- dots[arg_names]
  args[sapply(args, is.null)] <- NULL
  
  return(args)
}

#' Calculate outcome means per group per treatment level
#' 
#' @param Y name (character vector) of observed outcomes in data
#' @param A name of treatment variable
#' @param G name of grouping variable in data
#' @param a value of treatment level, defaults to NA.
#' @return data dataset to use
#' @export
#'

group_means <- function(Y, A, G, a = NA, data){
  
  N <- length(unique(data[ , G]))
  
  vals <- by(data, data[ , G], function(x){
    n <- length(x[ , Y])
    
    if(is.na(a)){
      sum(x[ , Y])/n
    } else {
      sum(x[ , Y] * (x[ , A] == a) * 1)/n
    }
  })
  
  out <- matrix(unlist(vals), nrow = N)
  
  return(out)
}