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