#-----------------------------------------------------------------------------#
#' Get arguments from a function
#' 
#' Extracts the names of the arguments from a function, and creates a list 
#' of those arguments where they exist in ... . 
#' 
#' @param FUN function for which to find arguments
#' @param args_list a list of arguments. Defaults to NULL.
#' @param ... any arguments. Those necessary for FUN must be named as appropriate for FUN
#' @return list of arguments for FUN
#' @export
#' @examples
#' myargs <- get_args(lm, formula = Sepal.Length ~ Sepal.Width, data = iris )
#' summary(do.call('lm', myargs))
#-----------------------------------------------------------------------------#

get_args <- function(FUN, args_list = NULL, ...){
  dots <- append(args_list, list(...))
  arg_names <- names(formals(match.fun(FUN)))
  
  args <- dots[arg_names]
  args[sapply(args, is.null)] <- NULL
  
  return(args)
}

#-----------------------------------------------------------------------------#
#' Calculate outcome means per group per treatment level
#' 
#' @param Y name (character vector) of observed outcomes in data
#' @param A name of treatment variable
#' @param G name of grouping variable in data
#' @param a value of treatment level, defaults to NA.
#' @param data data.frame with variables
#' @return data dataset to use
#-----------------------------------------------------------------------------#

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

#-----------------------------------------------------------------------------#
#' Create a dataset of arguments for effect estimates
#'
#' @param allocations vector of allocations 
#' @param treatments vector of treatments. defaults to \code{c(0 ,1)}
#' @return data.frame with arguments necessary for \code{\link{calc_effect}} to 
#' compute all outcome, direct, indirect, total, and overall effect estimates from
#' an object created from \code{\link{ipw_interference}} 
#' @export
#' @examples 
#' effect_grid(seq(0,1, by = .1), c(0,1))
#' 
#-----------------------------------------------------------------------------#

effect_grid <- function(allocations, treatments = c(0,1))
{
  marginal    <- c('TRUE', 'FALSE')
  
  # Outcomes
  g1 <- expand.grid(alpha1 = allocations, trt1 = treatments, 
                    alpha2 = NA, trt2 = NA,
                    marginal = marginal, 
                    effect_type = 'outcome', effect = 'outcome')
  
  # Direct Effects
  g2 <- expand.grid(alpha1 = allocations, trt1 = treatments, 
                    alpha2 = NA, trt2 = treatments,
                    marginal = FALSE, 
                    effect_type = 'contrast', effect = 'direct')
  g2$alpha2 <- g2$alpha1
  g2 <- g2[g2$trt1 != g2$trt2, ]
  
  # Indirect Effects
  g3 <- expand.grid(alpha1 = allocations, trt1 = treatments, 
                    alpha2 = allocations, trt2 = NA,
                    marginal = FALSE, 
                    effect_type = 'contrast', effect = 'indirect')
  g3 <- g3[g3$alpha1 != g3$alpha2, ]
  g3$trt2 <- g3$trt1
  
  # Total Effects
  g4 <- expand.grid(alpha1 = allocations, trt1 = treatments, 
                    alpha2 = allocations, trt2 = treatments,
                    marginal = FALSE, 
                    effect_type = 'contrast', effect = 'total')
  g4 <- g4[g4$alpha1 != g4$alpha2, ]
  g4 <- g4[g4$trt1 != g4$trt2, ]
  
  # Overall Effects
  g5 <- expand.grid(alpha1 = allocations, trt1 = NA, 
                    alpha2 = allocations, trt2 = NA,
                    marginal = TRUE, 
                    effect_type = 'contrast', effect = 'overall')
  g5 <- g5[g5$alpha1 != g5$alpha2, ]
  
  out <- rbind(g1, g2, g3, g4, g5)
  return(out)
}
