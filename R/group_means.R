#' Calculate outcome means per group per treatment level
#' 
#' @param Y name (character vector) of observed outcomes in data
#' @param A treatment levels 
#' @param G 
#' @param a 
#' @return data dataset to use
#' @export
#'

group_means <- function(Y, A, G, a, data){
  
  N <- length(unique(data[, G]))
  
  vals <- by(data, data[, G], function(x){
    n <- length(x[,Y])
    
    if(is.na(a)){
      sum(x[,Y])/n
    } else {
      sum(x[,Y]*(x[,A] == a)*1)/n
    }
  })
  
  out <- matrix(unlist(vals), nrow = N)
  
  return(out)
}
