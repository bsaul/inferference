#-----------------------------------------------------------------------------#
#' Simulate a dataset of N groups
#' 
#' Randomly sllocates n units to N groups and randomly samples a group distance
#' according to the chosen distribution
#'
#' @param N the number of groups
#' @param n the number of subjects
#' @param distrance.distr a character string defining the distribution from
#' which group level distances are sampled. Defaults to 'rlnorm'.
#' @param distr.mean the mean of the chosen distribution. Defaults to 0 as in 
#' Perez-Heydrich 2014.
#' @param distr.sd the standard deviation of the chosen distribution. Defaults
#' to 0.75.
#' @return a data frame with the size and distance for each group
#' @export
#-----------------------------------------------------------------------------#
sim_groups <- function(N, n, 
                       distance.distr = 'rlnorm', 
                       distr.mean = 0, 
                       distr.sd = 0.75){

  group <- sample(1:N, n, replace = T)
  grp.sizes <- as.data.frame(table(group))
   
  distance_group <- do.call(distance.distr, 
                            args = list(N, mean = distr.mean, sd = distr.sd))

  groups <- data.frame(group = as.character(grp.sizes[, 1]), 
                       group_size = grp.sizes[, 2],
                       distance_group,
                       stringsAsFactors = F)
  return(list(assignments = group, groups = groups))
}

#-----------------------------------------------------------------------------#
#' Simulate a dataset of n individuals into N groups
#' 
#' This is currently not very useful except to simulate a base dataset like in 
#' Step 0 of Perez et al's simulation
#'
#' @param N the number of groups
#' @param n the number of subjects
#' @param X1.args In Perez et al, X1 represents an individual's age in 
#' decades generated from an exponential distribution with mean 20
#' @param distr.mean the mean of the chosen distribution. Defaults to 0 as in 
#' Perez et al 2014.
#' @param distr.sd the standard deviation of the chosen distribution. Defaults
#' to 0.75.
#' @return a data frame with the size and distance for each group
#' @export
#-----------------------------------------------------------------------------#
sim_units <- function(N, n,
                      X1.distr = 'rexp',
                      X1.args  = list(rate = 1/20),
                      dist.distr = 'rnorm',
                      dist.args  = list(mean=0, sd=0.05),
                      grp.dist.distr = 'rlnorm', 
                      grp.dist.mean = 0, 
                      grp.dist.sd = 0.75)
{
  X1_ <- do.call(X1.distr, args = append(list(n = n), X1.args))
  
  ind_dist1 <- do.call(dist.distr, args = append(list(n = n), dist.args))
  ind_dist  <- ifelse(ind_dist1 < 0, 0, ind_dist1)

  groups.list <- sim_groups(N, n, 
                       distance.distr = grp.dist.distr, 
                       distr.mean = grp.dist.mean, 
                       distr.sd = grp.dist.sd)

  data <- merge(data.frame(id = 1:n, 
                           group = groups.list$assignments,
                           X1 = ifelse(X1_ > 100, 100, X1_/10)), 
                groups.list$groups, by="group")
  data$X2 <- ind_dist + data$distance_group
  
  return(data[, c('id', 'group', 'X1', 'X2')])
}

