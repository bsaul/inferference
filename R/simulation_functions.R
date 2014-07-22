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
  
  return(data[, c('id', 'group', 'X1', 'X2', 'group_size')])
}

#-----------------------------------------------------------------------------#
#' TBD
#'
#' @export
#-----------------------------------------------------------------------------#

po_y_i <- function(groupdt, theta){
  n <- groupdt$group_size[1]
  
  out <- sapply(0:1, function(a){
    t(sapply(0:(n-1), function(kk) {
    alpha <- (a + kk)/n
    X <- cbind(1, a, alpha, groupdt$X1, groupdt$X2, a*alpha)
    p  <- plogis(X %*% theta)
    sapply(p, function(pp) rbinom(1, 1, pp))
    }))
  }, simplify = 'array')
  
  dimnames(out) <- list(id = groupdt$id, k = 0:(n-1), A = 0:1)

  return(melt(out, value.name = 'Y'))
}

#-----------------------------------------------------------------------------#
#' TBD
#' 
#'  @export 
#-----------------------------------------------------------------------------#

sim_po <- function(base_dt, theta){
  po <- rbind.fill(lapply(split(base_dt, base_dt$group), po_y_i, theta))
  return(po)
}

#-----------------------------------------------------------------------------#
#' TBD
#' @export
#-----------------------------------------------------------------------------#

sim_data <- function(base_dt, potential.outcomes){

  N <- length(unique(base_dt$group))
  n <- nrow(base_dt)
  
  ## hood random effect ##
  # neighborhood level random effect from N(0, variance = 1.0859)
  group_effect <- data.frame(group = 1:N, 
                             b = rnorm(N, mean = 0, sd = sqrt(1.0859)))
  
  dt <- merge(base_dt, group_effect, by='group')
  dt <- within(dt, {
                    xsum = 0.2727 - 0.0387*X1 + 0.2179*X2 + b
                    B = rbinom(n, 1, p = exp(xsum)/(1 + exp(xsum)))
                    # Assign Treatments
                    A = ifelse(B == 0, 0, rbinom(n, 1, p=2/3))})
  # Calculate k = # treated in hood i excluding subject j 
  sim_dt <- ddply(dt, .(group), mutate, k = sum(A) - A)

  sim_dt <- merge(sim_dt, potential.outcomes, by=c('id', 'k', 'A'), all.x = T)
  #return(sim_dt[c('id', 'X1', "X2", "A", 'B', 'hood', 'alpha', 'k')])
  return(sim_dt[c('id', 'Y', 'X1', "X2", "A", 'B', 'group')])
}

#-----------------------------------------------------------------------------#
#' TBD
#' @export 
#-----------------------------------------------------------------------------#

replicate_sims <- function(base_dt, potential.outcomes, nsims){
  out <- list(nsims)
  for(sim in 1:nsims){
    out[[sim]] <- sim_data(base_dt, potential.outcomes)
  }
  return(out)
}

