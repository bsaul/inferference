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
                           X1 = ifelse(X1_ > 100, 100/10, X1_/10)), 
                groups.list$groups, by="group")
  data$X2 <- ind_dist + data$distance_group
  
  return(data[, c('id', 'group', 'X1', 'X2', 'group_size')])
}

#-----------------------------------------------------------------------------#
#' Individual Potential outcomes
#'
#' @export
#-----------------------------------------------------------------------------#

po_y_ij <- function(dt, a, k, theta){
  n <- nrow(dt)
  alpha <- (a + k)/n
  X <- cbind(1, a, alpha, dt$X1, dt$X2, a*alpha)
  p  <- plogis(X %*% theta)
  return(rbinom(n, 1, p))
}

#-----------------------------------------------------------------------------#
#' Group level potential outcomes
#' 
#'  @export 
#-----------------------------------------------------------------------------#

gen_po_ij <- function(dt, theta){
  n <- nrow(dt)
  ids <- dt$id
  po <- array(dim = c(n, n, 2), 
              dimnames = list(id = ids, k = 0:(n-1), trt = c(0,1)))
  
  for(k in 0:(n-1)){
    po[ , as.character(k), 1] <- po_y_ij(dt, 0, k, theta)
    po[ , as.character(k), 2] <- po_y_ij(dt, 1, k, theta)
  }
  
  return(po)
}

#-----------------------------------------------------------------------------#
#' Generate Potential Outcomes
#' @export
#-----------------------------------------------------------------------------#

make_po <- function(base, po.parameters){
  po <- lapply(split(base, base$group), gen_po_ij, theta = po.parameters)
  return(po)
}

#-----------------------------------------------------------------------------#
#' Individual level causal estimand Y(a, alpha)
#' 
#' Takes n X n - 1 x 2 array as input
#' @return n x 2 matrix with values of ybar_ij(a, alpha) 
#' @export
#-----------------------------------------------------------------------------#


ybar_ij_aalpha <- function(outcomes, alpha){
  
  n <- dim(outcomes)[1]
  
  # Keep the IDs along for the ride
  po_ybar <- array(dim = c(n, 2), dimnames = list(id = dimnames(outcomes)$id,
                                                  trt = dimnames(outcomes)$trt))
  
  #TODO: Speed up with apply() functions
  for(j in 1:n){
    y0 <- y1 <- 0
    for(k in 0:(n-1)){
      y0 <- y0 + outcomes[j, as.character(k), 1] * choose(n - 1, k) * alpha^k*(1-alpha)^(n - k -1)
      y1 <- y1 + outcomes[j, as.character(k), 2] * choose(n - 1, k) * alpha^k*(1-alpha)^(n - k -1)
    }
    po_ybar[j,1] <- y0
    po_ybar[j,2] <- y1
  }
  return(po_ybar)
}

#-----------------------------------------------------------------------------#
#' Individual level marginal causal estimand Y(alpha)
#' Takes the n x 2 matrix returned by ybar_ij_aalpha
#' 
#' @return n x 1 vector of ybar_ij(alpha) 
#' @export
#-----------------------------------------------------------------------------#
ybar_ij_alpha <- function(input, alpha){
  alpha*input[, 2] + (1 - alpha) * input[, 1]
}

#-----------------------------------------------------------------------------#
#' Causal estimands
#' @export
#-----------------------------------------------------------------------------#

estimands <- function(potential_outcomes, alpha){
  
  # --- Individual average potential outcomes --- #  
  ind_avg_po <- lapply(potential_outcomes, function(x) ybar_ij_aalpha(x, alpha))
  
  # --- Group average potential outcomes --- #                
  grp_avg_po <- lapply(ind_avg_po, function(x) apply(x, 2, mean))
  grp_avg_po <- matrix(unlist(grp_avg_po), ncol=2, byrow = T)
  # --- Population average potential outcomes --- # 
  pop_avg_po <- apply(grp_avg_po, 2, mean)
  
  # --- Marginal potential outcomes --- #
  ind_marg_po <- lapply(ind_avg_po, function(x) ybar_ij_alpha(x, alpha))
  grp_marg_po <- unlist(lapply(ind_marg_po, mean))
  pop_marg_po <- mean(grp_marg_po)
  return(list(alpha = alpha, 
              pop_avg_po = pop_avg_po, 
              pop_marg_po = pop_marg_po))
}

#-----------------------------------------------------------------------------#
#' Calculate estimands
#' @export
#-----------------------------------------------------------------------------#

calc_estimands <- function(potential_outcomes, alphas){
  true_effects <- lapply(alphas, estimands, 
                         potential_outcomes = potential_outcomes)
  
  nn <- length(alphas)
  
  # Create a data.frame of estimands
  truth <- data.frame(alpha = rep(NA, nn), 
                      a0 = rep(NA, nn), 
                      a1 = rep(NA, nn), 
                      marg = rep(NA, nn))
  
  for (i in 1:length(true_effects)){
    truth[i, 1] <- true_effects[[i]]$alpha
    truth[i, 2:3] <- true_effects[[i]]$pop_avg_po
    truth[i, 4] <- true_effects[[i]]$pop_marg_po
  }
  
  return(truth)
}


#-----------------------------------------------------------------------------#
#' Simulate date
#' @export
#-----------------------------------------------------------------------------#

sim_data <- function(base_dt, potential_outcomes, parameters){

  N <- length(unique(base_dt$group))
  n <- nrow(base_dt)
  theta_fix <- parameters[-length(parameters)]
  rvar <- parameters[length(parameters)]
  ## hood random effect ##
  # neighborhood level random effect from N(0, variance = 1.0859)
  group_effect <- data.frame(group = 1:N, 
                             b = rnorm(N, mean = 0, sd = sqrt(rvar)))
  
  dt <- merge(base_dt, group_effect, by='group')
  dt <- within(dt, {
                    xsum = theta_fix[1] + theta_fix[2]*X1 + theta_fix[3]*X2 + b
                    B = rbinom(n, 1, p = plogis(xsum))
                    # Assign Treatments
                    A = ifelse(B == 0, 0, rbinom(n, 1, p=2/3))})
  
  # Calculate k = # treated in hood i excluding subject j 
  sim_dt <- ddply(dt, .(group), mutate, k = sum(A) - A)
  
  # Observe the outcome
  y <- numeric(n)
  for(ii in 1:nrow(sim_dt)){
    subj <- sim_dt[ii, ]
    y[ii] <- potential_outcomes[[as.character(subj$group)]][as.character(subj$id), 
                                           as.character(subj$k), 
                                           as.character(subj$A)]
  }
  
  sim_dt$y <- y

  return(sim_dt[c('y', 'X1', "X2", "A", 'B', 'group')]) 
}

#-----------------------------------------------------------------------------#
#' TBD
#' @export 
#-----------------------------------------------------------------------------#

replicate_sims <- function(base_dt, 
                           potential_outcomes, 
                           nsims,
                           parameters){
  
  out <- list(nsims)
  for(sim in 1:nsims){
    out[[sim]] <- sim_data(base_dt, potential_outcomes, parameters)
  }
  return(out)
}

#-----------------------------------------------------------------------------#
#' Interference Simulations
#' @export 
#-----------------------------------------------------------------------------#

interference_sims <- function(n, N, nsims, base.parameters, parameters, alphas){
  basedt  <- sim_units(n = n, N = N)
  this.po <- make_po(basedt, base.parameters)
  
  out <- list()
  out$truth <- calc_estimands(this.po, alphas)  
  out$sims  <- replicate_sims(base_dt = basedt, 
                              potential_outcomes = this.po, 
                              nsims = nsims, 
                              parameters = parameters)
  out$base  <- basedt
  return(out)
}
