#-----------------------------------------------------------------------------#
#' Simulate a dataset of n individuals into N groups
#' 
#' This is currently not very useful except to simulate a dataset like in 
#' Step 0 of Perez et al's simulation from which potential outcomes 
#' will be generated.
#'
#' @param N.groups the number of groups
#' @param n.units the number of subjects (total)
#' @param X1.distr distribution function of \code{X1} variable. Defaults to
#' \code{\link{rexp}}
#' @param X1.args arguments passed to \code{X1.distr} function. In Perez et al, 
#' X1 represents an individual's age in decades generated from an exponential
#' distribution with mean 20
#' @param ind.dist.distr distribution function of individual distances. Defaults
#' to \code{\link{rnorm}}
#' @param ind.dist.args arguments passed to \code{ind.dist.distr} function
#' @param grp.dist.distr distribution function of group level distances. Defaults
#' to \code{\link{rlnorm}}.
#' @param grp.dist.args arguments passed to \code{grp.dist.distr} function
#' @return a data frame 
#' @export
#-----------------------------------------------------------------------------#
sim_base_data <- function(N.groups, 
                      n.units,
                      X1.distr = 'rexp',
                      X1.args  = list(rate = 1/20),
                      ind.dist.distr = 'rnorm',
                      ind.dist.args  = list(mean=0, sd=0.05),
                      grp.dist.distr = 'rlnorm', 
                      grp.dist.args  = list(meanlog=0, sdlog = 0.75))
{ 
  N <- N.groups
  n <- n.units
  
  ## Randomly generate group level distances ##
  groupdt <- data.frame(group = 1:N, 
                        group_dist = do.call(grp.dist.distr, 
                                             args = append(list(n = N), 
                                                           grp.dist.args)))
  ## Randomly assign subjects to groups ##
  out <- data.frame(id = 1:n, group = sample(1:N, n, replace = T))
  out <- merge(out, groupdt, by = 'group')
  
  ## Randomly generate individual level distances ##
  ind_dist <- do.call(ind.dist.distr, args = append(list(n = n), ind.dist.args))
  
  ## Randomly generate X1 (age) variable ##
  X1_ <- do.call(X1.distr, args = append(list(n = n), X1.args))
  
  ## Finish up ##
  out <- within(out, {
    X1 = ifelse(X1_ > 100, 10, X1_/10)
    X2 = group_dist + ifelse(ind_dist < 0, 0, ind_dist)})
  
  return(out[ , c('id', 'group', 'X1', 'X2')])
}

#-----------------------------------------------------------------------------#
#' Generate Potential Outcomes
#' 
#' Create a list of arrays representing potential outcomes \eqn{y_{ij}(a, k)}, 
#' (where \eqn{a} is the treatment level and \eqn{k} is the number of other
#' subjects in a group receiving treatment).
#' 
#' Potential outcomes are generated as bernoulli random variables with expectation:
#' \deqn{logit^{-1}(X\psi)}{plogis(X \%*\% parameters)}. \eqn{X_{ij}} is the vector:
#' \eqn{(1, a, \alpha, X_{1ij}, X_{2ij}, a * \alpha)}{(1, a, \alpha, X_{1ij}, 
#' X_{2ij}, a * \alpha)} where \eqn{\alpha = (a + k)/n_i}.
#' 
#' @param base_data simulated data from \code{\link{sim_base_data}}
#' @param parameters vector of parameters used to generate potential outcomes. 
#' Defaults to parameters used in Perez et. al 2014.
#' @return a list with length(# of groups) where each entry is a n_i x n_i x 2 
#' array representing potential outcomes for each subject by # of other subject's
#' treated by 2 treatment levels {0, 1}
#' @export
#-----------------------------------------------------------------------------#

generate_po <- function(base_data, 
                        parameters = c(.5, -.788, -2.953, -0.098, -0.145, 0.351)){
  
  ## Generates a vector of potential outcomes for a given 
  ## treatment level and number of other subjects in group treated
  po_y_ij <- function(grpdt, a, k, parameters){
    n_i <- nrow(grpdt)
    alpha <- (a + k)/n_i
    X <- cbind(1, a, alpha, grpdt$X1, grpdt$X2, a*alpha)
    p <- plogis(X %*% parameters)
    
    if(length(p) != n_i){stop('Problem in po_y_ij(): length(p) != n_i')}
    
    return(rbinom(n_i, 1, p))
  }
  
  ## Generates a matrix of potential outcomes for all treatment levels in 
  ## and all numbers of other subjects treated (from 0 to n_i - 1)
  gen_po_ij <- function(grpdt, parameters){
    n_i <- nrow(grpdt)
    ids <- grpdt$id
    po <- array(dim = c(n_i, n_i, 2), 
                dimnames = list(id = ids, k = 0:(n_i-1), trt = c(0,1)))
    
    for(k in 0:(n_i-1)){
      po[ , as.character(k), 1] <- po_y_ij(grpdt, 0, k, parameters)
      po[ , as.character(k), 2] <- po_y_ij(grpdt, 1, k, parameters)
    }
    
    return(po)
  }
  
  ## Generate potential outcomes for each group ##
  po <- lapply(split(base_data, base_data$group), 
               gen_po_ij, parameters = parameters)
  return(po)
}

#-----------------------------------------------------------------------------#
#' Simulate data
#' 
#' Simulate a single dataset as in Step 1 of 2014 Perez et al.
#' 
#' @param base_dt dataset created from \code{\link{sim_base_data}}
#' @param potential_outcomes array created from \code{\link{generate_po}}
#' @param parameters parameters used to generate simulated data.
#' @return single dataset with \code{nrow(base_dt)} observations.
#' @export
#-----------------------------------------------------------------------------#

sim_interference_data <- function(base_dt, potential_outcomes, parameters){

  N <- length(unique(base_dt$group))
  n <- nrow(base_dt)
  theta_fix <- parameters[-length(parameters)]
  rse <- parameters[length(parameters)]
  ## hood random effect ##
  # neighborhood level random effect from N(0, variance = 1.0859)
  group_effect <- data.frame(group = 1:N, 
                             b = rnorm(N, mean = 0, sd = rse))
  
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
#' Interference Simulations
#' 
#' Combines \code{\link{sim_base_data}}, \code{\link{generate_po}}, and 
#' \code{\link{sim_interference_data}} into one function.
#' 
#' @param n number of units
#' @param N number of groups
#' @param nsims number of simulations
#' @param base.parameters passed to \code{\link{sim_base_data}}
#' @param parameters passed to \code{\link{sim_interference_data}}
#' @param alphas allocation (coverage) levels to calculate estimand with 
#' \code{\link{calc_estimands}}
#' @return a list with 3 elements: (1) the base dataset, (2) sims, a list of 
#' \code{nsims} datasets (3) truth: the estimands computed from the sims.
#' 
#' @export 
#-----------------------------------------------------------------------------#

sim_interference <- function(n, 
                             N, 
                             nsims, 
                             base.parameters = c(.5, -.788, -2.953, -0.098, -0.145, 0.351), 
                             parameters = c(0.2727, -0.0387, 0.2179, 1.0859), 
                             alphas){
  
  basedt  <- sim_base_data(n = n, N = N)
  this.po <- generate_po(basedt, base.parameters)
  
  out <- list()
  out$base  <- basedt
  out$truth <- calc_estimands(this.po, alphas)  
  out$sims  <- replicate(nsims, sim_interference_data(basedt, 
                                                      this.po, 
                                                      parameters),
                         simplify = FALSE)

  return(out)
}
