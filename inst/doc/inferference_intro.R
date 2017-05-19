## ----setup, echo = FALSE, eval = TRUE, cache=FALSE-----------------------
library(knitr)
opts_chunk$set(prompt=TRUE)
options(continue ="+ ")

## ----library, echo = TRUE, eval = TRUE-----------------------------------
library(inferference)
head(vaccinesim)

## ----example1, echo = TRUE, eval = TRUE, results = 'hide', cache = FALSE----
example1 <- interference(
    formula = Y | A | B ~ X1 + X2 + (1|group) | group, 
    allocations = c(.3, .45,  .6), 
    data = vaccinesim, 
    randomization = 2/3,
    method = 'simple')

## ----example1_summary, echo = TRUE, eval = TRUE--------------------------
print(example1)

## ---- echo = TRUE, eval = TRUE-------------------------------------------
direct_effect(example1, .3)
ie(example1, .3)

## ----diagnostic, echo = TRUE, eval = TRUE, fig.width = 4.5, fig.height = 4.5----
diagnose_weights(example2, allocations = .5, breaks = 30)

## ----voters_data, echo = TRUE, eval = TRUE, cache = FALSE----------------
voters <-  within(voters, {
    treated     = (treatment == 1 & reached == 1) * 1 
    c_age       = (age - mean(age))/10 
   })
reach_cnt <- tapply(voters$reached, voters$family, sum)
voters <- voters[!(voters$family %in% names(reach_cnt[reach_cnt > 1])), ]
voters <- voters[voters$hsecontact == 1, ]

## ----voter_coef, echo=TRUE, eval = TRUE, cache = FALSE-------------------
voters1 <- do.call(rbind, by(voters, voters[, 'family'], function(x) x[1, ]))
coef.voters <- coef(glm(reached ~ c_age, data = voters1, 
		               family = binomial(link = 'logit')))

## ----propensity_fixed, echo = TRUE, eval = FALSE-------------------------
#  fixed_propensity <- function(b){
#  	return(0.5 * dnorm(b))
#  }

## ----household_propensity, echo =TRUE, eval = TRUE-----------------------
household_propensity <- function(b, X, A, 
                                 parameters, 
                                 group.randomization = .5){
  if(!is.matrix(X)){
    X <- as.matrix(X)
  }   
  if(sum(A) == 0){ 
    pr <-  group.randomization
  } else { 
    X.1 <- X[1 ,]; A.1 <- A[1] 
    h   <- plogis(X.1 %*% parameters)
    pr  <-  group.randomization * dbinom(A.1, 1, h)
  }   
  out <- pr * dnorm(b) 
  out
}

## ----example3, echo=TRUE, eval = TRUE, results = 'hide', cache = FALSE----
example3 <- interference(
  formula = voted02p | treated | reached ~ c_age | family,
  propensity_integrand = 'household_propensity',
  data = voters,
  model_method = 'oracle', 
  model_options = list(fixed.effects = coef.voters, random.effects = NULL),
  allocations   = c(0, .5),
  integrate_allocation = FALSE,
  causal_estimation_options = list(variance_estimation = 'robust'),
  conf.level = .9)

## ----example3_results, echo=TRUE, eval = TRUE----------------------------
ie(example3, .5, 0)[ , c('estimate', 'conf.low', 'conf.high')]

## ----example4, echo= TRUE, eval = TRUE, results = 'hide', cache = FALSE----
example4 <- interference(
  formula = voted02p | treated | reached ~ 1 | family,
  propensity_integrand = 'household_propensity',
  data = voters,
  model_method = 'oracle',
  model_options = list(fixed.effects = 0, random.effects = NULL),
  allocations   = c(0, .5),
  integrate_allocation = FALSE,
  causal_estimation_options = list(variance_estimation = 'naive'),
  conf.level = .9)

## ----example4_results, echo=TRUE, eval = TRUE----------------------------
ie(example4, .5, 0)[ , c('estimate', 'conf.low', 'conf.high')]

## ----example4_weights, echo = TRUE---------------------------------------
G <- tapply(voters[1:12, 'treated'], voters[1:12, 'family'], sum)
W <- head(example4[["weights"]])[, 2]
cbind(G, W)

## ----compare_weights, echo = TRUE, eval = TRUE---------------------------
compare_weights <- function(n, alpha = .5, h = .5){
  pi  <- rep(alpha, n)
  PrA <- rep(h, n)
  c(w1 = prod(pi)/prod(PrA),
    w2 = 1/prod(PrA/pi),
    w3 = 1/exp(sum(log(PrA/pi))))  
}
n <- c(50, 100, 500, 1074, 1075, 10000)
cbind(n, t(sapply(n, FUN = compare_weights)))

