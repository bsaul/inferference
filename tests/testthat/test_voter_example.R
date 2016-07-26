context("voters example")

test_that("Weights are computed correctly in voters example", {
  
  household_propensity <- function(b, ## b must the first argument ##
                                   X, A,
                                   parameters,
                                   group.randomization = .5){
    ## Necessary to be sure matrix multiplication works ##
    if(!is.matrix(X)){
      X <- as.matrix(X)
    }
    if(sum(A) == 0){ # No one in the household received treatment
      pr <-  group.randomization
    } else { # one subject received treatment
      X.1 <- X[1 ,]; A.1 <- A[1] # Use the first subject's values 
      h <- plogis(X.1 %*% parameters)
      pr <- group.randomization * dbinom(A.1, 1, h)
    }
    out <- pr * dnorm(b) # dnorm integrates to 1
    out
  }
  
  XXX  <- c(1, 1)
  AAA1 <- c(0, 0)
  AAA2 <- c(1, 0)
  
  
  expect_equal(wght_calc(integrand = household_propensity, 
            allocation = .5, 
            X = XXX, A = AAA2, 
            integrate_allocation = FALSE,
            parameters = c(0)), 1)
  
  expect_equal(wght_calc(integrand = household_propensity, 
                         allocation = .5, 
                         X = XXX, A = AAA1, 
                         integrate_allocation = FALSE,
                         parameters = c(0)), .5)
  
  # Using integrate_allocation = TRUE should return wrong answer
  expect_equal(wght_calc(integrand = household_propensity, 
                         allocation = .5, 
                         X = XXX, A = AAA1, 
                         integrate_allocation = TRUE,
                         parameters = c(0)), 2)
  
  
  
  
})