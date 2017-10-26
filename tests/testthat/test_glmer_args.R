

context("Pass args to glmer")

testdt  <- subset(vaccinesim, group %in% 1:10) 
allos  <- c(.35, .4)

set.seed(22)
inf_fit <- interference(
  data = testdt,
  allocations = allos,
  propensity_integrand = 'logit_integrand',
  formula = Y | A | B ~ X1 + (1|group) | group,
  causal_estimation_options = list(
    variance_estimation = "naive"
  ),
  model_method = 'glmer',
  method = 'simple')

set.seed(22)
glmer_fit <- lme4::glmer(
  data = testdt,
  family = binomial,
  formula =  B ~ X1 + (1|group)
)

inf_fit_model <- inf_fit$models$propensity_model

set.seed(33)
inf_fit_start <- interference(
  data = testdt,
  allocations = allos,
  propensity_integrand = 'logit_integrand',
  formula = Y | A | B ~ X1 + (1|group) | group,
  causal_estimation_options = list(
    variance_estimation = "naive"
  ),
  model_options = list(
    family = binomial,
    start = list(
      fixef = c(0.5, -0.15),
      theta=1
    )
  ),
  model_method = 'glmer',
  method = 'simple')

set.seed(33)
inf_fit_nostart <- interference(
  data = testdt,
  allocations = allos,
  propensity_integrand = 'logit_integrand',
  causal_estimation_options = list(
    variance_estimation = "naive"
  ),
  formula = Y | A | B ~ X1 + (1|group) | group,
  model_options = list(
    family = binomial 
  ),
  model_method = 'glmer',
  method = 'simple')

set.seed(33)
glmer_fit_start <- lme4::glmer(
  data = testdt,
  family = binomial,
  formula =  B ~ X1 + (1|group),
  start = list(
    fixef = c(0.5, -0.15),
    theta=1
  )
)

inf_fit_start_model <- inf_fit_start$models$propensity_model
inf_fit_nostart_model <- inf_fit_nostart$models$propensity_model


 
tolerance <-  1e-12

testthat::test_that( 
  "Starting values are appropriately passed to glmer",
  {
    
    testthat::expect_equal(
      lme4::getME(inf_fit_model, c("fixef", "theta")),
      lme4::getME(glmer_fit, c("fixef", "theta")),
      tol = tolerance
    )
    
    testthat::expect_equal(
      lme4::getME(inf_fit_start_model, c("fixef", "theta")),
      lme4::getME(glmer_fit_start, c("fixef", "theta")),
      tol = tolerance
    )
    
    ## The test should pass because these should be different
    testthat::expect_false(isTRUE(all.equal(
      lme4::getME(inf_fit_nostart_model, c("fixef", "theta")),
      lme4::getME(glmer_fit_start, c("fixef", "theta")),
      tol = tolerance
    )))
  }
)


set.seed(55)
inf_fit_nAGQ <- interference(
  data = testdt,
  allocations = allos,
  propensity_integrand = 'logit_integrand',
  formula = Y | A | B ~ X1 + (1|group) | group,
  causal_estimation_options = list(
    variance_estimation = "naive"
  ),
  model_options = list(
    family = binomial ,
    nAGQ =5
  ),
  model_method = 'glmer',
  method = 'simple')

set.seed(55)
inf_fit_nonAGQ <- interference(
  data = testdt,
  allocations = allos,
  propensity_integrand = 'logit_integrand',
  formula = Y | A | B ~ X1 + (1|group) | group,
  causal_estimation_options = list(
    variance_estimation = "naive"
  ),
  model_options = list(
    family = binomial
  ),
  model_method = 'glmer',
  method = 'simple')

set.seed(55)
glmer_fit_nAGQ <- lme4::glmer(
  data = testdt,
  family = binomial,
  formula =  B ~ X1 + (1|group),
  nAGQ =5
)


inf_fit_nAGQ_model <- inf_fit_nAGQ$models$propensity_model
inf_fit_nonAGQ_model <- inf_fit_nonAGQ$models$propensity_model
glmer_fit_nAGQ

tolerance <-  1e-12

testthat::test_that( 
  "Starting values are appropriately passed to glmer",
  {
    
    testthat::expect_equal(
      lme4::getME(inf_fit_nAGQ_model, c("fixef", "theta")),
      lme4::getME(glmer_fit_nAGQ, c("fixef", "theta")),
      tol = tolerance
    )
     
    
    ## The test should pass because these should be different
    testthat::expect_false(isTRUE(all.equal(
      lme4::getME(inf_fit_nAGQ_model, c("fixef", "theta")),
      lme4::getME(inf_fit_nonAGQ_model, c("fixef", "theta")),
      tol = tolerance
    )))
  }
)
