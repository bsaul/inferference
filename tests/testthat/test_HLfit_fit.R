

context("Fitting spaMM::HLfit() model works inside of interference()")

testdt  <- subset(vaccinesim, group %in% 1:10) 
allos  <- c(.35, .4)

set.seed(22)
inf_fit <- interference(
  data = testdt,
  allocations = allos,
  propensity_integrand = 'logit_integrand',
  formula = Y | A | B ~ X1 + (1|group) | group,
  causal_estimation_options = list(
    variance_estimation = "robust"
  ),
  model_options = list(
    family = binomial ,
    rand.family = gaussian()
  ),
  model_method = 'HLfit',
  method = 'Richardson'
)

inf_fit_model <- inf_fit$models$propensity_model

set.seed(22)
model_fit <- spaMM::HLfit(
  data = testdt,
  family = binomial(),
  rand.family = gaussian(),
  formula =  B ~ X1 + (1|group)
)


set.seed(22)
inf_fit_glmer <- interference(
  data = testdt,
  allocations = allos,
  propensity_integrand = 'logit_integrand',
  formula = Y | A | B ~ X1 + (1|group) | group,
  causal_estimation_options = list(
    variance_estimation = "robust"
  ),
  model_method = 'glmer',
  model_options = list(
    family = binomial ,
    nAGQ =5
  ),
  method = 'Richardson')

inf_fit_model_glmer <- inf_fit_glmer$models$propensity_model




tolerance <-  1e-12

testthat::test_that( 
  "HLfit model has the same parameters as when fit independently",
  {
    
    testthat::expect_equal(
      c(inf_fit_model$fixef, sqrt(inf_fit_model$lambda)),
      c(model_fit$fixef, sqrt(model_fit$lambda)),
      tol = tolerance
    )
  }
)


tolerance_fixef <-  1e-2
testthat::test_that( 
  "HLfit model has similar fixed effect parameters as glmer model",
  {
    
    testthat::expect_equal(
      inf_fit_model$fixef,
      lme4::getME(inf_fit_model_glmer, c("fixef")),
      tol = tolerance_fixef,
      check.attributes = FALSE
    )
  }
)

tolerance_ranef <- 1e-1
testthat::test_that( 
  "HLfit model has fairly similar random effect parameters as glmer model",
  {
    
    testthat::expect_equal(
      sqrt(inf_fit_model$lambda[[1]]),
      lme4::getME(inf_fit_model_glmer, c("theta")),
      tol = tolerance_ranef,
      check.attributes = FALSE
    )
  }
)

ests_HLfit <- inf_fit$estimates
ests_glmer <- inf_fit_glmer$estimates

testthat::test_that( 
  "HLfit model and glmer model have similar point estimates",
  {
    
    testthat::expect_equal(
      ests_HLfit[ests_HLfit$alpha1 == 0.35, ],
      ests_glmer[ests_glmer$alpha1 == 0.35, ],
      tol = 0.05,
      check.attributes = FALSE
    )
  }
)