context("Interference functions")

test_that("Interference() works in various situations", {
  testdt  <- subset(vaccinesim, group %in% 1:10)
  testdt2 <- data.frame(x = c(1, 1, 1, 0), y = c(1, 0, 1, 0), g = c(1,1,2,2))

  allos  <- c(.35, .4)
  
  # Using GLMER with random effect: should not give error
  expect_error(interference(data = testdt,
                           allocations = allos,
                           propensity_integrand = 'logit_integrand',
                           formula = Y | A | B ~ X1 + (1|group) | group,
                           model_method = 'glmer',
                           method = 'simple'), NA)
  
  # Using GLMER with no random effect: should fail
  expect_error(interference(data = testdt,
                            allocations = allos,
                            propensity_integrand = 'logit_integrand',
                            formula = Y | A | B ~ X1 | group,
                            model_method = 'glmer',
                            method = 'simple'))
  
  # Using GLMER when random effect == 0: should fail
  expect_error(interference(data = testdt2,
                            allocations = allos,
                            propensity_integrand = 'logit_integrand',
                            formula = y | x ~ (1|g) | g,
                            model_method = 'glmer',
                            method = 'simple'))
  
  # Using GLM: should pass
  expect_error(interference(data = testdt,
                           allocations = allos,
                           propensity_integrand = 'logit_integrand',
                           formula = Y | A | B ~ X1 | group,
                           model_method = 'glm',
                           method = 'simple'), NA)
  
  # Using GLM without intercept: should pass
  expect_error(interference(data = testdt,
                           allocations = allos,
                           propensity_integrand = 'logit_integrand',
                           formula = Y | A | B ~ 1 | group,
                           model_method = 'glm',
                           method = 'simple'), NA)
  
  # Using GLM without intercept with simple variance estimation: should pass
  expect_error(interference(data = testdt,
                            allocations = allos,
                            propensity_integrand = 'logit_integrand',
                            formula = Y | A | B ~ 1 | group,
                            model_method = 'glm',
                            causal_estimation_options = list(variance_estimation = 'naive'),
                            method = 'simple'), NA)
  
  # Using GLM without intercept with nonsense variance estimation: should fail
  expect_error(interference(data = testdt,
                            allocations = allos,
                            propensity_integrand = 'logit_integrand',
                            formula = Y | A | B ~ 1 | group,
                            model_method = 'glm',
                            causal_estimation_options = list(variance_estimation = 'jfklas'),
                            method = 'simple'))
  
  # Using Oracle without parameters defined: should fail
  expect_error(interference(data = testdt,
                           allocations = allos,
                           propensity_integrand = 'logit_integrand',
                           formula = Y | A | B ~ X1 | group,
                           model_method = 'oracle',
                           method = 'simple'))

  # Using Oracle with parameters defined: should pass
  expect_error(interference(data = testdt,
                           allocations = allos,
                           propensity_integrand = 'logit_integrand',
                           formula = Y | A | B ~ X1 | group,
                           model_method = 'oracle',
                           model_options = list(fixed.effects = c(1,.5), random.effects = NULL),
                           method = 'simple'), NA) 
})