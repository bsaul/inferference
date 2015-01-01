context("Interference functions")

test_that("Interfence() works in various situations", {
  testdt <- subset(interference_sample, group %in% 1:10)
  allos  <- c(.35, .4)
  
  # Using GLMER with random effect: should not give error
  expect_that(interference(data = testdt,
                             allocations = allos,
                             outcome = 'y',
                             treatment = 'A',
                             propensityB = 'B',
                             group = 'group',
                             propensity_formula = B ~ X1 + (1|group),
                             method = 'simple'), not(throws_error()))
  
  # Using GLMER with no random effect: should fail
  expect_error(interference(data = testdt,
                            allocations = allos,
                            outcome = 'y',
                            treatment = 'A',
                            propensityB = 'B',
                            group = 'group',
                            propensity_formula = B ~ X1,
                            method = 'simple'))
  # Using GLM: should pass
  expect_that(interference(data = testdt,
                           allocations = allos,
                           outcome = 'y',
                           treatment = 'A',
                           propensityB = 'B',
                           group = 'group',
                           propensity_formula = B ~ X1,
                           model_method = 'glm',
                           method = 'simple'), not(throws_error()))
  
  # Using GLM without intercept: should pass
  expect_that(interference(data = testdt,
                           allocations = allos,
                           outcome = 'y',
                           treatment = 'A',
                           propensityB = 'B',
                           group = 'group',
                           propensity_formula = B ~ 1,
                           model_method = 'glm',
                           method = 'simple'), not(throws_error()))
  
  # Using Oracle without parameters defined: should fail
  expect_that(interference(data = testdt,
                           allocations = allos,
                           outcome = 'y',
                           treatment = 'A',
                           propensityB = 'B',
                           group = 'group',
                           propensity_formula = B ~ X1,
                           model_method = 'oracle',
                           method = 'simple'), throws_error())

  # Using Oracle with parameters defined: should pass
  expect_that(interference(data = testdt,
                           allocations = allos,
                           outcome = 'y',
                           treatment = 'A',
                           propensityB = 'B',
                           group = 'group',
                           propensity_formula = B ~ X1,
                           model_method = 'oracle',
                           model_options = list(fixed.effects = c(1,.5), random.effects = NULL),
                           method = 'simple'), not(throws_error())) 
})