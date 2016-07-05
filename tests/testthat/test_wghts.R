context("Weight Calculations")

test_that("weights calculations equal", {
  XXX <- cbind(rep(1, 10), seq(0, 1, length = 10))
  AAA <- rep(c(0, 1), 5)
  fff <- c(.1, .1)
  aaa <- c(.1, .5, .9)
  GGG <- sort(rep(1:5, 2))

  mmm2 <- matrix(c(0.616547184192009, 1.712631, 0.616547184192009, 
                   0.614796738538913, 1.707769, 0.614796738538913, 
                   0.613057773666495, 1.702938, 0.613057773666495,
                   0.611330209096351, 1.698139, 0.611330209096351,
                   0.609613965009724, 1.693372, 0.609613965009724), byrow = T, ncol = 3,
                dimnames = list(unique(GGG), aaa))
  
  mmm3 <- array(c(-0.07902959, -0.01322103, 0.01576307,
                  -0.07851113, -0.03058825, 0.01598160,
                  -0.07799631, -0.04772718, 0.01619665,
                  -0.07748510, -0.06464023, 0.01640826,
                  -0.07697746, -0.08132978, 0.01661645,
                  -0.2195266, -0.03672508, 0.04378631,
                  -0.2180865, -0.08496735, 0.04439334,
                  -0.2166564, -0.13257550, 0.04499070,
                  -0.2152364, -0.17955620, 0.04557849,
                  -0.2138263, -0.22591605, 0.04615682,
                  -0.07902959, -0.01322103, 0.01576307,
                  -0.07851113, -0.03058825, 0.01598160,
                  -0.07799631, -0.04772718, 0.01619665,
                  -0.07748510, -0.06464023, 0.01640826,
                  -0.07697746, -0.08132978, 0.01661645),
                 dim = c(3, 5, 3),
                 dimnames = list(1:3, unique(GGG), aaa))
  mmm3 <- aperm(mmm3, c(2,1,3))
  
  # Checking weight calculations
  expect_equal(wght_calc(integrand = logit_integrand, 
                         allocation = .9, 
                         X = XXX[1:2, ], A = AAA[1:2], 
                         parameters = fff, 
                         randomization = .5), 0.462462731646267)
  
  expect_equal(wght_calc(integrand = logit_integrand, 
                         allocation = .9, 
                         X = XXX[1:2, ], A = AAA[1:2], 
                         parameters = c(fff, 1), 
                         randomization = .5), 0.492904177726425)

  expect_equal(wght_calc(integrand = logit_integrand, 
                         allocation = .9, 
                         X = XXX[1:2, ], A = AAA[1:2], 
                         parameters = c(fff, 1), 
                         randomization = .5), 0.492904177726425)
  
  # Checking derivative calculations
  expect_equal(wght_deriv_calc(integrand = logit_integrand, 
                               allocation = .3, 
                               X = XXX[1:2, ], A = AAA[1:2], 
                               parameters = c(1, .5, 1), 
                               randomization = .5), c(-0.1277388445291851,
                                                      -0.0280714662250347,
                                                      0.0947948595448454))
  
  # Checking weight matrix calculations
  expect_equal(wght_matrix(integrand = logit_integrand,
                           allocations = aaa,
                           X = XXX, A = AAA, G = GGG,
                           parameters = c(fff, 5),
                           randomization = .5,
                           runSilent = TRUE), mmm2,
               tolerance = 1e-6)
  
  # Checking weight derivative calculations
  expect_equal(wght_deriv_array(integrand = logit_integrand,
                           allocations = aaa,
                           X = XXX, A = AAA, G = GGG,
                           parameters = c(fff, 5),
                           randomization = .5,
                           runSilent = TRUE), mmm3,
               tolerance = 1e-6, check.attributes = FALSE)
  
})