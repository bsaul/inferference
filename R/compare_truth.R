#-----------------------------------------------------------------------------#
#' Compare simulation estimates to truth 
#'
#' @export
#-----------------------------------------------------------------------------#

compare_truth <- function(a1, 
                          trt1 = NA, 
                          a2 = NA, 
                          trt2 = NA, 
                          effect, 
                          marginal, 
                          list_of_simulation_estimates,
                          truth){
  
  # Pull together all the estimates
  estimates <- rbind.fill(lapply(list_of_simulation_estimates, function(x){
    calc_effect(x, a1, trt1, a2, trt2, effect=effect, marginal=marginal)
  }))
  
  # Find the truth and nothing but the truth
  tr <- truth
  
  this.truth <- ifelse(effect=='contrast', 
                       ifelse(marginal == TRUE, 
                              tr[tr$alpha == a1, 'marg'] - tr[tr$alpha == a2, 'marg'],
                              ifelse(trt2 == 1, 
                                     tr[tr$alpha == a1, 'a0'] - tr[tr$alpha == a2, 'a1'],
                                     tr[tr$alpha == a1, 'a0'] - tr[tr$alpha == a2, 'a0'])),
                       ifelse(marginal == TRUE, 
                              tr[tr$alpha == a1, 'marg'],
                              ifelse(trt1 == 0, tr[tr$alpha == a1, 'a0'],
                                     tr[tr$alpha == a1, 'a1'])))
  # Point estimate
  pe <- mean(estimates$point)
  
  # ASE
  ase <- mean(sqrt(estimates$variance))
  # ASE Coverage 
  ase_coverage <- mean(estimates$ll < this.truth & this.truth < estimates$ul)
  
  # ESE
  ese <- sd(estimates$point)

  # Bias
  bias <- mean(estimates$point) - this.truth
  
  # Output
  out <- list()
  out$summary <- data.frame(a1 = a1, trt1 = trt1, a2 = a2, trt2 = trt2, 
                            avg_estimate = pe, bias = bias, ase = ase,
                            ase_coverage = ase_coverage,
                            ese = ese)
  out$estimates <- estimates
  #write(paste(effect_name, round(pe, 3), round(bias, 3), round(ase, 4), signif(coverage, 2), 
  #            round(ese, 4), signif(coverage2, 2), sep = ' '), file = file, append= TRUE)
  return(out)
}

#-----------------------------------------------------------------------------#
#' Summarize comparisions to truth
#'
#' @export
#-----------------------------------------------------------------------------#

summarize_comparisons <- function(estimates, truth, alphas){
  
  ## Build list of arguments ##
  args1 <- expand.grid(a1 = alphas, t1 = c(0,1), a2 = NA, t2 = NA, 
              effect = c('outcome'), marginal = c(TRUE, FALSE))
  args1 <- args1[-which(args1$marginal == TRUE & args1$t1 == 1), ]
  args2 <- expand.grid(a1 = alphas, t1 = c(0,1), a2 = alphas, t2 = c(0,1), 
              effect = c('contrast'), marginal = c(TRUE, FALSE) )
  args2 <- args2[-which((args2$marginal == TRUE & (args2$t1 == 1 | args2$t2 == 1)) |
                          args2$a2 < args2$a1 | (args2$a1 == args2$a2 & args2$t1 == args2$t2) |
                          (args2$t1 == 1 & args2$t2 %in% c(0, 1))), ]
  args <- rbind(args1, args2)
  
  args <- within(args, {t1 <- ifelse(marginal == TRUE, NA, t1)
                        t2 <- ifelse(marginal == TRUE, NA, t2)})
  
  ## Output comparisons ##
  hold <- list()
#   out$summary_data <- data.frame(a1 = rep(NA, length = nrow(args)), 
#                                  trt1 = NA, a2 = NA, trt2 = NA, 
#                                  avg_estimate = NA, bias = NA, ase = NA,
#                                  ase_coverage = NA,
#                                  ese = NA)
  for(ii in 1:nrow(args)){
    hold[[ii]] <- compare_truth(a1 = args$a1[ii], trt1 = args$t1[ii], 
                               a2 = args$a2[ii], trt2 = args$t2[ii],
                               effect = args$effect[ii], 
                               marginal = args$marginal[ii], 
                               list_of_simulation_estimates = estimates,
                               truth = truth)$summary
    #out$final_data[[ii]] <- out[[ii]]$summary
  }
  
  out <- rbind.fill(hold)

  out <- within(out, {ef1 <- ifelse(is.na(a2), 'Y', 
                                  ifelse(is.na(trt2), 'OE',
                                        ifelse(a1 == a2, 'DE', 
                                               ifelse(trt1 == trt2, 'IE', 'TE'))))
                      ef2 <- paste0('(', a1, 
                                   ifelse(ef1 == 'Y' & !is.na(trt1), paste0(', ', trt1), ''), 
                                   ifelse(is.na(a2), '', paste0(', ', a2)), ')')
                      effect <- paste0(ef1, ef2)
                      ef1 <- factor(ef1, levels = c('Y', 'DE', 'IE', 'TE', 'OE'),
                                    ordered = TRUE)
                      })
  out <- out[order(out$ef1), ]

  return(out[ , c('effect', 'avg_estimate', 'bias', 'ase_coverage', 'ase', 'ese')])
}


  
