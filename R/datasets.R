#-----------------------------------------------------------------------------#
#' Vaccine Study Sample Data
#'
#' A sample dataset based on the simulations of a cholera vaccine trial 
#' in Heydrich-Perez et al. (2014) (\doi{10.1111/biom.12184})
#' except with 3000 individuals in 250 groups rather than 10000 in 500.
#' @docType data
#' @format a dataset with 6 variables and 3000 rows
#' \itemize{
#'    \item{Y}{the outcome (0 - no cholera; 1 - cholera)}
#'    \item{X1}{an individual's age (in decades)}
#'    \item{X2}{an individual's distance from river}
#'    \item{A}{an indicator of vaccination (0 - no vaccine; 1 - vaccine)}
#'    \item{B}{an indicator of participation (0 - did not participant in vaccine trial, 1 - did participate)}
#'    \item{group}{group membership}
#'  }
#' @references Perez-Heydrich, C., Hudgens, M. G., Halloran, M. E., Clemens, J. D., Ali, M., & Emch, M. E. (2014). 
#' Assessing effects of cholera vaccination in the presence of interference. Biometrics, 70(3), 731-741.
#' @name vaccinesim
#' @keywords datasets
#-----------------------------------------------------------------------------#
NULL

#-----------------------------------------------------------------------------#
#' Voting Contagion Experiment Data
#'
#' A dataset of a voting contagion experiment. See Nickerson (2008) for more
#' details. The variables used in the package vignette are documented here.
#'  
#' @docType data
#' @format a dataset with 21 variables and 7722 rows
#' \itemize{
#'    \item{family}{household ID}
#'    \item{denver}{1 = subject in Denver, 0 = Minneapolis}
#'    \item{treatment}{1 = voting encouragement, 2 = recycling message, 3 = not contacted}
#'    \item{reached}{1 = subject answered door, 0 = not}
#'    \item{hsecontact}{1 = household contacted by canvassers, 0 = not}
#'    \item{voted02p}{1 = voted in '02 primary, 0 = not}
#'    \item{party}{party affiliation}
#'    \item{age}{age}
#'    \item{gender}{gender}    
#'  }
#' @references Nickerson, D. W. (2008). Is voting contagious? Evidence from two field experiments. 
#'  American Political Science Review, 102(01), 49-57. \doi{10.1017/S0003055408080039}
#' @name voters
#' @keywords datasets
#-----------------------------------------------------------------------------#
NULL