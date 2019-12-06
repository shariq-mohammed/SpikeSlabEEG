#' Compute MAP estimate 
#'
#' Given a sequence of MCMC samples this function computes the MAP estimate of the posterior distribution.
#'
#' @param chain sequence of MCMC samples
#' @param burnin Number of initial samples to drop
#' @param thin Number of samples to be dropped between two consecutive selections in the chain
#' @keywords mapEst()
#' @export
#' @examples mapEst()

#############################################
######## Compute MAP estimate ###############
#############################################

mapEst = function(chain, # sequence of MCMC samples
                  burnin = 0, # Number of initial samples to drop
                  thin = 1 # Number of samples to be dropped between two consecutive selections in the chain
                  ){
  den = density(chain[seq(burnin+1, length(chain), by = thin)]) # Create density
  ind = which.max(den$y) # Find the maximum density
  return(den$x[ind])
}