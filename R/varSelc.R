#' Performing two-stage variable selection
#'
#' Performing two-stage variable selection using the MCMC estimates of local Bayesian modeling
#'
#' @param b \eqn{tau} x \eqn{p} matrix of estimates of beta
#' @param method thresholding method to be used 'Zcut' or 'Wcut', defaults to 'Zcut'
#' @param wEst estimates of complexity parameters required if method == 'Wcut'
#' @param bf.thres threshold to be used to compare Bayes factor
#' @keywords varSelc()
#' @export
#' @examples varSelc()

###############################################
### Performing two-stage variable selection ###
###############################################
varSelc = function(b,
                   method = 'Zcut',
                   wEst = NULL,
                   bf.thres){
  ######################################
  ### First-stage variable selection ###
  ######################################
  if(method == 'Zcut') n0b = SpikeSlabEEG::n0beta(b, method = "Zcut")
  if(method == 'Wcut') n0b = SpikeSlabEEG::n0beta(b, method = "Wcut", w = wEst)
  
  #######################################
  ### Second stage variable selection ###
  #######################################
  b.bin = matrix(1, nrow = tau, ncol = p)
  b.bin[n0b==0] = 0
  
  maxBF = numeric()
  for(l in 1:p){
    if(sum(b.bin[,l])==0){
      maxBF[l] = -10
    } else if(sum(b.bin[,l])==tau){
      maxBF[l] = bf.thres
    } else maxBF[l] = BayesTestStreak::bayes.factor.function(b.bin[,l],
                                                             log.K = seq(0, 3, by = 0.1))$max.log.BF
  }
  
  bInd = which(maxBF>=bf.thres)
  
  list(n0b = n0b,
       bInd = bInd)
}