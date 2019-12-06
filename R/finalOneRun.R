#' Building local models in parallel for all the time points
#'
#' Builds local models in parallel by detecting the number of cores (\eqn{k}) on the machine. This function builds models at \eqn{tau} time points parllely using all avialable cores.
#'
#' @param y binary response correspoding to \eqn{n} subjects
#' @param X \eqn{n} x \eqn{p} x \eqn{tau} - tensor (EEG data) for \eqn{n} subjects, \eqn{p} locations and \eqn{tau} time points
#' @param c_init initial value for \eqn{c}
#' @param zeta_init initial values for \eqn{\zeta}
#' @param nusq_init initial values for \eqn{\nu^2}
#' @param w_init initial value for \eqn{w} - the complexity parameter
#' @param v0 hyperparameter \eqn{v0}
#' @param a1 hyperparameter \eqn{a1}
#' @param a2 hyperparameter \eqn{a2}
#' @param q hyperparameter \eqn{q}
#' @param Nmcmc number of MCMC samples to generate
#' @param ind indices of MCMC samples to use after burnin ad thinning
#' @keywords finalOneRun()
#' @export
#' @examples finalOneRun()

#################################################################
### Building local models in parallel for all the time points ###
#################################################################

finalOneRun = function(y,
                       X, # n x p x tau - tensor (EEG data) for n subjects, p locations and tau time points  
                       c_init,
                       zeta_init,
                       nusq_init,
                       w_init,
                       v0,
                       a1,
                       a2,
                       q,
                       Nmcmc,
                       ind){
  p = dim(X)[2]
  tau = dim(X)[3]
  
  nCores = detectCores() # detect the number of cores
  cl = makeCluster(nCores) # specify cluster with the specified number of cores
  registerDoParallel(cl, cores = nCores) # initiate the cluster
  # Estimating local models in parallel
  Nsim = foreach(t = 1:tau, .combine=cbind, .inorder=TRUE, .errorhandling = 'pass',
                 .packages = c('coda', 'truncnorm', 'mvtnorm', 'MASS')) %dopar% {
                   
                   #fit = glm(y~0+X[, ,t], family=binomial(link="logit"))
                   #beta_init = t(fit$coefficients)
                   beta_init = t(rep(0, p))
                   
                   # estimating local model at time point t
                   sim = SpikeSlabEEG::priorT(y,
                                              X[, ,t],
                                              c_init,
                                              beta_init,
                                              zeta_init,
                                              nusq_init,
                                              w_init,
                                              v0,
                                              a1,
                                              a2,
                                              q,
                                              Nmcmc,
                                              ind)
                   
                   # computing estimates using the MCMC samples
                   sim$zeta[sim$zeta==v0] = 0
                   zeta = apply(sim$zeta, 2, sum)
                   
                   bTemp = sim$b
                   b = apply(bTemp, 2, mean)
                   
                   wTemp = sim$w
                   wMAP = mapEst(wTemp)
                   wMED = median(wTemp)
                   wMEAN = mean(wTemp)
                   
                   cTemp = sim$c
                   cMAP = mapEst(cTemp)
                   cMED = median(cTemp)
                   cMEAN = mean(cTemp)
                   
                   list(zeta = zeta,
                        b = b,
                        wMAP = wMAP,
                        wMED = wMED,
                        wMEAN = wMEAN,
                        cMAP = cMAP,
                        cMED = cMED,
                        cMEAN = cMEAN)
                 }
  stopCluster(cl) # Stop cluster
  
  zeta = do.call(rbind, Nsim['zeta',])
  b = do.call(rbind, Nsim['b',])
  wMAP = unlist(Nsim['wMAP',])
  wMED = unlist(Nsim['wMED',])
  wMEAN = unlist(Nsim['wMEAN',])
  cMAP = unlist(Nsim['cMAP',])
  cMED = unlist(Nsim['cMED',])
  cMEAN = unlist(Nsim['cMEAN',])
  
  # return the estimates of parameters obtained from the MCMC samples...
  # ... after building all the local models.
  list(zeta = zeta,
       b = b,
       wMAP = wMAP,
       wMED = wMED,
       wMEAN = wMEAN,
       cMAP = cMAP,
       cMED = cMED,
       cMEAN = cMEAN)
}