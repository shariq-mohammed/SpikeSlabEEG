#' MCMC algorithm for the estimation of local model
#'
#' Gibbs sampling algorithm for the estiamtion of a local model
#'
#' @param y binary response correspoding to \eqn{n} subjects
#' @param X data at time \eqn{t} - dimension \eqn{n} x \eqn{L}
#' @param c_init initial value for \eqn{c}
#' @param beta_init initial values for \eqn{\beta}
#' @param zeta_init initial values for \eqn{\zeta}
#' @param nusq_init initial values for \eqn{\nu^2}
#' @param w_init initial value for \eqn{w} - the complexity parameter
#' @param v0 hyperparameter \eqn{v0}
#' @param a1 hyperparameter \eqn{a1}
#' @param a2 hyperparameter \eqn{a2}
#' @param q hyperparameter \eqn{q}
#' @param Nmcmc number of MCMC samples to generate
#' @param ind indices of MCMC samples to use after burnin ad thinning
#' @keywords priorT()
#' @export
#' @examples priorT()

########################################################
### MCMC algorithm for the estimation of local model ###
########################################################

priorT = function(y, # binary response
                  X, # data at time t
                  c_init, # initial values for c
                  beta_init, # initial values for beta
                  zeta_init, # initial values for zeta
                  nusq_init, # initial values for nusq
                  w_init, # initial values for w
                  v0, # v0 hyperparameter
                  a1, # a1 hyperparameter
                  a2, # a2 hyperparameter
                  q, # q hyperparameter
                  Nmcmc, # number of MCMC smaples to generate
                  ind # indices of MCMC samples to use after burnin and thinning
                  ){
  nc = ncol(X)
  nr = nrow(X)
  
  c = c_init
  b = beta_init
  zeta = zeta_init
  nusq = nusq_init
  w = w_init
  
  # Storing all the required parameters from the MCMC samples
  c_store = matrix(nrow = Nmcmc, ncol = 1)
  zeta_store = matrix(nrow = Nmcmc, ncol = nc)
  b_store = matrix(nrow = Nmcmc, ncol = nc)
  #nusq_store = matrix(nrow = Nmcmc, ncol = nc)
  w_store = matrix(nrow = Nmcmc, ncol = 1)
  
  # start Gibbs sampling
  for(i in 1:Nmcmc){
    if(i %% 1000 == 0) print(i)
    
    # update z
    zMean = crossprod(t(X), t(b))
    z = (y*rtruncnorm(1, 0, Inf, zMean, sd = sqrt(c))+
           (1-y)*rtruncnorm(1, -Inf, 0, zMean, sd = sqrt(c)))
    
    # update beta
    symX = crossprod(X, X/c)
    GMinv = diag(1/(nusq*zeta))
    bVar = chol2inv(chol(symX+GMinv))
    temp1 = crossprod(X, z/c)
    bMean = crossprod(bVar, temp1)
    
    b = rmvnorm(1, mean = bMean, sigma = bVar)
    b_store[i,] = b
    
    # update zeta
    tempExp = exp(-(b^2)/(2*nusq))
    pr1 = (1-w)*(tempExp^(1/v0))/sqrt(v0)
    pr2 = w*tempExp
    temp.probs = pr2/(pr1+pr2)
    temp.probs[is.na(temp.probs)] = 0.5
    zeta = rbinom(nc, 1, prob = temp.probs)
    zeta[zeta==0] = v0
    zeta_store[i,] = zeta
    
    # update nusq
    nusq = 1/rgamma(nc, shape = a1+0.5, rate = a2+((b^2)/(2*zeta)))
    #nusq_store[i,] = nusq
    
    # update w
    w = rbeta(1, 1+sum(zeta==1), 1+sum(zeta==v0))
    w_store[i,] = w
    
    # update c
    c = 1/rgamma(1, shape = (nr+q)/2, rate = (q+sum((z-zMean)^2))/2)
    c_store[i,] = c
  }
  
  # return samples of beta, zeta, w and c as they will be needed for ...
  # ... variable selection and prediction
  list(zeta = zeta_store[ind,],
       c = c_store[ind,],
       b = b_store[ind,],
       #nusq = nusq_store[ind,],
       w = w_store[ind,])
}