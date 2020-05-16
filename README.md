# SpikeSlabEEG
R package for "Mohammed, S., Dey, D.K. and Zhang, Y., 2019. Bayesian variable selection using spike‐and‐slab priors with application to high dimensional electroencephalography data by local modelling. Journal of the Royal Statistical Society: Series C (Applied Statistics), 68(5), pp.1305-1326."

```

################################################################
### Code to run the local Bayesian modeling on ...           ###
### ... the EEG data using the package SpikeSlabEEG.         ###
### Contains:                                                ###
### 0. Loads EEG data and divides it into training and test  ###
### 1. Builds local models in parallel                       ###
### 2. Uses MCMC estimates to do variable selection          ###
### 3. Predicts the alcoholic status of subjects in test set ###
################################################################

setwd('./Downloads') # set working directory where the package folder is stored
devtools::install('SpikeSlabEEG') # install the package SpikeSlabEEG
require(SpikeSlabEEG)

#loading required packages
require(coda)
require(truncnorm)
require(mvtnorm)
require(MASS)
require(doParallel)
require(BayesTestStreak) # install from https://github.com/bayesball/BayesTestStreak


#################################################################
### 0. Loads EEG data and divides it into training and test   ###
#################################################################
X = SpikeSlabEEG::X # EEG data
ind64to57 = SpikeSlabEEG::ind64to57 # location indices for which distance matrix is available
X = X[,ind64to57,] # update EEG data for locations above
y = t(SpikeSlabEEG::y) # binary response

# dimensions of EEG data
tau = dim(X)[1]
p = dim(X)[2]
n = dim(X)[3]

# scale data at subject level by its Frobenius norm
X_int = sapply(1:n, function(i) sqrt(sum(X[,,i]^2)))
X_sc = array(NA, dim = dim(X))
for(i in 1:n) X_sc[,,i] = X[,,i]/X_int[i]

trn.ind = sample.int(n, size = 100) # sample indices from training data

###################################################
### 1. Builds local models in parallel          ###
###################################################
c_init = 1
w_init = 0
zeta_init = rep(0, p)
nusq_init = rep(1, p)

q = 10
a1 = 5
a2 = 20
v0 = 0.005

Nmcmc = 2000
burnin = 1000
ind = seq(burnin+1, Nmcmc, by=10)

modelFit = finalOneRun(y[trn.ind], aperm(X_sc[,,trn.ind], c(3,2,1)),
                       c_init, zeta_init, nusq_init, w_init,
                       v0, a1, a2, q,
                       Nmcmc, ind)

##############################################################
### 2. Uses MCMC estimates to do variable selection        ###
##############################################################
vSelc = varSelc(modelFit$b, method = 'Wcut', wEst = modelFit$wMED, bf.thres = 1)

#################################################################
### 3. Predicts the alcoholic status of subjects in test set  ###
#################################################################
y.pred = pred(n0b = vSelc$n0b, bInd = vSelc$bInd,
              X.test = aperm(X_sc[,,-trn.ind], c(3,2,1)), cEst = modelFit$cMED)

# function to estimates TPR, FPR and Prediction error
rates = function(x, x.hat){
  tpr = sum((x+x.hat) == 2)/sum(x)
  fpr = sum((x-x.hat) == -1)/(length(x) - sum(x))
  pe = sum((x+x.hat) == 1)/length(x)
  
  return(c(tpr, fpr, pe))
}

# compute prediction error
err.rates = rates(y[-trn.ind], y.pred)


```
