# SpikeSlabEEG

R package for "[Mohammed, S.](shariq-mohammed.github.io), Dey, D.K. and Zhang, Y., 2019. Bayesian variable selection using spike-and-slab priors with application to high dimensional electroencephalography data by local modelling. Journal of the Royal Statistical Society: Series C (Applied Statistics), 68(5), pp.1305-1326. [https://doi.org/10.1111/rssc.12369](https://doi.org/10.1111/rssc.12369)"

Code to run the local Bayesian modeling on the EEG data using the package SpikeSlabEEG.   

## Contents
1. Load EEG data and divide it into training and test.
2. Build local models in parallel.
3. Use MCMC estimates to do variable selection.
4. Predict the alcoholic status of subjects in test set.

```
# install the package (devtools package needed)
if(!require(devtools)) install.packages("devtools")
devtools::install_github('shariq-mohammed/SpikeSlabEEG')
library(SpikeSlabEEG)

#loading required packages
library(coda)
library(truncnorm)
library(mvtnorm)
library(MASS)
library(doParallel)
library(BayesTestStreak) # install from https://github.com/bayesball/BayesTestStreak
```

### Load EEG data and divide it into training and test
Load EEG Data
```
X = SpikeSlabEEG::X
```
Load location indices for which distance matrix is available
```
ind64to57 = SpikeSlabEEG::ind64to57
```

Update EEG data for locations above
```
X = X[,ind64to57,]
```

Load binary response
```
y = t(SpikeSlabEEG::y)
```
Dimensions of EEG data
```
tau = dim(X)[1]
p = dim(X)[2]
n = dim(X)[3]
```

Scale data at subject level by its Frobenius norm
```
X_int = sapply(1:n, function(i) sqrt(sum(X[,,i]^2)))
X_sc = array(NA, dim = dim(X))
for(i in 1:n) X_sc[,,i] = X[,,i]/X_int[i]
```

Sample indices for training data
```
trn.ind = sample.int(n, size = 100)
```


### Build local models in parallel
Initialize and choose hyperparamters
```
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
```

Fit the model
```
modelFit = finalOneRun(y[trn.ind], aperm(X_sc[,,trn.ind], c(3,2,1)),
                       c_init, zeta_init, nusq_init, w_init,
                       v0, a1, a2, q, Nmcmc, ind)
```

### Use MCMC estimates to do variable selection
```
vSelc = varSelc(modelFit$b, method = 'Wcut', wEst = modelFit$wMED, bf.thres = 1)
```

### Predict the alcoholic status of subjects in test set
```
y.pred = pred(n0b = vSelc$n0b, bInd = vSelc$bInd,
              X.test = aperm(X_sc[,,-trn.ind], c(3,2,1)), cEst = modelFit$cMED)
```

Function to estimates true positive rate, false positive rate, and prediction error
```
rates = function(x, x.hat){
  tpr = sum((x+x.hat) == 2)/sum(x)
  fpr = sum((x-x.hat) == -1)/(length(x) - sum(x))
  pe = sum((x+x.hat) == 1)/length(x)
  
  return(c(tpr, fpr, pe))
}
```
Compute prediction error
```
err.rates = rates(y[-trn.ind], y.pred)
```
