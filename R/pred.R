#' Prediction for subjects
#'
#' Prediction of the alcoholic status using the estimates from local Bayesian modeling and variable selection
#'
#' @param n0b \eqn{tau} x \eqn{p} matrix of estimates of beta as returned by varSelc function
#' @param bInd indices of selected locations after two-stage variable selection as returned by varSelc function
#' @param X.test test data with dimensions \eqn{m} x \eqn{p} x \eqn{tau}
#' @param cEst \eqn{tau} estimates of \eqn{c} for each local models
#' @keywords pred()
#' @export
#' @examples pred()

###############################
### Prediction for subjects ###
###############################
pred = function(n0b,
                bInd,
                X.test,
                cEst){
  m = dim(X.test)[1]
  tau = dim(X.test)[3]
  
  if(length(bInd) == 1){
    Xb = sapply(1:tau, function(t) X.test[,bInd,t]*n0b[t,bInd]/cEst[t])
  } else{
    Xb = sapply(1:tau, function(t) X.test[,bInd,t]%*%n0b[t,bInd]/cEst[t])
  }
  
  temp = matrix(nrow = m, ncol = tau)
  for(t in 1:tau) temp[,t] = rbinom(m, 1, pnorm(Xb[,t]))
  y.p = sapply(1:m, function(i){
    g = rle(temp[i,])
    g$values[which.max(g$lengths)]
  })
  
  y.p
}