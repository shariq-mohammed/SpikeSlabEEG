#' First-stage variable selection
#'
#' First-stage variable selection using the 'MCMC estimates of local Bayesian modeling 'Zcut', 'Wcut' or 'Ccut' method
#'
#' @param beta \eqn{tau} x \eqn{p} matrix of estimates of beta
#' @param method thresholding method to be used 'Zcut', 'Wcut' or 'Ccut'
#' @param w estimates of complexity parameters required if method == 'Wcut'
#' @param const threshold required if method == 'Ccut'
#' @keywords n0beta()
#' @export
#' @examples n0beta()

######################################
### First-stage variable selection ###
######################################
n0beta = function(beta,
                  method,
                  w=NULL,
                  const=NULL){
  if(method=="Zcut"){
    beta[abs(beta)<qnorm(0.975)]=0
    return(beta)
  }
  
  if(method=="Wcut"){
    if(is.null(w)) return("For Wcut method, w is required.")
    else{
      L = ncol(beta)
      u = sapply(1:nrow(beta),
                 function(t){
                   o = order(abs(beta[t,]), decreasing = TRUE)
                   beta[t, -o[1:floor(w[t]*L)]] = 0
                   beta[t,]
                 })
      return(t(u))
    }
  }
  
  if(method=="Ccut"){
    if(is.null(const)) return("For Ccut method, const is required.")
    else{
      beta[abs(beta)<const]=0
      return(beta)
    }
  }
}