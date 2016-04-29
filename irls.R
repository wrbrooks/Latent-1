irls <- function(data, params, event) {
  #Basic constants:
  n = nrow(data)
  p = ncol(data)
  d = length(unique(event))
  
  #Extract parameters from the parameter vector:
  alpha = params[1:(p*d)]
  beta = params[(p*d+1):(p*d+p)]
  gamma = params[(p*d+1+p):(p*d+p+n)]
  sigma = params[(p*d+1+p+n)]
  
  #Event-level alphas:
  alpha.local = matrix(0, n, p)
  
  sumLogLikNormal = 0
  for (k in 1:d) {
    indx = which(event==unique(event)[k])
    alpha.local[indx,] = matrix(alpha[(p*(k-1)+1):(p*k)], length(indx), p, byrow=TRUE)
    
  }

  eta = as.vector(alpha.local + gamma%*%t(beta))
  mu.local = exp(eta)
  
  z = eta - as.vector(alpha.local) + (as.vector(as.matrix(data))-mu.local)*1/mu.local
  response = c(z, rep(0,n))
  wt = c(mu.local, rep(1,n))
  
  
  designMat = rbind(diag(n)*beta[1],diag(n)*beta[2],diag(n)*beta[3],diag(n)*beta[4],diag(n)*beta[5],diag(n)*1/sqrt(2*sigma))

  pois = lsfit(designMat, response, wt=wt, intercept=FALSE)
  gamma.new = pois$coefficients
  
  
  gamma.new
}