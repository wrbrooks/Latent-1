irls <- function(data, params, event) {
  #Basic constants:
  n <- nrow(data)
  p <- ncol(data)
  d <- length(unique(event))
  
  # Extract parameters from the parameter vector.
  # sigma is length 2
  # lambda is the lagrange multiplier
  alpha <- params[1:(p*d)]
  beta <- params[p*d + 1:(2*p)]
  gamma <- params[(2+d)*p + 1:(2*n)]
  sigma1 <- params[(2+d)*p + 2*n + 1]
  sigma2 <- params[(2+d)*p + 2*n + 2]
  lambda <- tail(params, d)
  
  # read the separate elements of beta and gamma
  # beta is 2 per column
  # gamma is 2 per row
  beta1 <- beta[1:p]
  beta2 <- tail(beta, p)
  gamma1 <- gamma[1:n]
  gamma2 <- tail(gamma, n)
  
  #Make structure for alpha:
  alpha.local = matrix(0, n, p)
  
  #Make the lagrange multipliers here
  lambda.local = vector()
  for (k in 1:d) {
    indx <- which(event==unique(event)[k])
    alpha.local[indx,] = matrix(alpha[(p*(k-1)+1):(p*k)], length(indx), p, byrow=TRUE)
    lambda.local <- c(lambda.local, lambda[k]*gamma2[indx])
  }
  
  #Add zeros to fit column width for the design matrix
  lambda.local = c(lambda.local, rep(0,n))
  
  #Make eta into vector, should be 1x(5*38) = 1x190
  eta = as.vector(gamma1 %*% t(beta1) + gamma2 %*% t(beta2) + alpha.local)
  mu.local = exp(eta)
  wt = c(mu.local,rep(1,n),rep(1,n),1)
  
  #Deduct alpha.local here as the design matrix will not include it
  z = eta - as.vector(alpha.local) + (as.vector(as.matrix(data)) - mu.local) / mu.local
  
  #Create response vector, consisting of vector z, 2 vectors of n length for the gaussian priors, and 1 entry for the lagrange multiplier
  response = c(z, rep(0, n), rep(0, n), 0)
  
  
  #Make design matrix
  #Block matrices of [(beta1, beta2),(sigma1, 0),(0, sigma2)]
  designMat = rbind(cbind(diag(n)*beta1[1], diag(n)*beta2[1]), cbind(diag(n)*beta1[2], diag(n)*beta2[2]),
                    cbind(diag(n)*beta1[3], diag(n)*beta2[3]), cbind(diag(n)*beta1[4], diag(n)*beta2[4]),
                    cbind(diag(n)*beta1[5], diag(n)*beta2[5]), cbind(diag(n)/sqrt(2*sigma1[1]), matrix(0,n,n)),
                    cbind(matrix(0, n, n), diag(n)/sqrt(2*sigma2[1])), lambda.local)
  
  pois = lsfit(x=designMat, y=response, intercept=FALSE, wt=wt)
  
  gamma.new = pois$coefficients
  
  gamma.new
}