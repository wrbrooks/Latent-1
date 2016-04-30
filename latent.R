#' Estimate parameters of a latent variable model
#' 
#' Estimate the parameters of a latent variable model. Estimation is by the
#' method of maximum likelihood, with conjugate gradient descent used to 
#' maximize likelihood. An EM algorithm is used to impute the values for
#' censored observations.
#' 
#' @param data matrix of observed counts, number of rows equal to the number of
#'     observations and columns equal to the number of categories.
#' @param min.detect minimum limit of detection for the FIB assay - values
#'     below this threshold are censored
#' @param event vector of event assignment labels, one entry per row of data.
#'     Indicates the event to which each row of data belongs.
#' @param specific vector of TRUE/FALSE values, one entry per column of data.
#'     Each entry indicates whether the corresponding FIB species is
#'     human-specific.
#' @param verbose should the function provide verbose output about its
#'     progress? Defaults to TRUE.
#' 
#' @return A list containing the following elements:
#'     \item{data}{final version of the data, where censored values have been
#'         imputed by the EM algorithm (see \code{\link{E.step}})}
#'     \item{min.detect}{a vector where each element is the minimum level of
#'         detection for the corresponding column of data (copied from function
#'         input)}
#'     \item{event}{vector listing the event ID associated with the
#'         corresponding row of the data (copied from the input)}
#'     \item{specific}{vector indicating whether the corresponding column of
#'         the data is a human-specific indicator bacteria - in which case its
#'         intercept is set to zero (copied from input)}
#'     \item{alpha}{matrix of estimated intercept parameters for each indicator
#'         bacteria for each event - each row is for one event, and each column
#'         is for one species of bacteria}
#'     \item{beta}{vector of estimated slope parameters for the indicator
#'         bacteria, representing the expected marginal increase in log-
#'         concentration when the contamination index increases by one unit}
#'     \item{gamma}{vector of estimated contamination indices, one for each row
#'         of the data}
#' 
#' @examples # This script is an example of using the latent package to 
#' # estimate the parameters of a latent variable model for the
#' # FIB counts in the storm sewer data set (with event 3 removed).
#' 
#' # Import data, drop event 3, change mei[4] from "TNTC" to 0, and then
#' # convert all FIB to numerics:
#' data("dfOptAnalysisDataSSJan2015")
#' indx = which(dfOptAnalysisDataSSJan2015$Event != "03")
#' fib = dfOptAnalysisDataSSJan2015[indx, c("mei", "modmtec", "FC",
#'     "Bac.human", "Lachno.2")]
#' fib$mei[4] = 0
#' for (n in names(fib))
#'     fib[[n]] = as.numeric(fib[[n]])
#' 
#' # Set the censoring values for each of the FIB (these are guesstimates)
#' min.detect = c('mei'=1, 'modmtec'=1, 'FC'=1, 'Bac.human'=225,
#'     'Lachno.2'=225)
#' 
#' # The human-specific FIB are Bac.human and Lachno.2, which are the fourth
#' # and fifth columns
#' specific = c(FALSE, FALSE, FALSE, TRUE, TRUE)
#' 
#' # Get the event IDs
#' event = as.integer(dfOptAnalysisDataSSJan2015$Event[indx])
#' 
#' # Now estimate the model parameters:
#' latent(fib, min.detect, event, specific)
#' 
#' @export
latent <- function(data, min.detect, event, specific=NULL, verbose=TRUE) {
  #Initial parameters:
  n = nrow(data)
  p = ncol(data)
  d = length(unique(event))
  #1 alpha value per column per event
  #2 beta value per column
  #Note that beta2 will have 0s for human specific FIB
  #Victor, 2 gammas
  #2 gamma value per row/reading
  #2 sigma values, 1 for each gamma
  #1 lambda per event for orthogonality constraint
#   xx = c(5.137,-2.658,1.265,0,0,4.648,-2.069,2.701,0,0,
#          4.507,-1.891,2.414,0,0,1.193,3.849,2.670,5.019,4.763,
#          2.413,3.837,3.054,0,0,
#          2.056,2.104,2.127,2.033,2.193,2.255,2.200,2.066,1.533,0.952,1.795,1.278,1.909,1.669,1.381,1.166,1.176,1.375,1.426,1.294,1.185,0.596,1.045,1.129,1.232,1.195,1.226,0.939,0.788,1.267,0.519,0.558,1.027,1.073,1.979,2.355,2.012,2.065,
#          1.214,0.950,1.075,0.208,1.101,1.045,1.127,0.111,0.056,-0.052,1.377,1.854,0.787,1.004,1.651,1.309,1.396,1.070,1.072,0.998,0.980,1.963,1.280,1.067,1.243,1.117,1.164,1.497,1.838,1.187,2.248,2.214,1.905,1.209,1.385,0.858,0.641,0.595,
#          1,1,1,1,1)
  cens <- sapply(1:p, function(k) ifelse(data[,k]<=min.detect[k], min.detect[k], data[,k]))
  xx = c(rep(as.integer(!specific), d), log(colMeans(cens)), as.integer(!specific), rep(1, n), rowMeans(sweep(log(cens), 2, log(colMeans(cens)), '-')), 1, 1, rep(15, d))
  finished = FALSE
  
  f.old = -Inf
  f.proposed = 0
  f.best = 0
  dir.new = vector()
  cdir.old = vector()
  dir.old = vector()
  tol = sqrt(.Machine$double.eps)
  tol = 1e-5
  check=Inf
  
  
  while (!finished) {
    #These iterations restart conjugacy:
    converged.cg = FALSE
    i = 0
    
    f.old = log.lik(data, xx, event)
    while (!converged.cg) {
      i = i+1
      t = 1
      
      dir.new = score(data, xx, event, specific = specific)
      
      dir.new = dir.new/sqrt(sum(dir.new^2))
      
      #Use conjugate gradient here
      if (i > 1) {
        cdir.old = cdir.old*max(0, sum(dir.new*(dir.new-dir.old))/(sum(dir.old^2)))
        dir.old = dir.new
        dir.new = dir.new + cdir.old
      } else {
        dir.old = dir.new
      }
      
      newxx = xx + dir.new * t
      
      foundStep = FALSE
      f.proposed = log.lik(data, newxx, event)
      
      #First make sure that we are at least improving
      while(!foundStep && !converged.cg){
        if(f.proposed > f.old)
          foundStep = TRUE
        else {
          t <- t*0.9
          if(t <= .Machine$double.eps){
            #Step size too small, just take original parameters
            converged.cg = TRUE
            newxx = xx
            cat("Step size too small, converged") 
          } else {
            newxx <- xx + dir.new * t
            f.proposed = log.lik(data, newxx, event)
          }
        }
      }
      #Now we try to find the optimal step
      foundStep = FALSE
      f.best = f.proposed
      while(!foundStep && !converged.cg){
        
        t <- t * 0.9
        if(t <= .Machine$double.eps){
          converged.cg = TRUE
          cat("Step size too small, converged") 
        }
        newxx <- xx + dir.new * t
        
        f.proposed <- log.lik(data, newxx, event)
        
        if(f.proposed > f.best){
          f.best = f.proposed
        }
        else{
          t <- t/0.9
          newxx <- xx + dir.new * t
          foundStep <- TRUE
        }
      }
      
      xx = newxx

      
      f.proposed = log.lik(data, xx, event)
      
      
      if (i%%10 == 0) {
        cat(paste(i, " iterations of CG, step size ", t, "likelihood at ",f.proposed ,"\n"))
        if (i%%1000 == 0) {
          print.table(xx)
          if (i%%10000 == 0) {
            cat("Taken the maximum amount of steps, treat as converged")
            converged.cg = TRUE
          }
        }
      }
      
      if ((f.proposed - f.old) < f.old * tol){ 
        converged.cg = TRUE
        cat(paste("Converged, new step: ", f.proposed, " old step = ", f.old, "\n"))
      }
      
      f.old = f.proposed
      cdir.old = dir.new
    }
    
    converged.irls = FALSE
    gamma.old <- xx[(2+d)*p + 1:(2*n)]
    f.old = log.lik(data, xx, event)
    i = 0
    while (!converged.irls) {
      i <- i+1
      #Do IRLS here
      gamma.new = irls(data, xx, event)
      
      xx[(2+d)*p + 1:(2*n)] <- gamma.new
      
      if (i%%10 == 0) {
        cat(paste(i, " iterations of IRLS, likelihood = ", log.lik(data, xx, event), "\n"))
      }
      
      if (sum((gamma.old - gamma.new)^2 / sum(gamma.old^2)) < tol)
        converged.irls <- TRUE
      else{
        gamma.old <- gamma.new
      }
    }
    
    #2 betas per col
    #2 gammas per row
    #2 mus/sigmas per event
    alpha <- xx[1:(p * d)]
    beta <- xx[(p * d + 1):(p * (2 + d))]
    gamma <- xx[(p * (2 + d) + 1):(p * (2 + d) + 2 * n)]
    sigma1 <- xx[p * (2 + d) + 2 * n + 1]
    sigma2 <- xx[p * (2 + d) + 2 * n + 2]
    lambda <- tail(xx, d)
    
    data.new = E.step(alpha, beta, gamma, data, min.detect, event)
    
    check.old = check
    check = sum((data.new - data)^2, na.rm=TRUE) / sum(data^2, na.rm=TRUE)
    data = data.new
    
    cat(paste("E-Step, check: ", check.old, ", ", check, "\n"))
    if (check < tol) {
      if (tol<=1e-7)
        finished = TRUE
      tol = tol/2
      cat(paste("Iterating with tol=", tol, "\n", sep=""))
      
    }
  }
  
  cat(log.lik(data, xx, event))
  #Compile the results and return
  alpha <- xx[1:(p * d)]
  beta <- xx[(p * d + 1):(p * (2 + d))]
  gamma <- xx[(p * (2 + d) + 1):(p * (2 + d) + 2 * n)]
  sigma1 <- xx[p * (2 + d) + 2 * n + 1]
  sigma2 <- xx[p * (2 + d) + 2 * n + 2]
  lambda <- tail(xx, d)
  beta1 <- beta[1:p]
  beta2 <- tail(beta, p)
  gamma1 <- gamma[1:n]
  gamma2 <- tail(gamma, n)
  
  result = list()
  result$data = data
  result$min.detect = min.detect
  result$event = event
  result$specific = specific
  result$alpha = matrix(alpha, d, p, byrow=TRUE)
  colnames(result$alpha) = names(result$beta)
  result$beta1 = beta1
  result$gamma1 = gamma1
  result$sigma1 = sigma1
  result$beta2 = beta2
  result$gamma2 = gamma2
  result$sigma2 = sigma2
  result$lambda = lambda
  
  options(scipen=999)
  return(result)
}