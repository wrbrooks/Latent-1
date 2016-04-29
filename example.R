#' An  example illustrating the function of packace latent
#' 
#' Estimates a latent contamination variable for the storm sewer data, using EM to impute censored values
#' 
#' @return latent a list as returned by function \code{\link{latent}}
#' 
#' @export
example <- function() {
  #---------
  # This script is an example of using the latent package to 
  # estimate the parameters of a latent variable model for the
  # FIB counts in the storm sewer data set (with event 3 removed).
  # Note that the data is automatically lazy-loaded, so no call
  # to data(...) is necessary.
  
  # Drop event 3, change mei[4] from "TNTC" to 0, and then convert all FIB to numerics:
  indx = which(dfOptAnalysisDataSSJan2015$Event != "03")
  fib = dfOptAnalysisDataSSJan2015[indx, c("mei", "modmtec", "FC", "Bac.human", "Lachno.2")]
  fib$mei[4] = 0
  for (n in names(fib))
    fib[[n]] = as.numeric(fib[[n]])
  
  # Set the censoring values for each of the FIB (these are guesstimates)
  min.detect = c('mei'=1, 'modmtec'=1, 'FC'=1, 'Bac.human'=225, 'Lachno.2'=225)
  
  # The human-specific FIB are Bac.human and Lachno.2, which are the fourth and fifth columns
  specific = c(FALSE, FALSE, FALSE, TRUE, TRUE)
  
  # Get the event IDs
  event = as.integer(dfOptAnalysisDataSSJan2015$Event[indx])
  
  # Now estimate the model parameters:
  latent(fib, min.detect, event, specific)
}