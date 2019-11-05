#' @param x the input matrix, the columns of which were centered on quantiles
#' @param probs the quantiles used to controls/align. If probs is a length 1 numerical 
#'   vector, the quantiles will aligned. If probs' length is > 1, the quantiles are calculated
#'   and linear models are used to estimated slope and intersect between quantiles.
#' @param sharedProtein normalization based on the shared proteins between reference 
#'   experiment each individual experiment
#' @param ref the columns name or index to specify the reference experiment, only used 
#'   when sharedProtein = TRUE
#' @importFrom matrixStats colQuantiles

require(matrixStats)

normalize.softQuantile <- function(x, probs = 0.5, sharedProtein = FALSE, ref = 1) {
  
  if (is.data.frame(x))
    x <- apply(x, 2, as.numeric)
  
  if (sharedProtein) {
    if (length(probs) == 1) {
      fac <- sapply(1:ncol(x), function(i) {
        ir <- !is.na(x[, ref]) & !is.na(x[, i])
        quantile(x[ir, ref], probs = probs) - quantile(x[ir, i], probs = probs)
      })
      x <- sweep(x, 2, fac, "+")
    } else {
      for (i in 1:ncol(x)) {
        ir <- !is.na(x[, ref]) & !is.na(x[, i])
        refquant <- quantile(x[ir, ref], probs = probs, na.rm = TRUE)
        colquant <- quantile(x[ir, i], probs = probs, na.rm = TRUE)
        mod <- lm(refquant ~ colquant)
        x[, i] <- mod$coefficients[[2]] * x[, i] + mod$coefficients[[1]]
      }
    }
  } else {
    if (length(probs) == 1) {
      colquant <- colQuantiles(x, probs = probs, na.rm = TRUE)
      grandquant <- quantile(x, probs = probs, na.rm = TRUE)
      off <- colquant - grandquant
      x <- sweep(x, 2, off, "-")
    } else {
      colquant <- t(colQuantiles(x, probs = probs, na.rm = TRUE))
      grandquant <- quantile(x, probs = probs, na.rm = TRUE)
      for (i in 1:ncol(x)) {
        mod <- lm(grandquant ~ colquant[, i])
        x[, i] <- mod$coefficients[[2]] * x[, i] + mod$coefficients[[1]]
      }
    }
  }
  x
}
