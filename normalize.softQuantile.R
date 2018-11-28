#' @param x the input matrix, the columns of which were centered on quantiles
#' @param probs the quantiles used to controls/align. If probs is a length 1 numerical 
#'   vector, the quantiles will aligned. If probs' length is > 1, the quantiles are calculated
#'   and linear models are used to estimated slope and intersect between quantiles.
#' @importFrom matrixStats colQuantiles

normalize.softQuantile <- function(x, probs = 0.5) {
  if (is.data.frame(x))
    x <- apply(x, 2, as.numeric)
  
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
  x
}
