# this file contains functions used in outlier analysis of an expression matrix
# To load functions in this file into your R working space, run:
# source("https://raw.githubusercontent.com/mengchen18/RFunctionCollection/master/outliers.R")


#' @title identify reliable outlier cut
#' @param x a log transformed matrix
#' @param nbreaks the number of value used to check in the range of row max of x
#' @param window the intensity window 
#' @param pvalue the p value threshold
#' @importForm matrixStats rowMaxs

outlierThresh <- function(x, nbreaks = 100, window = 0.05, pvalue = 0.05) {
  
  nna <- rowSums(is.na(x))
  nna <- as.factor(nna)
  rmax <- rowMaxs(x, na.rm = TRUE)
  rmax[is.infinite(rmax)] <- NA
  mi <- min(rmax, na.rm = TRUE)
  ma <- max(rmax, na.rm = TRUE)
  vec <- seq(mi, ma, length.out = nbreaks)
  
  max.na <- ncol(x) - 1
  rat <- sapply(vec, function(v) {
    i <- which(rmax < v + window & rmax > v - window)
    tb <- table(nna[i])
    round(tb[as.character(max.na)]/sum(tb), digits = 2)
  })
  
  le <- lowess(vec, rat, delta = 0.5, f = 1)
  cut <- vec[which(le$y < pvalue)[1]]
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))

  layout(matrix(1:3, 1, 3))
  plot(vec, rat, xlab = "intensity", ylab = "empirical p value (# of non-missing = 1)")
  abline(h = pvalue, v = cut, col = "green", lwd = 2)
  lines(le)
  
  hist(rmax, breaks = 50, main = 'Histogram of row max of x')
  abline(v = cut, col = "green", lwd = 2)
  
  hist(x, breaks = 50)
  abline(v = cut, col = "green", lwd = 2)
  
  cut
}

#' @title find outliers of every column
#' @param x a log transformed matrix
#' @param foldthresh
#' @param logbase the base of log transformation
#' @param reachLowBound logical. Whether the lowest intensity of a row should be lower than the threshould 
#'  determined by \code{outlierThresh}
#' @param ... other arguments passed to \code{outlierThresh}
#' @importFrom matrixStats rowMins

findOutlier <- function(x, foldthresh = 5, logbase = c(10, 2, "e")[1], reachLowBound = TRUE, ...) {
  
  bd <- switch(logbase, 
               "10" = log10(foldthresh), 
               "2" = log2(foldthresh), 
               "e" = log(foldthresh))
  
  sx <- apply(x, 1, sort, decreasing = TRUE, na.last = TRUE)
  v <- sx[1, ] - sx[2, ]
  i1 <- v > bd
  i2 <- !is.na(sx[1, ]) & is.na(sx[2, ])
  
  cc <- outlierThresh(x, ...)
  i3 <- sx[1, ] > cc  
  if (reachLowBound)
    i3[which(sx[2, ] > cc)] <- FALSE
  ii <- which(i3 & (i1 | i2))
  
  om <- x[ii, ]
  rmax <- rowMaxs(om, na.rm = TRUE)
  ct <- lapply(1:ncol(om), function(p) {
    ii[which(om[, p] == rmax)]
  })
  names(ct) <- colnames(x)
  list(outlierIndex = ii, outlierMatrix = om, outlierIndexColumns = ct)
}
