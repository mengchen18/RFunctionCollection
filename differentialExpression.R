require(limma)

#' row-wise pair wise t test
#' @param x the log10 transformed expression matrix where rows are proteins and columns are samples
#' @param mainFactor the factor of interest, it can only have two unique values for now
#' @param sideFactor other factors, such as batch factors, a vector or matrix where each column is a factor to be considered in the model
#' @param comparisons a list consists of vectors of lenght 2 to indicate which comparisons should be done. 
#'   If is NULL, all pairwise comparisons would be done. 
#' @param ... other parameters passed to limmaAddMod

pairwise.limmaAddMod <- function(x, mainFactor, sideFactor=NULL, comparisons=NULL, ...) {
  
  lv <- unique(mainFactor)
  if (length(lv) > 5)
    warning("More than 5 groups, might be too many comparisoins. Please consider use other method.")
  
  if (is.vector(sideFactor))
    sideFactor <- matrix(sideFactor, ncol = 1)
  
  if (is.null(comparisons)) {
    cbn <- combn(lv, 2)
    comparisons <- split(cbn, col(cbn))
  }
  
  ll <- lapply(comparisons, function(gp) {
    ii <- mainFactor %in% gp
    x <- x[, ii]
    mf <- mainFactor[ii]
    sf <- sideFactor[ii, ]
    limmaAddMod(x=x, mainFactor=mf, sideFactor=sf, ...)
  })
  
  names(ll) <- sapply(comparisons, paste, collapse="_")
  ll
}


#' @param x the log10 transformed expression matrix where rows are proteins and columns are samples
#' @param mainFactor the factor of interest, it can only have two unique values for now
#' @param sideFactor other factors, such as batch factors, a vector or matrix where each column is a factor to be considered in the model
#' @param impute logical; whether impute the matrix, half min of the row will be used.
#' @param imputeMin the minimum imputed value
#' @param nCompleteOneGroup integer; for a protien, the number of measured values at least in one group. If less than this number, the row will be ignored
#' @param nCompletePerGroup integer; for a protien, the number of measured values at least in all groups. If less than this number, the row will be ignored
#' @import limma

limmaAddMod <- function(x, mainFactor, sideFactor=NULL, impute=FALSE, imputeMin=NULL, nCompleteOneGroup=0, nCompletePerGroup=0) {
  
  # check
  # main facotr only allow two groups
  # dimension, etc
  
  # creating design matrix
  df <- data.frame(mainFactor, sideFactor)
  df[1:length(df)] <- lapply(df, as.factor)
  design <- model.matrix(~ ., df)
  
  # filter by nComplete
  xmat <- x
  nComplete <- data.frame(
    rowSums(!is.na(xmat[, df$mainFactor == levels(df$mainFactor)[1]])), 
    rowSums(!is.na(xmat[, df$mainFactor == levels(df$mainFactor)[2]]))
  )
  colnames(nComplete) <- levels(df$mainFactor)
  i <- (nComplete[[1]] >= nCompleteOneGroup | nComplete[[2]] >= nCompleteOneGroup) & 
    (nComplete[[1]] >= nCompletePerGroup & nComplete[[2]] >= nCompletePerGroup)
  excl <- rep(TRUE, nrow(xmat))
  excl[i] <- FALSE
  
  for (i in 1:nrow(xmat)) {
    if (excl[i]) {
      xmat[i, ] <- NA
    } else {
      # only do something if impute
      if (impute) {
        if (is.null(imputeMin))
          imputeMin <- min(x, na.rm = TRUE)
        ir <- is.na(xmat[i, ])
        xmat[i, ir] <- max(min(xmat[i, ], na.rm = TRUE) - log10(2), imputeMin)
      }
    }
  }
  
  # model fitting
  fit <- lmFit(xmat, design)
  fit <- eBayes(fit)
  tab <- topTable(fit, coef = paste0('mainFactor', levels(df$mainFactor)[2]), sort.by = "none", number = Inf)
  cbind(tab, n_measured = nComplete, excluded = excl, X=x, limmaX=xmat)
}
