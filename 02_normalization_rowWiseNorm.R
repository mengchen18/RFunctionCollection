#' @param x a matrix where rows are genes/proteins and columns are samples
#' @param batch a factor or vecter has the same length as ncol(x) to indicate
#'   the batch assignment of samples
#' @param ref a logical vector has the same length as ncol(x) to indicated which
#'   columns are the common references among batches. If it is NULL (by default), 
#'   the mean of all channels will be used as batch reference.
#' @importFrom matrixStats rowMeans

rowshift <- function(x, batch, ref=NULL) {
  
  if (is.data.frame(x))
    x <- apply(x, 2, as.numeric)
  if (is.null(ref))
    ref <- 1:ncol(x)
  
  b_ref <- batch[ref]
  expr_ref <- x[, ref, drop=FALSE]
  grandmeans <- rowMeans(expr_ref, na.rm = TRUE)
  for (i in unique(batch)) {
    off <- rowMeans(expr_ref[, b_ref == i, drop = FALSE], na.rm = TRUE) - grandmeans
    x[, batch == i] <- x[, batch == i] - off
  }
  x
  
}
