#' @title linear transform of a vector or a matrix
#' @param x the vector/matrix to be transformed
#' @param min the minimum value after transformation
#' @param max the maximum value after transformation
#' @note if x contains negative value, the 
#' 
transform_linear <- function(x, min, max) {
  if (min < 0 || max < 0)
    stop("min and max should be greater than 0.")
  
  sg <- sign(x)
  x <- abs(x)
  x <- x - min(x)
  sc <- (max-min)/max(x)
  x <- x * sc + min
  x * sg
}