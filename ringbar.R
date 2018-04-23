#' @title barplot of matrix in ring
#' @param x a matrix, top row is for outer ring
#' @param col a matrix has the same dimension as \code{x} to indicate the color of bars
#' @param gap.degree degree of gap
#' @param start.degree start degree of rings
#' @param y.axis logical; whether the y axis should be drawn
#' @param x.axis.inner logical; whether the x axis of the inner ring should be drawn
#' @param x.axis logical; whether the x axis of the rings should be drawn
#' @param track.margin.exp the expansion factor of track.margin
#' @param track.height passed to \code{circos.trackPlotRegion}
#' @import circlize
#' @examples 
#'  x <- matrix(rnorm(600), nrow = 6)
#'  col <- row(x)
#'  ringbar(x)

ringbar <- function(x, col = 1, gap.degree = 40, start.degree = 90, track.margin.exp = 1, 
  y.axis = TRUE, x.axis.inner = TRUE, x.axis = FALSE, track.height = 0.1) {

  n <- ncol(x)
  if (length(col) == 1)
    col <- matrix(col, nrow(x), ncol(x))

  if (any(dim(x) != dim(col)))
    stop("The dimension of color should be the same as x.")
  
  circos.clear()
  circos.par(gap.degree = gap.degree, start.degree = start.degree)
  circos.initialize(factors = factor(rep("x", n+1)), x = 0:n) #
  
  for (i in 1:nrow(x)) {
    x1 <- x[i, ]
    col1 <- col[i,]
    
    bot <- 0.95 * min(x1, na.rm = TRUE)
    top <- 1.01 * max(x1, na.rm = TRUE)

    circos.trackPlotRegion(ylim = c(bot, top), 
                           track.margin = c(0.01, 0.01)*track.margin.exp,
                           cell.padding = c(0,0,0,0), 
                           track.height = track.height, 
                           bg.border  = FALSE)
    circos.rect(1:n-1, ybottom = rep(bot, n), 1:n, x1, col = col1, border = "white")
    if (y.axis)
      circos.yaxis(side = "left", labels.niceFacing = T, 
        tick.length = convert_x(0.1, "mm", get.cell.meta.data("sector.index"), get.cell.meta.data("track.index")), 
        labels.cex = 0.3, at = round(max(x1, na.rm = TRUE), digits = 2))
    
    if (x.axis)
      circos.axis(h = "bottom", major.at = 1:n-0.5, labels = "", labels.cex = 0.5,
                  minor.ticks = 0, major.tick = FALSE, direction = "inside", col = "gray",
                  labels.facing = "reverse.clockwise")
    
  }
  if (x.axis.inner)
    circos.axis(h = "bottom", major.at = 1:n-0.5, labels = "", labels.cex = 0.5,
                minor.ticks = 0, major.tick = FALSE, direction = "inside", 
                labels.facing = "reverse.clockwise")
  
}

