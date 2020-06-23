plotQC <- function(x, group, labelPCA = FALSE, fillNA = TRUE, checkPC = 1:2, type = c("protein", "metabolites")[1]) {
  
  require(matrixStats)
  require(randomcoloR)
  
  ord <- order(group)
  group <- group[ord]
  
  group <- as.factor(group)
  n <- nlevels(group)
  pal <- structure(sort(distinctColorPalette(n)), names = levels(group))
  emat <- x[, ord]
  
  layout(matrix(1:3, 3, 1))
  # id plot
  idmat <- apply(!is.na(emat), 2, as.integer)
  idna <- idmat
  idna[idna == 0] <- NA
  shareid <- colSums(!is.na(rowCumsums(idna)))
  
  bp <- barplot(colSums(idmat, na.rm = TRUE), col = pal[as.character(group)], 
                las = 2, ylim = ceiling(c(0, nrow(emat)*1.05)), ylab = sprintf("# %s IDs", type))
  lines(bp, colSums(rowCumsums(idmat) > 0), col = 1)
  points(bp, colSums(rowCumsums(idmat) > 0), col = 1, pch = 19)
  lines(bp, shareid, col = 1, lty = 2)
  points(bp, shareid, col = 1, pch = 15)
  legend("topleft", col = pal, pch = 15, legend = levels(group), bty = "n", pt.cex = 2)
  
  # boxplot
  logemat <- emat
  logemat[is.infinite(logemat)] <- NA
  boxplot(logemat, ylab = "Intensity (log10)", col = pal[as.character(group)], 
          las = 2)
  
  # pca 
  if (fillNA) 
    logemat[is.na(logemat)] <- min(logemat, na.rm = TRUE) - log10(2) else
      logemat <- logemat[!is.na(rowSums(logemat)), ]
  
  pc <- prcomp(t(logemat))
  vars <- signif(pc$sdev^2/sum(pc$sdev^2), 3)
  plot(pc$x[, checkPC], col = pal[as.character(group)], pch = 19, cex = 2, 
       xlab = paste0("PC", checkPC[1], " (", vars[checkPC[1]]*100, "%)"),
       ylab = paste0("PC", checkPC[2], " (", vars[checkPC[2]]*100, "%)")
  )
  if (labelPCA)
    text(pc$x[, checkPC], labels = as.character(group))
  
  list(ID = shareid, 
       pc = pc$x,
       pc_sdev = pc$sdev
      )
}
