#' @title basic QC for MQ output
#' @description basic QC, including barplot for IDs, boxplot and PCA
#' @param x the input matrix, usually a proteingroups table from maxquant
#' @param cols the columns to be included in the QC, usually the intensity, iBAQ or LFQ columns
#' @param pheno the pheno type vector, should be the same length as cols
#' @param xmpq whether x is an maxquant protein group input
#' @param log whether x should be log transformed
#' @import matrixStats
#' @import randomcoloR

basicMQC <- function(x, cols, pheno, xmpg = TRUE, log = TRUE) {
  
  require(matrixStats)
  require(randomcoloR)
  
  ord <- order(pheno)
  pheno <- pheno[ord]
  
  pheno <- as.factor(pheno)
  n <- nlevels(pheno)
  pal <- structure(distinctColorPalette(n), names = levels(pheno))
  
  if (xmpg) {
      cols <- cols[ord]
      i1 <- grepl("^CON", x$Majority.protein.IDs)
      i2 <- grepl("^REV", x$Majority.protein.IDs)
      i3 <- x$Only.identified.by.site == "+"
      i4 <- rowSums(emat, na.rm = TRUE) == 0
      i <- !(i1 | i2 | i3 | i4)
      emat <- x[, cols]
      colnames(emat) <- gsub("LFQ.intensity.|iBAQ.|Intensity.", "", colnames(emat))
      emat <- apply(emat[i, ], 2, as.numeric)
      emat[is.na(emat)] <- 0
  } else
      emat <- apply(x, 2, as.numeric)
  emat[is.infinite(emat)] <- NA
  # id plot
  layout(matrix(1:3, 3, 1))
  
  bp <- barplot(colSums(emat > 0, na.rm = TRUE), col = pal[as.character(pheno)], 
                las = 2, ylim = ceiling(c(0, nrow(emat)*1.05)), ylab = "# protein IDs")
  lines(bp, colSums(rowCumsums(emat) > 0), col = 1)
  points(bp, colSums(rowCumsums(emat) > 0), col = 1, pch = 19)
  legend("topleft", col = pal, pch = 15, legend = levels(pheno), bty = "n", pt.cex = 2)
  
  # boxplot
  if (log)
    logemat <- apply(emat, 2, log10) else
      logemat <- emat
  logemat[is.infinite(logemat)] <- NA
  boxplot(logemat, ylab = "Intensity (log10)", col = pal[as.character(pheno)], 
          las = 2)
  
  # pca 
  logemat[is.na(logemat)] <- min(logemat, na.rm = TRUE) - log10(2)
  pc <- prcomp(t(logemat))
  vars <- signif(pc$sdev^2/sum(pc$sdev^2), 3)
  plot(pc$x[, 1:2], col = pal[as.character(pheno)], pch = 19, cex = 2, 
       xlab = paste0("PC1 (", vars[1]*100, "%)"),
       ylab = paste0("PC2 (", vars[2]*100, "%)")
  )
  
}
