#' Read protein groups output of maxquant output and split
#'   it to columns
#' @param file Maxquant proteinGroup.txt file path

read.proteinGroups <- function(file) {
  pg <- read.delim(file, stringsAsFactors = FALSE)
  df <- data.frame(val = c("iBAQ.",
                           "LFQ.intensity.", 
                           "Peptides.", 
                           "Razor...unique.peptides.",
                           "Unique.peptides.",
                           "Sequence.coverage.",
                           "Intensity.",
                           "MS.MS.Count.",
                           "MS.MS.count."
  ), 
  log = c(T, T, F, F, F, F, F, F, F), 
  stringsAsFactors = FALSE)
  
  vi <- sapply(df$val, function(x) length(grep(x, colnames(pg))) > 0)
  df <- df[vi, ]             
  
  i <- ! (grepl("^REV_", pg$Majority.protein.IDs) |
            grepl("^CON_", pg$Majority.protein.IDs) |
            pg$Only.identified.by.site == "+")
  
  annot <- pg[i, -grep(paste(df$val, collapse = "|"), colnames(pg))]
  
  getExpr <- function(x, type = "iBAQ.", log = TRUE, keep.row = NULL) {
    val <- apply(pg[, grep(type, colnames(x), ignore.case = TRUE)], 2, as.numeric)
    if (log) {
      val <- log10(val)
      val[is.infinite(val)] <- NA
    }
    if (!is.null(keep.row))
      val <- val[keep.row, ]
    colnames(val) <- gsub(type, "", colnames(val))
    val
  }
  
  ml <- mapply(function(val, log) getExpr(pg, type = val, log = log, keep.row = i), 
               val = df$val, log = df$log)
  names(ml) <- gsub("\\.$", "", df$val)
  
  if (!is.null(ml$iBAQ))
    ml$iBAQ_mc <- sweep(ml$iBAQ, 2, colMedians(ml$iBAQ, na.rm = TRUE), "-") + median(ml$iBAQ, na.rm = TRUE)
  
  ml$annot <- annot
  ml
}

read.proteinGroups.tmt <- function(file) {
  ab <- read.delim(file, stringsAsFactors = FALSE)
  
  ir <- c(grep("^CON_", ab$Majority.protein.IDs), 
         grep("^REV_", ab$Majority.protein.IDs), 
         which(ab$Only.identified.by.site == "+"))
  
  eSum <- c("Fraction", "Reporter.intensity.corrected", "Reporter.intensity", "Reporter.intensity.count")
  ls <- list()
  for (i in eSum) {
    gb <- grep(paste0(i, ".[0-9]*$"), colnames(ab), value = TRUE)
    ls[[i]] <- apply(ab[-ir, gb], 2, as.numeric)
    ab[gb] <- NULL
  }
  
  lsind <- list()
  eInd <- c("Reporter.intensity.corrected", "Reporter.intensity.count", "Reporter.intensity")
  for (i in eInd) {
    gb <- grep(i, colnames(ab), value = TRUE)
    lsind[[i]] <- apply(ab[-ir, gb], 2, as.numeric)
    ab[gb] <- NULL
  }
  
  lsind$Reporter.intensity.corrected.log2 <- log2(lsind$Reporter.intensity.corrected)
  lsind$Reporter.intensity.corrected.log2[is.infinite(lsind$Reporter.intensity.corrected.log2)] <- NA
  lsind$annot <- ab[-ir, ]
  lsind$Summed <- ls
  
  lsind
}



plotQC <- function(x, group, labelPCA = FALSE, fillNA = TRUE, checkPC = 1:2) {
  
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
                las = 2, ylim = ceiling(c(0, nrow(emat)*1.05)), ylab = "# protein IDs")
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
}
