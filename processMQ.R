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
  pal <- structure(sort(distinctColorPalette(n)), names = levels(pheno))
  
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
  } else
      emat <- apply(x, 2, as.numeric)
  emat[is.infinite(emat)] <- NA
  
  layout(matrix(1:3, 3, 1))
  # id plot
  idmat <- apply(!is.na(emat), 2, as.integer)
  bp <- barplot(colSums(idmat, na.rm = TRUE), col = pal[as.character(pheno)], 
                las = 2, ylim = ceiling(c(0, nrow(emat)*1.05)), ylab = "# protein IDs")
  lines(bp, colSums(rowCumsums(idmat) > 0), col = 1)
  points(bp, colSums(rowCumsums(idmat) > 0), col = 1, pch = 19)
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


#' Read protein groups output of maxquant output and split
#'   it to columns
#' @param file Maxquant proteinGroup.txt file path

read.proteinGroups <- function(file) {
  pg <- read.delim(file)
  df <- data.frame(val = c("iBAQ.",
                           "LFQ.intensity.", 
                           "Peptides.", 
                           "Razor...unique.peptides.",
                           "Unique.peptides.",
                           "Sequence.coverage.",
                           "Intensity.",
                           "MS.MS.Count."), 
                   log = c(T, T, F, F, F, F, F, F))
  
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
  
  ml$iBAQ_mc <- sweep(ml$iBAQ, 2, colMedians(ml$iBAQ, na.rm = TRUE), "-") + median(ml$iBAQ, na.rm = TRUE)
  
  ml$annot <- annot
  ml
}


#' Do multiple t-test comparisons
#'   it to columns
#' @param file x expression matrix to be tested
#' @param label the label of columns
#' @compare the comparison to be done, a list of length-2 character vectors
#'   to indicate which groups should be compared

multi.t.test <- function(x, label, compare = NULL) {
  if (!is.list(compare))
    compare <- list(compare)
  
  lc <- lapply(compare, function(c1) {
    if (length(c1) == 2) {
      tv <- apply(x, 1, function(xx) {
        t <- try(t.test(xx[label == c1[1]], xx[label == c1[2]]), silent = TRUE)
        if (class(t) != "htest")
          return(c(md = NA, tstat = NA, pval = NA, df = NA))
        c(md = t$estimate[1] - t$estimate[2],
          tstat = t$statistic[["t"]], 
          pval = t$p.value,
          df = t$parameter[["df"]])
      })
      tv <- data.frame(
        md = as.numeric(tv[1, ]),
        df = as.numeric(tv[4, ]),
        tstat = as.numeric(tv[2, ]),
        pval = as.numeric(tv[3, ]),
        fdr = NA
      )
      tv$fdr[!is.na(tv$pval)] <- p.adjust(tv$pval[!is.na(tv$pval)], method = "fdr")
      
    } else {
      warning("Only allow comparing two groups, ignored.")
      tv <- NULL
    }
    tv
  })
  names(lc) <- sapply(compare, paste, collapse = "_")
  do.call(cbind, lc)
}

#' @title Basic QC for MQ output, used after calling "read.proteinGroups"
#' @description basic QC, including barplot for IDs, boxplot and PCA
#' @param x the input matrix, usually a proteingroups table from maxquant
#' @param group the group type vector, should be the same length as cols
#' @import matrixStats
#' @import randomcoloR

plotQC <- function(x, group) {
  
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
  bp <- barplot(colSums(idmat, na.rm = TRUE), col = pal[as.character(group)], 
                las = 2, ylim = ceiling(c(0, nrow(emat)*1.05)), ylab = "# protein IDs")
  lines(bp, colSums(rowCumsums(idmat) > 0), col = 1)
  points(bp, colSums(rowCumsums(idmat) > 0), col = 1, pch = 19)
  legend("topleft", col = pal, pch = 15, legend = levels(group), bty = "n", pt.cex = 2)
  
  # boxplot
  logemat <- emat
  logemat[is.infinite(logemat)] <- NA
  boxplot(logemat, ylab = "Intensity (log10)", col = pal[as.character(group)], 
          las = 2)
  
  # pca 
  logemat[is.na(logemat)] <- min(logemat, na.rm = TRUE) - log10(2)
  pc <- prcomp(t(logemat))
  vars <- signif(pc$sdev^2/sum(pc$sdev^2), 3)
  plot(pc$x[, 1:2], col = pal[as.character(group)], pch = 19, cex = 2, 
       xlab = paste0("PC1 (", vars[1]*100, "%)"),
       ylab = paste0("PC2 (", vars[2]*100, "%)")
  )
  
}




