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
  lsind$annot <- ab[-i, ]
  lsind$Summed <- ls
  
  lsind
}

#' Do multiple t-test comparisons
#'   it to columns
#' @param file x expression matrix to be tested
#' @param label the label of columns
#' @param compare the comparison to be done, a list of length-2 character vectors
#'   to indicate which groups should be compared
#' @param ... other parameters passed to the t.test function

#' Do multiple t-test comparisons
#'   it to columns
#' @param file x expression matrix to be tested
#' @param label the label of columns
#' @param compare the comparison to be done, a list of length-2 character vectors
#'   to indicate which groups should be compared
#' @param xlsx.file the xlsx output file
#' @param other.sheets other tables that should be write into the xlsx.file, it has to
#'   be a named list. Only used when xlsx.file is not NULL.
#' @param ... other parameters passed to the t.test function
#' @import openxlsx

multi.t.test <- function(x, label, compare = NULL, xlsx.file = NULL, other.sheets = list(), ...) {
  if (!is.list(compare))
    compare <- list(compare)
  
  lc <- lapply(compare, function(c1) {
    if (length(c1) == 2) {
      
      m1 <- rowMeans(x[, label == c1[1]], na.rm = TRUE)
      m2 <- rowMeans(x[, label == c1[2]], na.rm = TRUE)
      df <- data.frame(m1 = m1, m2 = m2, q1 = NA, q2 = NA)
      df$n1 <- rowSums(!is.na(x[, label == c1[1]]))
      df$n2 <- rowSums(!is.na(x[, label == c1[2]]))
      nom1 <- na.omit(df$m1)
      nom2 <- na.omit(df$m2)
      df$q1[!is.na(df$m1)] <- rank(nom1)/length(nom1)
      df$q2[!is.na(df$m2)] <- rank(nom2)/length(nom2)
      df$md <- m1 - m2
      
      tv <- apply(x, 1, function(xx) {
        t <- try(t.test(xx[label == c1[1]], xx[label == c1[2]], ...), silent = TRUE)
        if (class(t) != "htest")
          return(c(tstat = NA, pval = NA, df = NA))
        
        c(tstat = t$statistic[["t"]], 
          pval = t$p.value,
          df = t$parameter[["df"]])
      })
      tv <- data.frame(
        df = as.numeric(tv[3, ]),
        tstat = as.numeric(tv[1, ]),
        pval = as.numeric(tv[2, ]),
        fdr = NA
      )
      tv <- cbind(df, tv)
      tv$fdr[!is.na(tv$pval)] <- p.adjust(tv$pval[!is.na(tv$pval)], method = "fdr")
      
    } else {
      warning("Only allow comparing two groups, ignored.")
      tv <- NULL
    }
    tv
  })
  names(lc) <- sapply(compare, paste, collapse = "_")
  
  if (!is.null(xlsx.file)) {
    wb <- createWorkbook("BayBioMS")
    
    addWorksheet(wb, "Groups")
    writeData(wb, sheet = "Groups", x = cbind(sample = colnames(x), group = label))
    
    x1 <- rbind(
      c("Description of headers in the t_[xxx] sheets (results of t-tests)", ""),
      c("Header", "Description"),
      c("m_[xxx]", "the average intensity of the group [xxx]."),
      c("q_[xxx]", "the rank of average intensity in the group [xxx], it is normalized to a range from 0 (low intensity) to 1 (high intensity)."),
      c("n_[xxx]", "the number of detected values in the group [xxx]."),
      c("md", "Mean difference"),
      c("tstat", "The t-statistics of t-test"),
      c("df", "The degree of freedom of t-test"),
      c("pval", "The p-value returned by t-test"),
      c("fdr", "Benjamini and Hochberg corrected p-value, which controls the false discovery rate (FDR)")
    )
    
    addWorksheet(wb, "ttest_description")
    writeData(wb, sheet = "ttest_description", x = x1, colNames = FALSE)
    
    for (i in 1:length(lc)) {
      dt <- lc[[i]]
      colnames(dt) <- gsub("1$", paste0("_", compare[[i]][1]), colnames(dt))
      colnames(dt) <- gsub("2$", paste0("_", compare[[i]][2]), colnames(dt))
      addWorksheet(wb, paste0("t_", names(lc)[i]))
      writeData(wb, sheet = paste0("t_", names(lc)[i]), x = dt)
    }
    
    if (length(other.sheets) > 0) {
      if (is.null(names(other.sheets)))
        stop("other.sheets has to be named list!")
      for (i in names(other.sheets)) {
        addWorksheet(wb, i)
        writeData(wb, sheet = i, x = other.sheets[[i]])
      }
    }
    
    addWorksheet(wb, "Intensity")
    writeData(wb, sheet = "Intensity", x = x)
    saveWorkbook(wb, xlsx.file, overwrite = TRUE)
  }
  lc
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
  logemat[is.na(logemat)] <- min(logemat, na.rm = TRUE) - log10(2)
  pc <- prcomp(t(logemat))
  vars <- signif(pc$sdev^2/sum(pc$sdev^2), 3)
  plot(pc$x[, 1:2], col = pal[as.character(group)], pch = 19, cex = 2, 
       xlab = paste0("PC1 (", vars[1]*100, "%)"),
       ylab = paste0("PC2 (", vars[2]*100, "%)")
  )
  
}




