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
                           "MS.MS.count.",
                           "Identification.type."
  ), 
  log = c(T, T, F, F, F, F, F, F, F, F), 
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
    ml$iBAQ_mc <- sweep(ml$iBAQ, 2, matrixStats::colMedians(ml$iBAQ, na.rm = TRUE), "-") + median(ml$iBAQ, na.rm = TRUE)
  
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

#' Read modificationSpecificPeptides output of maxquant output and split
#'   it to columns
#' @param file Maxquant proteinGroup.txt file path
#' @param modifications the modification to be included
read.modSepPep <- function(file, modifications = "Phospho \\(STY\\)") {
  
  dat <- read.delim(file, stringsAsFactors = FALSE)
  df <- data.frame(val = c("Intensity.", "Fraction.", "Experiment."), 
                   log = c(T, F, F), 
                   stringsAsFactors = FALSE)
  
  vi <- sapply(df$val, function(x) length(grep(x, colnames(dat))) > 0)
  df <- df[vi, ]             
  i <- dat$Potential.contaminant != "+" & dat$Reverse != "+" & grepl("Phospho \\(STY\\)", dat$Modifications)
  
  annot <- dat[i, -grep(paste(df$val, collapse = "|"), colnames(dat))]
  nn <- grep("Intensity.", colnames(dat))
  getExpr <- function(x, type = "iBAQ.", log = TRUE, keep.row = NULL) {
    if (length(grep(type, colnames(x), ignore.case = TRUE)) < length(nn))
      return(NULL)
    val <- apply(dat[, grep(type, colnames(x), ignore.case = TRUE)], 2, as.numeric)
    if (log) {
      val <- log10(val)
      val[is.infinite(val)] <- NA
    }
    if (!is.null(keep.row))
      val <- val[keep.row, ]
    colnames(val) <- gsub(type, "", colnames(val))
    val
  }
  ml <- mapply(function(val, log) getExpr(dat, type = val, log = log, keep.row = i), 
               val = df$val, log = df$log)
  names(ml) <- gsub("\\.$", "", df$val)
  
  ml$annot <- annot 
  ml
}

