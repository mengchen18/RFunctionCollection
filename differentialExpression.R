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