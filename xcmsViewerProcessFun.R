

multi.t.test2 <- function(x, pheno, compare = NULL, log10 = FALSE, median.center = TRUE, fillNA = TRUE, ...) {
  
  if (is.null(rownames(x)))
    rownames(x) <- 1:nrow(x)
  
  halfValue <- function(x) x/2
  if (log10) {
    x <- apply(x, 2, log10)
    halfValue <- function(x) x - log10(2)
  }
  
  tl <- lapply(unique(compare[, 1]), function(x) {
    x <- compare[x == compare[, 1], -1]
    unique(unlist(split(x, row(x))))
  })
  names(tl) <- unique(compare[, 1])
  
  if (median.center)
    x <- sweep(x, 2, matrixStats::colMedians(x, na.rm = TRUE), "-") + median(x, na.rm = TRUE)
  
  df <- data.frame(metabolite = rownames(x), stringsAsFactors = FALSE)
  for ( i in names(tl) ) {
    for ( j in tl[[i]] ) {
      m <- x[, pheno[[i]] == j]
      rv <- rowSums(!is.na(m))
      rm <- rowMeans(m, na.rm = TRUE)
      rq <- rank(rm, na.last = TRUE)/sum(!is.na(rm))
      rq[is.na(rm)] <- NA
      df[[paste("mean", j)]] <- rm
      df[[paste("n value", j)]] <- rv
      df[[paste("quantile", j)]] <- rq
    }
  }
  
  if (fillNA) {
    x <- apply(x, 1, function(xx) {
      x3 <- xx
      x3[is.na(x3)] <- halfValue(min(x3, na.rm = TRUE))
      x3
    })
    x <- t(x)
  }
  
  for ( i in 1:nrow(compare) ) {
    v <- compare[i, ]
    i1 <- pheno[[v[1]]] == v[2]
    i2 <- pheno[[v[1]]] == v[3]
    
    tv <- apply(x, 1, function(xx) {
      t <- try(t.test(xx[i1], xx[i2], ...), silent = TRUE)
      if (class(t) != "htest")
        return(c(pvalue = NA, mean.diff = NA))
      c(pvalue = t$p.value, mean.diff = t$estimate[1] - t$estimate[2])
    })
    
    df[[paste("ttest", "pvalue", v[2], v[3], sep = "|")]] <- tv[1, ]
    df[[paste("ttest", "log.pvalue", v[2], v[3], sep = "|")]] <- -log10(tv[1, ])
    df[[paste("ttest", "fdr", v[2], v[3], sep = "|")]] <- p.adjust(tv[1, ], method = "fdr")
    df[[paste("ttest", "log.fdr", v[2], v[3], sep = "|")]] <- -log10(df[[paste("ttest", "fdr", v[2], v[3], sep = "|")]])
    df[[paste("ttest", "mean.diff", v[2], v[3], sep = "|")]] <- tv[2, ]
  }
  
  if (log10)
    x.ori <- 10^x else
      x.ori <- x
      
    
  list(ttest = df, mat = x, mat.rawscale = x.ori)
}

#####
pca <- function(x, n = 6) {
  
  writePC <- function(x, n) {
    var <- round(x$sdev[1:n]^2/(sum(x$sdev^2))*100, digits = 1)
    xx <- x$x[, 1:min(n, ncol(x$x))]
    colnames(xx) <- paste0(colnames(xx), " (", var, "%", ")")
    pp <- x$rotation[, 1:min(n, ncol(x$x))]
    colnames(pp) <- paste0(colnames(pp), " (", var, "%", ")")
    list(samples = xx, features = pp)
  }
  
  pc <- prcomp(t(x))
  writePC(pc, n = n)
}

###
phenoFeatureData <- function(
  object, compare, pheno=NULL, log10 = TRUE, median.center = TRUE, fillNA = TRUE, ...
) {
  
  if (is.null(pheno))
    pheno <- Biobase::pData(object@featureSet)
  
  ts <- multi.t.test2(
    x = Biobase::exprs(object@featureSet), 
    compare = compare, 
    pheno = pheno,
    log10 = log10, 
    median.center = median.center, 
    fillNA = fillNA, ...)
  
  mat <- ts$mat
  mat.ori <- ts$mat.rawscale
  ts <- ts$ttest
  pc <- pca(mat)
  
  fd <- cbind(fData(object@featureSet), ts[-1], pc$features)
  
  pd <- cbind(pd, "n value" = colSums(!is.na(exprs(object@featureSet))), pc$samples)
  
  pData(object@featureSet) <- pd
  fData(object@featureSet) <- fd
  object@featureSet@assayData$exprs.normalized <- mat.ori
  
  object
}

