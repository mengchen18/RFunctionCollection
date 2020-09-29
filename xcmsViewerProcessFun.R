

multi.t.test2 <- function(x, pheno, compare = NULL, log10 = FALSE, median.center = TRUE, fillNA = TRUE, ...) {
  
  if (is.null(rownames(x)))
    rownames(x) <- 1:nrow(x)
  
  halfValue <- function(x) x/2
  if (log10) {
    x <- apply(x, 2, log10)
    halfValue <- function(x) x - log10(2)
  }
  
  if (is.null(compare)) {
    tl <- NULL
  } else {
    tl <- lapply(unique(compare[, 1]), function(x) {
      x <- compare[x == compare[, 1], -1]
      unique(unlist(split(x, row(x))))
    })
    names(tl) <- unique(compare[, 1])
  }
  
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
      df[[paste("mean", j, sep = "|")]] <- rm
      df[[paste("n value", j, sep = "|")]] <- rv
      df[[paste("quantile", j, sep = "|")]] <- rq
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
pca <- function(x, n = 6, prefix = "") {
  
  writePC <- function(x, n) {
    var <- round(x$sdev[1:n]^2/(sum(x$sdev^2))*100, digits = 1)
    xx <- x$x[, 1:min(n, ncol(x$x))]
    colnames(xx) <- paste0(prefix, "|", colnames(xx), " (", var, "%", ")")
    pp <- x$rotation[, 1:min(n, ncol(x$x))]
    colnames(pp) <- paste0(prefix, "|", colnames(pp), " (", var, "%", ")")
    list(samples = xx, features = pp)
  }
  
  pc <- prcomp(t(x))
  writePC(pc, n = n)
}


multi.pca <- function(x, pheno, compare, n = 6) {
  
  cn <- setdiff(names(compare), c("_all_", colnames(pheno)))
  if (length(cn) > 0)
    stop(sprintf("These columns do not exist in pheno: %s", paste(cn, collapse = ", ")))
  
  if (!inherits(x, "matrix"))
    x <- apply(x, 2, as.numeric)
  
  temp_sample <- matrix(NA, nrow = ncol(x), ncol = n)
  rownames(temp_sample) <- rownames(pheno)
  temp_feature <- matrix(NA, nrow = nrow(x), ncol = n)
  rownames(temp_feature) <- rownames(x)
  
  t <- lapply(names(compare), function(i, temp_sample, temp_feature) {
    if (i == "_all_") {
      v <- pca(x, n = n, prefix = "PCA|AllSamples")
    } else {
      ph <- pheno[[i]]
      
      sd <- setdiff(compare[[i]], ph)
      if (length(sd) > 0)
        stop(sprintf("Column '%s' in 'pheno' do not contain these vaues: %s", i, paste(sd, collapse = ", ")))
      
      ii <- which(ph %in% compare[[i]])
      x0 <- x[, ii]
      ir <- which(matrixStats::rowVars(x0, na.rm = TRUE) > 0)
      x0 <- x0[ir, ]
      
      pc <- pca(x0, n = n, prefix = paste(c("PCA", compare[[i]]), collapse = "|"))
      temp_sample[ii, ] <- pc$samples
      colnames(temp_sample) <- colnames(pc$samples)
      temp_feature[ir, ] <- pc$features
      colnames(temp_feature) <- colnames(pc$features)
      v <- list(samples = temp_sample, features = temp_feature)  
    }
    v
  },  temp_sample = temp_sample, temp_feature = temp_feature)
  
  list(
    samples = data.frame(do.call(cbind, lapply(t, "[[", "samples")), stringsAsFactors = FALSE, check.names = FALSE),
    features = data.frame(do.call(cbind, lapply(t, "[[", "features")), stringsAsFactors = FALSE, check.names = FALSE)
  )
}


  

###
phenoFeatureData <- function(
  object, compare.t.test = NA, compare.pca = list("_all_" = TRUE),
  pheno=NULL, log10 = TRUE, median.center = TRUE, fillNA = TRUE, 
  nf = 6, ...
) {
  
  if (is.null(pheno))
    pheno <- Biobase::pData(object@featureSet)
  
  ts <- multi.t.test2(
    x = Biobase::exprs(object@featureSet), 
    compare = compare.t.test, 
    pheno = pheno,
    log10 = log10, 
    median.center = median.center, 
    fillNA = fillNA, ...)
  
  mat <- ts$mat
  mat.ori <- ts$mat.rawscale
  ts <- ts$ttest
  
  pc <- multi.pca(mat, pheno = pheno, compare = compare.pca, n = nf)
  
  fd <- cbind(fData(object@featureSet), ts[-1], pc$features)
  pd <- cbind(pd, "n value" = colSums(!is.na(exprs(object@featureSet))), pc$samples)
  
  updateDF <- function(x, d) {
    ii <- intersect(colnames(x), colnames(d))
    ia <- setdiff(colnames(d), colnames(x))
    x[, ii] <- d[, ii]
    x <- cbind(x, d[, ia])
    x
  }
  
  pData(object@featureSet) <- updateDF(pData(object@featureSet), pd)
  fData(object@featureSet) <- updateDF(fData(object@featureSet), fd)
  object@featureSet@assayData$exprs.normalized <- mat.ori
  
  object
}

