multi.t.test2 <- function(x, pheno, compare = NULL, log10 = FALSE, median.center = TRUE, fillNA = TRUE, ...) {
  
  x.raw <- x
  if (is.null(rownames(x)))
    rownames(x) <- 1:nrow(x)
  
  halfValue <- function(x) x/2  
  if (log10) {
    x <- apply(x, 2, log10)
    halfValue <- function(x) x - log10(2)
  }
  
  if (median.center)
    x <- sweep(x, 2, matrixStats::colMedians(x, na.rm = TRUE), "-") + median(x, na.rm = TRUE)
  
  if (fillNA) {
    x <- apply(x, 1, function(xx) {
      x3 <- xx
      x3[is.na(x3)] <- halfValue(min(x3, na.rm = TRUE))
      x3
    })
    x <- t(x)
  }
  
  if (log10)
    x.ori <- 10^x else
      x.ori <- x
  
  if (is.null(compare)) {
    return( list(ttest = NULL, mat = x, mat.rawscale = x.ori) )
  } else {
    tl <- lapply(unique(compare[, 1]), function(x) {
      x <- compare[x == compare[, 1], -1, drop = FALSE]
      unique(unlist(split(x, row(x))))
    })
    names(tl) <- unique(compare[, 1])
  }
  
  df <- data.frame(metabolite = rownames(x), stringsAsFactors = FALSE)
  for ( i in names(tl) ) {
    for ( j in tl[[i]] ) {
      m <- x[, pheno[[i]] == j]
      rv <- rowSums(!is.na(x.raw[, pheno[[i]] == j]))
      rm <- rowMeans(m, na.rm = TRUE)
      rq <- rank(rm, na.last = TRUE)/sum(!is.na(rm))
      rq[is.na(rm)] <- NA
      df[[paste("Stats|Mean", j, sep = "|")]] <- rm
      df[[paste("Stats|N value", j, sep = "|")]] <- rv
      df[[paste("Stats|Quantile", j, sep = "|")]] <- rq
    }
  }  
  for ( i in 1:nrow(compare) ) {
    v <- compare[i, ]
    i1 <- pheno[[v[1]]] == v[2]
    i2 <- pheno[[v[1]]] == v[3]
    
    tv <- apply(x, 1, function(xx) {
      t <- try(t.test(xx[i1], xx[i2], ...), silent = TRUE)
      if (class(t) != "htest")
        return(c(pvalue = NA, mean.diff = NA))
      if (length(t$estimate) == 1)
        md <- t$estimate[[1]] else
          md <- t$estimate[1] - t$estimate[2]
      c(pvalue = t$p.value, mean.diff = md)
    })
    
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "pvalue", sep = "|")]] <- tv[1, ]
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "log.pvalue", sep = "|")]] <- -log10(tv[1, ])
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "fdr", sep = "|")]] <- p.adjust(tv[1, ], method = "fdr")
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "log.fdr", sep = "|")]] <- -log10(df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "fdr", sep = "|")]])
    df[[paste("ttest", paste(v[2], v[3], sep = "_vs_"), "mean.diff", sep = "|")]] <- tv[2, ]
  }  
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
      v <- pca(x, n = n, prefix = "PCA|All")
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

prepViewerData <- function(
  object, pheno = NULL, log10 = TRUE, median.center = FALSE, fillNA = TRUE,
  compare.t.test = NULL, compare.pca = list("_all_" = TRUE), nf = 6) {
  
  # phenotype data
  if (is.null(pheno))
    pd <- Biobase::pData(object@featureSet) else
      pd <- pheno
    if (!is.null(pd$label))
      pd$file <- pd$label
    colnames(pd) <- paste("General|All", colnames(pd), sep = "|")
    
    # feature data
    fd <- fData(object@featureSet)
    gac <- c("ID","rtmed","mzmed","annot_ms1", "annot_ms2")
    ga <- fd[, gac]
    ea <- fd[, setdiff(colnames(fd), gac)]
    colnames(ga) <- paste("General|All", colnames(ga), sep = "|")
    colnames(ea) <- paste("General|Extended", colnames(ea), sep = "|")
    fdx <- cbind(ga, ea)
    rownames(fdx) <- rownames(fd)
    
    ######## t-test #######
    ts <- multi.t.test2(
      x = Biobase::exprs(object@featureSet), 
      compare = compare.t.test, 
      pheno = pd,
      log10 = log10, 
      median.center = median.center, 
      fillNA = fillNA)
    
    mat <- ts$mat
    mat.ori <- ts$mat.rawscale
    ts <- ts$ttest
    
    ######## PCA #######
    pc <- multi.pca(mat, pheno = pd, compare = compare.pca, n = nf)
    fd <- cbind(fdx, pc$features)
    if (!is.null(ts))
      fd <- cbind(fd, ts[-1])
    
    pd <- cbind(pd, "Stats|All|n value" = colSums(!is.na(exprs(object@featureSet))), pc$samples)
    
    pData(object@featureSet) <- pd
    fData(object@featureSet) <- fd
    exprs(object@featureSet) <- mat
    
    attr(d, "fx") <- "General|All|rtmed"
    attr(d, "fy") <- "General|All|mzmed"
    attr(d, "sx") <- "PCA|AllSample|PC1("
    attr(d, "sy") <- "PCA|AllSample|PC2("
    
    object
}
