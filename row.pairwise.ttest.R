#' row-wise pair wise t test
#' @param x expression matrix where rows are features and columns are samples
#' @param g group vector, has the same length as column number of x
#' @param cut if all value in a comparison lower than cut, the test won't be performed
#' @param p.adjust p.adjust method passed to p.adjust function
#' @param ... passed to t.test

row.pairwise.ttest <- function(x, g, cut=-Inf, p.adjust = p.adjust.methods, comb.adjust = FALSE, ...) {
  ord <- order(g)
  g <- g[ord]
  x <- x[, ord]
  
  lv <- unique(g)
  
  if (length(lv) > 5)
    warning("More than 5 groups, there are 10 pairwise comparisoins. Please consider use other method.")
    
  cbn <- combn(lv, 2)
  
  l1 <- lapply(1:ncol(cbn), function(i) {
    ii <- g %in% cbn[, i]
    cat(paste("Testing", cbn[1, i], "vs", cbn[2, i], "...\n"))
    gg <- g[ii]
    tres <- apply(x[, ii], 1, function(xx) {
      if (max(xx, na.rm = TRUE) < cut)
        r <- c(NA, NA) else {
          tt <- try(t.test(xx ~ gg, ...), silent = TRUE)
          if (class(tt) == "try-error")
            r <- c(NA, NA) else
              r <- c(p.val = tt$p.value, fc = tt$estimate[1] - tt$estimate[2])
        }
      r
    })
    data.frame(pval = tres[1, ], 
               fc = tres[2, ], 
               padj = p.adjust(tres[1, ], method = p.adjust))
  })
  names(l1) <- paste(cbn[1, ], cbn[2, ], sep = "_")
  
  ntests <- sum(sapply(l1, function(x) sum(!is.na(x$pval))))
  l1 <- lapply(l1, function(x) {
    x$padj.comb <- p.adjust(x$pval, method = p.adjust, n = ntests)
    x
  })
  la <- do.call(cbind, l1)
  rownames(la) <- rownames(x)
  la
}
