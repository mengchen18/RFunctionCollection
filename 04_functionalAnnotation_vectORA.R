#' three columns in 
#' @param gene.id the gene id, should be the same length as p.value. It should be unique and not delimited for multiple genes.
#' @param p.val the p.val from a test, usually a test for differential expression
#' @param geneset a character vector of indicates which geneset are associated with a specific gene.id, should be delimited by semicolon (;).
#'   It should have the same length as p.value
#' @param p.thresh the threshold of p value, lower than this value would be considered as significant
#' @param min.size the minimum size of geneset that should be considered
#' 
#' @import fastmatch

vectORA <- function(p.val, gene.id, geneset, p.thresh = 0.05, min.size = 3) {
  
  library(fastmatch)
  if (any(duplicated(gene.id)))
    stop("duplicated gene IDs")
  
  go <- geneset
  names(go) <- gene.id
  go <- strsplit(go, split = " ")
  go <- lapply(names(go), function(x) cbind(rep(x, length(go[[x]])), go[[x]]))
  go <- do.call(rbind, go)
  go <- go[!go[, 2] %in% c("available", "GO:not"), ]
  go <- split(go[, 1], go[, 2])
  go <- go[sapply(go, length) >= min.size]
  
  i <- p.val < p.thresh
  
  sigID <- gene.id[i]
  nsigID <- length(sigID)
  nid <- length(gene.id)
  
  rmat <- sapply(names(go), function(x) {
    gs1 <- go[[x]]
    ii <- gs1 %in% sigID
    inin <- sum(ii)
    ctab <- rbind(c(inin, nsigID - inin), 
                  c(length(gs1) - inin, nid - length(gs1) - nsigID + inin))
    t <- fisher.test(ctab, alternative = "greater")
    
    c(geneset = x,
      p.value = t$p.value,
      OR = t$estimate[[1]],
      overlap_size = inin,
      overlap_id = paste(gs1[ii], collapse = ";")
      )
  })
  df <- data.frame(
    geneset = rmat["geneset" , ],
    p.value = as.numeric(rmat["p.value", ]),
    fdr = NA,
    OR = as.numeric(rmat["OR", ]),
    overlap_size = as.numeric(rmat["overlap_size", ]),
    overlap_id = rmat["overlap_id", ],
    row.names = NULL,
    stringsAsFactors = FALSE)
  df$fdr <- p.adjust(df$p.value, method = "fdr")
  df <- df[order(df$OR, decreasing = TRUE), ]
  df
}
