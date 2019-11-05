#' three columns in 
#' @param gene.id the gene id, should be the same length as p.value. It should be unique and not delimited for multiple genes.
#' @param p.val the p.val from a test, usually a test for differential expression
#' @param geneset a character vector of indicates which geneset are associated with a specific gene.id, should be delimited by semicolon (;).
#'   It should have the same length as p.value
#' @param p.thresh the threshold of p value, lower than this value would be considered as significant
#' @param min.size the minimum size of geneset that should be considered
#' 
#' @import fastmatch


library(fastmatch)

vectORA <- function(pathways, pvec, pvec.cut, pvecAnnot=NULL, 
                    minOverlap = 3, minSize=5, maxSize=Inf, size_background = NULL,
                    unconditional.or = TRUE, mtc.method = "fdr") {
  
  if (is.null(pvecAnnot))
    pvecAnnot <- names(pvec)
  
  allid <- unique(unlist(pvecAnnot))
  allid.pathways <- unique(unlist(pathways))
  bkgn <- sum(allid %fin% allid.pathways)
  if (!is.null(size_background)) {
    bkgn <- max(size_background, bkgn)
    if (size_background < bkgn)
      message(paste0("size_background is set to ", bkgn, "!" ))
  }
  
  deid <- unique(unlist(pvecAnnot[pvec < pvec.cut]))
  overlap <- lapply(pathways, function(x) x[x %fin% deid])
  
  nol <- sapply(overlap, length)
  ngs <- sapply(pathways, length)
  i <- nol >=  minOverlap & ngs >= minSize & ngs <= maxSize
  
  pathways <- pathways[i]
  overlap <- overlap[i]
  nol <- nol[i]
  ngs <- ngs[i]
  
  
  bdf <- vectORA.core(
    n.overlap = nol, 
    n.de = length(deid), 
    n.gs = ngs, 
    n.bkg = bkgn, 
    unconditional.or = unconditional.or, mtc.method = mtc.method)
  cbind(
    pathway = names(pathways),
    bdf, 
    overlap_ids = sapply(overlap, paste0, collapse = ";")
  )
}


vectORA.core <- function(n.overlap, n.de, n.gs, n.bkg, unconditional.or = TRUE, mtc.method = "fdr") {
  
  pval <- phyper(q = n.overlap-1, m = n.gs, n = n.bkg - n.gs, k = n.de, lower.tail = FALSE)
  
  if (unconditional.or)
    or <- (n.de/(n.de - n.overlap))/((n.gs-n.overlap)/(n.bkg - n.gs - n.de + n.overlap)) else {
      stop("not implemneted yet!")
    }
  
  data.frame(
    p.value = pval,
    p.adjusted = p.adjust(pval, method = mtc.method),
    OR = or,
    size_overlap = n.overlap,
    size_geneset = n.gs,
    size_input = n.de,
    size_backgroung = n.bkg,
    stringsAsFactors = FALSE
  )
}

# # example
# xq <- rbind(c(4, 2, 4), 
#             c(20, 40, 10),
#             c(11, 234, 10),
#             c(200, 1000, 100))
# 
# vectORA(xq[1, ], xq[2, ], xq[3, ], xq[4, ])
# 
# apply(xq, 2, function(x1) {
#   fisher.test(rbind(c(x1[1], x1[2]-x1[1]), c(x1[3]-x1[1], x1[4] - x1[2] - x1[3] + x1[1])), alternative = "greater")$p.value
# })



