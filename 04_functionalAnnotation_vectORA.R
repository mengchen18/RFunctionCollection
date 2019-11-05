#' three columns in 
#' @param pathways List of gene sets to check.
#' @param genelist A list of gene/proteins you want to annotate, often from differential expression analysis.
#'  The ID needs to be the same as in pathways
#' @param background An integer to indicate the size of background; or a character vector 
#'  stores all background genes, the ID needs to be the same as in pathways
#' @param trimPathway remove ID from pathways that are not present in the background, 
#'  only used if background is a gene/protein name list.
#' @param minOverlap the minimum required overlap between pathway and gene list, if the overlap is lower
#'  than this value, no test would be done on this pathway
#' @param minSize the minimum size of pathway should be tested
#' @param maxSize the maximum size of pathways should be tested
#' @param pathway_desc description of pathways, a name character vector storing the description of
#'  of the pathway, the names should be the same as names in "pathways" argument.
#' @param unconditional.or calculate odds ratio using Maximum Likelihood Estimate (the sample odds ratio). 
#'  Note that the conditional Maximum Likelihood Estimate (MLE) is used in fisher.test. 
#' @param mtc.method multiple test correction methods, passed to p.adjust function
#' @import fastmatch
require(fastmatch)

vectORA <- function(pathways, genelist, background, trimPathway = FALSE,
                    minOverlap = 3, minSize=5, maxSize=Inf, pathway_desc = NULL,
                    unconditional.or = TRUE, mtc.method = "fdr") {
  # check id, duplicates
  genelist <- unique(genelist)
  if (length(background) == 1 && is.integer(background)) {
  	bkgn <- background 
  } else if (length(background) > 1 && is.character(background)) {
  	background <- unique(background)
  	bkgn <- length(background)
  	if (trimPathway)
  		pathways <- lapply(pathways, function(x) x[x %fin% background]) 
  } else 
  	stop("Unknown background type!")
  
  overlap <- lapply(pathways, function(x) x[x %fin% genelist])
  
  nol <- sapply(overlap, length)
  ngs <- sapply(pathways, length)
  i <- nol >=  minOverlap & ngs >= minSize & ngs <= maxSize
  pathways <- pathways[i]
  overlap <- overlap[i]
  nol <- nol[i]
  ngs <- ngs[i]
  
  if (!is.null(pathway_desc))
  	pathway_annot_x <- pathway_desc[names(pathways)]

  bdf <- vectORA.core(
    n.overlap = nol, 
    n.de = length(genelist), 
    n.gs = ngs, 
    n.bkg = bkgn, 
    unconditional.or = unconditional.or, mtc.method = mtc.method)
  cbind(
    pathway = names(pathways),
    desc = pathway_annot_x,
    bdf, 
    overlap_ids = sapply(overlap, paste0, collapse = ";")
  )
}


vectORA.core <- function(n.overlap, n.de, n.gs, n.bkg, unconditional.or = TRUE, mtc.method = "fdr") {
  
  pval <- phyper(q = n.overlap-1, m = n.gs, n = n.bkg - n.gs, k = n.de, lower.tail = FALSE)
  
  if (unconditional.or)
    or <- (n.overlap/(n.de - n.overlap))/((n.gs-n.overlap)/(n.bkg - n.gs - n.de + n.overlap)) else {
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



