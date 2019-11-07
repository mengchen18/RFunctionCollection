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
#' @param sort could be one of c("none", "p.value", "OR") to indicate how the result should be sorted. 
#' @import fastmatch

require(fastmatch)

vectORA <- function(pathways, genelist, background, trimPathway = FALSE,
                    minOverlap = 3, minSize=5, maxSize=Inf, pathway_desc = NULL,
                    unconditional.or = TRUE, mtc.method = "fdr", 
                    sort = c("none", "p.value", "OR")[1]) {
  # check id, duplicates
  genelist <- unique(genelist)

  if (length(background) == 1 && is.integer(background)) {
    bkgn <- background 
  } else if (length(background) > 1 && is.character(background)) {
    if (!all(genelist %in% background)) {
      background <- c(background, genelist)
      message("Some IDs in genelist is not in background, background is expanded!")
    }
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
  
  pathway_annot_x <- NULL
  if (!is.null(pathway_desc))
    pathway_annot_x <- pathway_desc[names(pathways)]
  
  bdf <- vectORA.core(
    n.overlap = nol, 
    n.de = length(genelist), 
    n.gs = ngs, 
    n.bkg = bkgn, 
    unconditional.or = unconditional.or, mtc.method = mtc.method)
  rs <- cbind(
    pathway = names(pathways),
    desc = pathway_annot_x,
    bdf, 
    overlap_ids = sapply(overlap, paste0, collapse = ";")
  )
  sort <- sort[1]
  if (sort == "p.value") {
    rs <- rs[order(rs$p.value, decreasing = FALSE), ]
  } else if (sort == "OR") {
    rs <- rs[order(rs$OR, decreasing = TRUE), ]
  } else if (sort != "none") 
    warning("Unknown sort method, the results are not sorted!")
  rs
}

#' @param n.overlap the number of overlap between de and gs. The number of white balls drawn 
#'  without replacement from an urn which contains both black and white balls (compared to hyper).
#' @param n.de the number of DE gene. The number of balls drawn from the urn (compared to hyper).
#' @param n.gs the size of gene set. The number of white balls in the urn (compared to hyper).
#' @param n.bkg the background size.
#' @param unconditional.or calculate odds ratio using Maximum Likelihood Estimate (the sample odds ratio). 
#'  Note that the conditional Maximum Likelihood Estimate (MLE) is used in fisher.test. 
#' @param mtc.method multiple test correction methods, passed to p.adjust function
#' @import fastmatch
#' @examples
#' xq <- rbind(c(4, 2, 4),
#'             c(20, 40, 10),
#'             c(11, 234, 10),
#'             c(200, 1000, 100))
#' 
#' vectORA.core(xq[1, ], xq[2, ], xq[3, ], xq[4, ])
#' vectORA.core(xq[1, ], xq[2, ], xq[3, ], xq[4, ], unconditional.or = TRUE)
#' 
#' # fisher's test
#' t(apply(xq, 2, function(x1) {
#'   v <- fisher.test(rbind(c(x1[1], x1[2]-x1[1]), c(x1[3]-x1[1], x1[4] - x1[2] - x1[3] + x1[1])), alternative = "greater")
#'   c(p.value = v$p.value, v$estimate)
#' }))


vectORA.core <- function(n.overlap, n.de, n.gs, n.bkg, unconditional.or = TRUE, mtc.method = "fdr") {
  
  pval <- phyper(q = n.overlap-1, m = n.gs, n = n.bkg - n.gs, k = n.de, lower.tail = FALSE)
  
  if (unconditional.or)
    or <- (n.overlap/(n.de - n.overlap))/((n.gs-n.overlap)/(n.bkg - n.gs - n.de + n.overlap)) else {

      or <- function(n.overlap, n.gs, n.de, n.bkg) {        
        m <- n.gs
        n <- n.bkg - n.gs
        k <- n.de
        x <- n.overlap
        lo <- pmax(0L, k - n)
        hi <- pmin(k, m)
        
        supportl <- mapply(":", lo, hi, SIMPLIFY = FALSE)
        
        sapply(1:length(x), function(i) {
          support <- supportl[[i]]
          logdc <- dhyper(support, m[i], n[i], k[i], log = TRUE)
          
          dnhyper <- function(ncp) {
            d <- logdc + log(ncp) * support
            d <- exp(d - max(d))
            d/sum(d)
          }
          mnhyper <- function(ncp) {
            if (ncp == 0) 
              return(lo[i])
            if (ncp == Inf) 
              return(hi[i])
            sum(support * dnhyper(ncp))
          }
          mle <- function(x) {
            if (x == lo[i]) 
              return(0)
            if (x == hi[i]) 
              return(Inf)
            
            mu <- mnhyper(1)
            if (mu > x) 
              uniroot(function(t) mnhyper(t) - x, c(0, 1))$root
            else if (mu < x) 
              1/uniroot(function(t) mnhyper(1/t) - x, c(.Machine$double.eps, 1))$root
            else 1
          }
          mle(x[i])
        })  
      }
      or <- or(n.overlap=n.overlap, n.gs=n.gs, n.de=n.de, n.bkg=n.bkg)
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





