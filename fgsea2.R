library(fgsea)
library(BiocParallel)
library(fastmatch)
library(stats)
library(ggplot2)
library(gridExtra)
library(grid)
library(data.table)

#' list matching
#'
#' returns a vector of the positions of all matches of a list
#' it is used by fgsea2
#' @param x the values to be matched
#' @param list a list or vector to be matched
#' @import fastmatch
#' @examples
#'  listmatch(c("x", "z"), letters)
#'  listmatch(c("x", "z"), list(letters[1:10], letters))
listmatch <- function(x, list) {
    which(sapply(list, function(xx) any(xx %fin% x)))
}

#' Runs preranked gene set enrichment analysis on variables linked to multiple IDs.
#'
#' The function is modified from function \code{fgsea} in package "fgsea", 
#' but a new argument "statsAnnot" is
#' added to support running GSEA on variables linked to multiple IDs, such as proteomics or 
#' PTM data.
#' @param pathways List of gene sets to check.
#' @param stats Named vector of gene-level stats. Names should be the same as in 'pathways'
#' @param nperm Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
#' @param statsAnnot a list has the same length as stats. If this argument is given the IDs in this list
#' 	will be mapped to pathways, so one gene name in pathways could be mapped to multiple 
#'  variables in stats. 
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param nproc If not equal to zero sets BPPARAM to use nproc workers (default = 0).
#' @param gseaParam GSEA parameter value, all gene-level statis are raised to the power of `gseaParam`
#'                  before calculation of GSEA enrichment scores.
#' @param BPPARAM Parallelization parameter used in bplapply.
#'  Can be used to specify cluster to run. If not initialized explicitly or
#'  by setting `nproc` default value `bpparam()` is used.
#' @return A table with GSEA results. Each row corresponds to a tested pathway.
#' The columns are the following:
#' \itemize{
#'  \item pathway -- name of the pathway as in `names(pathway)`;
#'  \item pval -- an enrichment p-value;
#'  \item padj -- a BH-adjusted p-value;
#'  \item ES -- enrichment score, same as in Broad GSEA implementation;
#'  \item NES -- enrichment score normalized to mean enrichment of random samples of the same size;
#'  \item nMoreExtreme` -- a number of times a random gene set had a more
#'      extreme enrichment score value;
#'  \item size -- size of the pathway after removing genes not present in `names(stats)`.
#'  \item leadingEdge -- vector with indexes of leading edge genes that drive the enrichment, see \url{http://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_a_Leading}.
#' }
#'
#' @export
#' @import data.table
#' @import BiocParallel
#' @import fastmatch
#' @import stats
#' @import fgsea
#' @examples
#' ############ example in fgsea package ###############
#' data(examplePathways)
#' data(exampleRanks)
#' 
#' fgseaRes2 <- fgsea2(examplePathways[1:10], 
#'                     exampleRanks, 
#'                     statsAnnot = names(exampleRanks), 
#'                     nperm=1000, maxSize=500)
#' 
#' plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]], exampleRanks) + 
#'   labs(title="Programmed Cell Death")
#' 
#' plotEnrichment2(examplePathways[["5991130_Programmed_Cell_Death"]], exampleRanks, 
#'                 statsAnnot = names(exampleRanks)) + 
#'   labs(title="Programmed Cell Death")
#' 
#' plotGseaTable2(pathways = examplePathways[1:10], stats = exampleRanks, 
#'                statsAnnot=names(exampleRanks), fgseaRes = fgseaRes2)
#' 
#' 
#' 
#' ############ example simulated data ###############
#' stats <- 1:20 + 0.5
#' names(stats) <- letters[1:20]
#' stannot <- strsplit(names(stats), "")
#' stannot[1:5] <- lapply(stannot[1:5], c, "a")
#' stannot[16:20] <- lapply(stannot[16:20], c, "t")
#' gs <- list(top = c("a"),   # positive control - ES > 0
#'            bottom = c("t"), # positive control - ES < 0
#'            r1 = "q", # negative control 
#'            r2 = "h" # negative control
#' )
#' res <- fgsea2(gs, stats, statsAnnot = stannot, nperm=10000, minSize = 1)
#' 
#' plotEnrichment2(pathway = gs[["top"]], stats = stats, statsAnnot = stannot)
#' plotEnrichment2(pathway = gs[["bottom"]], stats = stats, statsAnnot = stannot)
#' plotEnrichment2(pathway = gs[["r1"]], stats = stats, statsAnnot = stannot)
#' plotEnrichment2(pathway = gs[["r2"]], stats = stats, statsAnnot = stannot)
#' 
#' plotGseaTable2(pathways = gs, stats = stats, statsAnnot=stannot, fgseaRes = res)
fgsea2 <- function(pathways, stats, nperm,
				   statsAnnot=NULL, 
                   minSize=1, maxSize=Inf,
                   nproc=0,
                   gseaParam=1,
                   BPPARAM=NULL) {

    # Warning message for ties in stats
    ties <- sum(duplicated(stats[stats != 0]))
    if (ties != 0) {
        warning("There are ties in the preranked stats (",
                paste(round(ties * 100 / length(stats), digits = 2)),
                "% of the list).\n",
                "The order of those tied genes will be arbitrary, which may produce unexpected results.")
    }

    # Warning message for duplicate gene names
    if (any(duplicated(names(stats)))) {
        warning("There are duplicate gene names, fgsea may produce unexpected results")
    }

    granularity <- 1000
    permPerProc <- rep(granularity, floor(nperm / granularity))
    if (nperm - sum(permPerProc) > 0) {
        permPerProc <- c(permPerProc, nperm - sum(permPerProc))
    }
    seeds <- sample.int(10^9, length(permPerProc))

    if (is.null(BPPARAM)) {
        if (nproc != 0) {
            if (.Platform$OS.type == "windows") {
                # windows doesn't support multicore, using snow instead
                BPPARAM <- SnowParam(workers = nproc)
            } else {
                BPPARAM <- MulticoreParam(workers = nproc)
            }
        } else {
            BPPARAM <- bpparam()
        }
    }

    minSize <- max(minSize, 1)
    ### ===============================================
    if (is.null(statsAnnot))
        statsAnnot <- names(stats)
    # stats <- sort(stats, decreasing = TRUE)
    # stats <- abs(stats)^gseaParam
    sord <- order(stats, decreasing = TRUE)
    stats <- stats[sord]
    stats <- abs(stats)^gseaParam 
    statsAnnot <- statsAnnot[sord]   
    # pathwaysFiltered <- lapply(pathways, function(p) {
    #     as.vector(na.omit(fmatch(p, names(stats))))
    # })
    pathwaysFiltered <- lapply(pathways, function(p) {
        as.vector(listmatch(p, statsAnnot))
    })
    ### ===============================================
    pathwaysSizes <- sapply(pathwaysFiltered, length)

    toKeep <- which(minSize <= pathwaysSizes & pathwaysSizes <= maxSize)
    m <- length(toKeep)

    if (m == 0) {
        return(data.table(pathway=character(),
                          pval=numeric(),
                          padj=numeric(),
                          ES=numeric(),
                          NES=numeric(),
                          nMoreExtreme=numeric(),
                          size=integer(),
                          leadingEdge=list()))
    }

    pathwaysFiltered <- pathwaysFiltered[toKeep]
    pathwaysSizes <- pathwaysSizes[toKeep]

    K <- max(pathwaysSizes)

    gseaStatRes <- do.call(rbind,
                lapply(pathwaysFiltered, calcGseaStat,
                       stats=stats,
                       returnLeadingEdge=TRUE))


    leadingEdges <- mapply("[", list(names(stats)), gseaStatRes[, "leadingEdge"], SIMPLIFY = FALSE)
    pathwayScores <- unlist(gseaStatRes[, "res"])



    universe <- seq_along(stats)

    counts <- bplapply(seq_along(permPerProc), function(i) {
        nperm1 <- permPerProc[i]
        leEs <- rep(0, m)
        geEs <- rep(0, m)
        leZero <- rep(0, m)
        geZero <- rep(0, m)
        leZeroSum <- rep(0, m)
        geZeroSum <- rep(0, m)
        if (m == 1) {
            for (i in seq_len(nperm1)) {
                randSample <- sample.int(length(universe), K)
                randEsP <- calcGseaStat(
                    stats = stats,
                    selectedStats = randSample,
                    gseaParam = 1)
                leEs <- leEs + (randEsP <= pathwayScores)
                geEs <- geEs + (randEsP >= pathwayScores)
                leZero <- leZero + (randEsP <= 0)
                geZero <- geZero + (randEsP >= 0)
                leZeroSum <- leZeroSum + pmin(randEsP, 0)
                geZeroSum <- geZeroSum + pmax(randEsP, 0)
            }
        } else {
            aux <- fgsea:::calcGseaStatCumulativeBatch(
                stats = stats,
                gseaParam = 1,
                pathwayScores = pathwayScores,
                pathwaysSizes = pathwaysSizes,
                iterations = nperm1,
                seed = seeds[i])
            leEs = get("leEs", aux)
            geEs = get("geEs", aux)
            leZero = get("leZero", aux)
            geZero = get("geZero", aux)
            leZeroSum = get("leZeroSum", aux)
            geZeroSum = get("geZeroSum", aux)
        }
        data.table(pathway=seq_len(m),
                   leEs=leEs, geEs=geEs,
                   leZero=leZero, geZero=geZero,
                   leZeroSum=leZeroSum, geZeroSum=geZeroSum
                   )
    }, BPPARAM=BPPARAM)

    counts <- rbindlist(counts)

    # Getting rid of check NOTEs
    leEs=leZero=geEs=geZero=leZeroSum=geZeroSum=NULL
    pathway=padj=pval=ES=NES=geZeroMean=leZeroMean=NULL
    nMoreExtreme=nGeEs=nLeEs=size=NULL
    leadingEdge=NULL
    .="damn notes"


    pvals <- counts[,
        list(pval=min((1+sum(leEs)) / (1 + sum(leZero)),
                 (1+sum(geEs)) / (1 + sum(geZero))),
             leZeroMean = sum(leZeroSum) / sum(leZero),
             geZeroMean = sum(geZeroSum) / sum(geZero),
             nLeEs=sum(leEs),
             nGeEs=sum(geEs)
             )
        ,
        by=.(pathway)]
    pvals[, padj := p.adjust(pval, method="BH")]
    pvals[, ES := pathwayScores[pathway]]
    pvals[, NES := ES / ifelse(ES > 0, geZeroMean, abs(leZeroMean))]
    pvals[, leZeroMean := NULL]
    pvals[, geZeroMean := NULL]

    pvals[, nMoreExtreme :=  ifelse(ES > 0, nGeEs, nLeEs)]
    pvals[, nLeEs := NULL]
    pvals[, nGeEs := NULL]

    pvals[, size := pathwaysSizes[pathway]]
    pvals[, pathway := names(pathwaysFiltered)[pathway]]

    pvals[, leadingEdge := .(leadingEdges)]


    # Makes pvals object printable immediatly
    pvals <- pvals[]
    pvals$leadingEdge <- sapply(pvals$leadingEdge, paste, collapse = ";")
    pvals
}


#' Plots GSEA enrichment plot for function fgsea2.
#' @param pathway Gene set to plot.
#' @param stats Gene-level statistics.
#' @param statsAnnot a list has the same length as stats. If this argument is given the IDs in this list
#' 	will be mapped to pathways, so one gene name in pathways could be mapped to multiple 
#'  variables in stats. 
#' @param gseaParam GSEA parameter.
#' @param ticksSize width of vertical line corresponding to a gene (default: 0.2)
#' @return ggplot object with the enrichment plot.
#' @export
#' @examples
#' ############ example in fgsea package ###############
#' data(examplePathways)
#' data(exampleRanks)
#' 
#' fgseaRes2 <- fgsea2(examplePathways[1:10], 
#'                     exampleRanks, 
#'                     statsAnnot = names(exampleRanks), 
#'                     nperm=1000, maxSize=500)
#' 
#' plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]], exampleRanks) + 
#'   labs(title="Programmed Cell Death")
#' 
#' plotEnrichment2(examplePathways[["5991130_Programmed_Cell_Death"]], exampleRanks, 
#'                 statsAnnot = names(exampleRanks)) + 
#'   labs(title="Programmed Cell Death")
#' 
#' plotGseaTable2(pathways = examplePathways[1:10], stats = exampleRanks, 
#'                statsAnnot=names(exampleRanks), fgseaRes = fgseaRes2)
#' 
#' 
#' 
#' ############ example simulated data ###############
#' stats <- 1:20 + 0.5
#' names(stats) <- letters[1:20]
#' stannot <- strsplit(names(stats), "")
#' stannot[1:5] <- lapply(stannot[1:5], c, "a")
#' stannot[16:20] <- lapply(stannot[16:20], c, "t")
#' gs <- list(top = c("a"),   # positive control - ES > 0
#'            bottom = c("t"), # positive control - ES < 0
#'            r1 = "q", # negative control 
#'            r2 = "h" # negative control
#' )
#' res <- fgsea2(gs, stats, statsAnnot = stannot, nperm=10000, minSize = 1)
#' 
#' plotEnrichment2(pathway = gs[["top"]], stats = stats, statsAnnot = stannot)
#' plotEnrichment2(pathway = gs[["bottom"]], stats = stats, statsAnnot = stannot)
#' plotEnrichment2(pathway = gs[["r1"]], stats = stats, statsAnnot = stannot)
#' plotEnrichment2(pathway = gs[["r2"]], stats = stats, statsAnnot = stannot)
#' 
#' plotGseaTable2(pathways = gs, stats = stats, statsAnnot=stannot, fgseaRes = res)
plotEnrichment2 <- function(pathway, stats, statsAnnot = NULL,
                          gseaParam=1,
                          ticksSize=0.2) {

    rnk <- rank(-stats)
    ord <- order(rnk)

    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
    statsAdj <- statsAdj / max(abs(statsAdj))

    ### ===============================================
    # pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    if (is.null(statsAnnot))
        statsAnnot <- names(stats)
    statsAnnot <- statsAnnot[ord]    
    pathway <- listmatch(pathway, statsAnnot)
    ### ===============================================
    pathway <- sort(pathway)

    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                            returnAllExtremes = TRUE)

    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops

    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))

    diff <- (max(tops) - min(bottoms)) / 8

    # Getting rid of NOTEs
    x=y=NULL
    g <- ggplot(toPlot, aes(x=x, y=y)) +
        geom_point(color="green", size=0.1) +
        geom_hline(yintercept=max(tops), colour="red", linetype="dashed") +
        geom_hline(yintercept=min(bottoms), colour="red", linetype="dashed") +
        geom_hline(yintercept=0, colour="black") +
        geom_line(color="green") + theme_bw() +
        geom_segment(data=data.frame(x=pathway),
                     mapping=aes(x=x, y=-diff/2,
                                 xend=x, yend=diff/2),
                     size=ticksSize) +

        theme(panel.border=element_blank(),
              panel.grid.minor=element_blank()) +

        labs(x="rank", y="enrichment score")
    g
}



#' Plots table of enrichment graphs using ggplot and gridExtra for fgsea2.
#' @param pathways Pathways to plot table, as in `fgsea` function.
#' @param stats Gene-level stats, as in `fgsea` function.
#' @param fgseaRes Table with fgsea results.
#' @param statsAnnot a list has the same length as stats. If this argument is given the IDs in this list
#' 	will be mapped to pathways, so one gene name in pathways could be mapped to multiple 
#'  variables in stats. 
#' @param gseaParam GSEA-like parameter. Adjusts displayed statistic values,
#'      values closer to 0 flatten plots. Default = 1, value of 0.5 is a good
#'      choice too.
#' @param colwidths Vector of five elements corresponding to column width for
#'      grid.arrange. If column width is set to zero, the column is not drawn.
#' @return TableGrob object returned by grid.arrange.
#' @import ggplot2
#' @import gridExtra
#' @import grid
#' @export
#' @examples
#' @examples
#' ############ example in fgsea package ###############
#' data(examplePathways)
#' data(exampleRanks)
#' 
#' fgseaRes2 <- fgsea2(examplePathways[1:10], 
#'                     exampleRanks, 
#'                     statsAnnot = names(exampleRanks), 
#'                     nperm=1000, maxSize=500)
#' 
#' plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]], exampleRanks) + 
#'   labs(title="Programmed Cell Death")
#' 
#' plotEnrichment2(examplePathways[["5991130_Programmed_Cell_Death"]], exampleRanks, 
#'                 statsAnnot = names(exampleRanks)) + 
#'   labs(title="Programmed Cell Death")
#' 
#' plotGseaTable2(pathways = examplePathways[1:10], stats = exampleRanks, 
#'                statsAnnot=names(exampleRanks), fgseaRes = fgseaRes2)
#' 
#' 
#' 
#' ############ example simulated data ###############
#' stats <- 1:20 + 0.5
#' names(stats) <- letters[1:20]
#' stannot <- strsplit(names(stats), "")
#' stannot[1:5] <- lapply(stannot[1:5], c, "a")
#' stannot[16:20] <- lapply(stannot[16:20], c, "t")
#' gs <- list(top = c("a"),   # positive control - ES > 0
#'            bottom = c("t"), # positive control - ES < 0
#'            r1 = "q", # negative control 
#'            r2 = "h" # negative control
#' )
#' res <- fgsea2(gs, stats, statsAnnot = stannot, nperm=10000, minSize = 1)
#' 
#' plotEnrichment2(pathway = gs[["top"]], stats = stats, statsAnnot = stannot)
#' plotEnrichment2(pathway = gs[["bottom"]], stats = stats, statsAnnot = stannot)
#' plotEnrichment2(pathway = gs[["r1"]], stats = stats, statsAnnot = stannot)
#' plotEnrichment2(pathway = gs[["r2"]], stats = stats, statsAnnot = stannot)
#' 
#' plotGseaTable2(pathways = gs, stats = stats, statsAnnot=stannot, fgseaRes = res)
plotGseaTable2 <- function(pathways, stats, fgseaRes,statsAnnot=NULL,
                          gseaParam=1,
                          colwidths=c(5, 3, 0.8, 1.2, 1.2)) {

    rnk <- rank(-stats)
    ord <- order(rnk)

    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
    statsAdj <- statsAdj / max(abs(statsAdj))

    ### ===============================================
    # pathways <- lapply(pathways, function(p) {
    #     unname(as.vector(na.omit(match(p, names(statsAdj)))))
    # })
    if (is.null(statsAnnot))
        statsAnnot <- names(stats)
    statsAnnot <- statsAnnot[ord]
    pathways <- lapply(pathways, function(p) {
        listmatch(p, statsAnnot)
    })
    ### ===============================================

    ps <- lapply(names(pathways), function(pn) {
        p <- pathways[[pn]]
        annotation <- fgseaRes[match(pn, fgseaRes$pathway), ]
        list(
            textGrob(pn, just="right", x=unit(0.95, "npc")),
            ggplot() +
                geom_segment(aes(x=p, xend=p,
                                 y=0, yend=statsAdj[p]),
                             size=0.2) +
                scale_x_continuous(limits=c(0, length(statsAdj)),
                                   expand=c(0, 0)) +
                scale_y_continuous(limits=c(-1, 1),
                                   expand=c(0, 0)) +
                xlab(NULL) + ylab(NULL) +
                theme(panel.background = element_blank(),
                      axis.line=element_blank(),
                      axis.text=element_blank(),
                      axis.ticks=element_blank(),
                      panel.grid = element_blank(),
                      axis.title=element_blank(),
                      plot.margin = rep(unit(0,"null"),4),
                      panel.spacing = rep(unit(0,"null"),4)
                ),
            textGrob(sprintf("%.2f", annotation$NES)),
            textGrob(sprintf("%.1e", annotation$pval)),
            textGrob(sprintf("%.1e", annotation$padj))
            )
    })


    rankPlot <-
        ggplot() +
        geom_blank() +
        scale_x_continuous(limits=c(0, length(statsAdj)),
                           expand=c(0, 0)) +
        scale_y_continuous(limits=c(-1, 1),
                           expand=c(0, 0)) +
        xlab(NULL) + ylab(NULL) +
        theme(panel.background = element_blank(),
              axis.line=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid = element_blank(),
              axis.title=element_blank(),
              plot.margin = unit(c(0,0,0.5,0), "npc"),
              panel.spacing = unit(c(0,0,0,0), "npc")
        )

    grobs <- c(
        lapply(c("Pathway", "Gene ranks", "NES", "pval", "padj"), textGrob),
        unlist(ps, recursive = FALSE),
        list(nullGrob(),
             rankPlot,
             nullGrob(),
             nullGrob(),
             nullGrob()))

    # not drawing column if corresponding colwidth is set to zero
    grobsToDraw <- rep(colwidths != 0, length(grobs)/length(colwidths))


    grid.arrange(grobs=grobs[grobsToDraw],
                 ncol=sum(colwidths != 0),
                 widths=colwidths[colwidths != 0])
}
