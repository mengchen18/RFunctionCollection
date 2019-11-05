#' @title Reading .gmt file
#' @description Frequently the .gmt files are downloaed from MSigDB database
#' @param x the name/path of the gmt file to be read
#' @param id the id used in gene sets, either "SYMBOL" or "ENTREZ"

read.gmt <- function(x, id = NA) {
  
  x <- readLines(x)
  x <- strsplit(x, "\t")
  names <- sapply(x, "[", 1)
  x <- lapply(x, function(xx) {
    structure(xx[-(1:2)], 
              name = xx[1],
              link = xx[2])
  })
  names(x) <- names
  
  if (!is.na(id)) {
    id <- match.arg(id, c("SYMBOL", "ENTREZ"))
  } else {
    if (any(grepl("[A-Z]", unlist(x[1:10]))))
      id = "SYMBOL" else
        id = "ENTREZ"
  }
  attr(x, "ID") <- id
  x
} 
