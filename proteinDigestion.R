# Digest proteins 
trypsinDigest <- function(x, maxMissed = 2, minLen = 5, maxLen = 52, annot = NA) {
  si <- str_locate_all(x, pattern = c("R|K", "RP|KP"))
  rk <- unique(c(0, setdiff(si[[1]][, 1], si[[2]][, 1]), nchar(x)))
  
  seqs <- sapply(1:(length(rk)-1), function(i) {
    str_sub(x, rk[i]+1, rk[(i+1):(i+1+maxMissed)])
  })
  ss <- na.omit(c(seqs))
  ss <- ss[nchar(ss) <= maxLen & nchar(ss) >= minLen]
  if (length(ss) < 1)
    return(NULL)
  if (is.na(annot) && !is.null(attr(x, "Annot")))
    annot <- attr(x, "Annot")
  data.frame(seq = ss, annot = annot, stringsAsFactors = FALSE)
}
# Digest proteins list
trypsinDigestList <- function(x, ..., mc.cores) {
  qq <- parallel::mclapply(x, trypsinDigest, mc.cores = mc.cores, ...)
  do.call(rbind, qq)
}
