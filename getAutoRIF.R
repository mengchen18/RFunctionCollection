getAutoRIF <- function(term, filter = TRUE) {
  # https://amp.pharm.mssm.edu/geneshot/
  term <- gsub(" ", "%20", term)
  term <- paste(term, collapse = ",")
  qry <- sprintf("http://amp.pharm.mssm.edu/geneshot/api/search/auto/%s", term)
  r <- jsonlite::read_json(qry)
  v <- sapply(r$gene_count, unlist)
  df <- data.frame(
    gene = colnames(v),
    n = v[1, ],
    perc = v[2, ],
    stringsAsFactors = FALSE
  )
  df$rank <- df$n * df$perc
  if (filter)
    df <- df[df$n > min(ceiling(nrow(df)/200), 3), ]
  attr(df, "term") <- r$search_term
  attr(df, "pubmedID_count") <- r$pubmedID_count
  df
}
