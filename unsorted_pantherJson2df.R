
pantherJson2df <- function(js) {

  getfixwidth <- function(x) {
    if (is.null(x))
      x <- NA
    x
  }
  
  einfo <- function(v, cat) {
      c(GO_group = cat,
        GO_level = getfixwidth(v$term$level),
        GO_ID = getfixwidth(v$term$id),
        GO_Desc = getfixwidth(v$term$label),
        N_reference = getfixwidth(v$number_in_reference),
        N_observed = getfixwidth(v$input_list$number_in_list),
        N_expected = getfixwidth(v$input_list$expected),
        Fold_enrichment = getfixwidth(v$input_list$fold_enrichment),
        Direction = getfixwidth(v$input_list$plus_minus),
        Q_value = getfixwidth(v$input_list$pValue),
        Genes = getfixwidth(paste(unlist(v$input_list$mapped_id_list), collapse = ";"))
        )
  }
  
  ii <- js$overrepresentation$group
  lt <- lapply(seq_along(ii), function(i) {
    x <- ii[[i]][[1]]
    if (!is.null(names(x)))
      r <- einfo(x, cat = i) else {
        r <- sapply(x, einfo, cat = i)
        r <- t(r)
      }
    r
  })
  
  sapply(lt, dim)
  tab <- do.call(rbind, lt)
  tab <- as.data.frame(tab, stringsAsFactors = FALSE)
  numcol <- c("GO_group", "GO_level", "N_reference", "N_observed", "N_expected",
    "Fold_enrichment", "Q_value")
  tab[numcol] <- lapply(tab[numcol], as.numeric)
  
  attr(t, "tool_release_date") <- js$overrepresentation$tool_release_date
  attr(t, "data_version_release_date") <- js$overrepresentation$data_version_release_date
  attr(t, "test_type") <- js$overrepresentation$test_type
  attr(t, "correction") <- js$overrepresentation$correction
  attr(t, "annotation_type") <- js$overrepresentation$annotation_type
  tab
}

