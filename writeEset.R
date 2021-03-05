writeEset <- function(eset, file, creator) {
  
  td <- function(tab) {
    ic <- which(sapply(tab, is.list))
    if (length(ic) > 0) {
      for (ii in ic) {
        tab[, ii] <- sapply(tab[, ii], paste, collapse = ";")
      }
    }
    tab
  }
  
  expr <- exprs(eset)
  pd <- pData(eset)
  fd <- fData(eset)
  
  wb <- createWorkbook(creator = creator)
  addWorksheet(wb, sheetName = "Phenotype info")
  addWorksheet(wb, sheetName = "Feature info")
  addWorksheet(wb, sheetName = "Expression")
  id <- paste0("ID", 1:nrow(expr))
  writeData(wb, sheet = "Expression", data.frame(ID = id, expr))
  writeData(wb, sheet = "Feature info", td(cbind(ID = id, fd)))
  writeData(wb, sheet = "Phenotype info", td(pd))
  saveWorkbook(wb, file = file, overwrite = TRUE)
}
