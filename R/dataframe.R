#' Select subset of rows of data.frame based on
#' the order of defined row keys.
#' Keys not in the data.frame will be ignored.
#' @param df data.frame
#' @param rowKeys vector, row keys for df
#' @param selectKeys vector, keys we want to select in order
#' @return data.frame
#' @export
getMatchDataFrame <- function(df, rowKeys, selectKeys) {
  if (length(unique(rowKeys)) != length(rowKeys)) {
    stop("Duplicates found in rowKeys.")
  }
  if (nrow(df) != length(rowKeys)) {
    stop(paste("nrow of df is ", nrow(df),
      "; while length of rowKeys is", length(rowKeys)))
  }
  if (length(unique(selectKeys)) != length(selectKeys)) {
    stop("Duplicates found in selectKeys")
  }
  message("nrow of df is ", nrow(df))
  rowIndex <- base::match(x = selectKeys, table = rowKeys, nomatch = 0)
  if (length(which(rowIndex > 0)) == 0) {
    stop("No matched rows found.")
  }
  message(paste(length(which(rowIndex > 0)), "of",
    length(selectKeys), " matched the df."))
  rowIndex <- rowIndex[rowIndex > 0]
  with(df, {
    df[rowIndex, , drop = FALSE]
  })
}
