#' Get the field-th column of each string splitted by given symbols.
#' 
#' @param strs, 1d vector of characters
#' @param split characters, default is "_"
#' @param fixed bool, default is TRUE
#' @param field integer, default is 1
#' @return 1d vector of characters
#' @export
splitSelect <- function(strs,
                        split = "_",
                        fixed = TRUE,
                        field = 1) {
  t <- strs[1]
  if(length(base::strsplit(t, split = split,
                           fixed = fixed)[[1]]) < field){
    stop("Cannot get the filed: ", field, " from ", t,
         " with split ", "split")
  }
  d <- vapply(strs, function(s) {
    return(base::strsplit(x = s, split = split,
                          fixed = fixed)[[1]][field])
  }, FUN.VALUE = "3C")
  return(d)
}
