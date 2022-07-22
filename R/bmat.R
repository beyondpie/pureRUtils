#' Shuffle bmat in both row and column.
#' @param bmat sparseMatrix, cell by bins
#' @return sparseMatrix
#' @export
shuffleBmat <- function(bmat) {
  cnt <- sum(bmat)
  nrow <- nrow(bmat)
  ncol <- ncol(bmat)
  message("Counts from bmat: ", cnt)
  message("bmat dim: ", nrow, " x ", ncol)
  r <- tryCatch(
    expr = {
      i <- sample(length(bmat), size = cnt)
      jth <- i %% ncol
      jth[jth == 0] <- ncol
      ith <- ceiling(i / ncol)
      Matrix::sparseMatrix(i = ith, j = jth, x = 1)
    }, error = function(cond) {
      message("ShuffleBmat faces errors.")
      message(cond)
      message("Use sample row and column way.")
      ith <- sample(seq(nrow), size = cnt, replace = TRUE)
      jth <- sample(seq(ncol), size = cnt, replace = TRUE)
      r <- Matrix::sparseMatrix(i = ith, j = jth, x = 1)
      return(r)
    }, finally = {
      message("SuffleBmat is done.")
    })
  return(r)
}

#' Get enriched bin's percentage cutoff based on shuffledBmat.
#' @param sBmat sparseMatrix, cell by bins, randomly shuffled.
#' @param scale double, default is 3.0, scale the standard deviation.
#' @param noLessThan double, default is 0.04.
#' cutoff should be no less than this.
#' @return double
#' @export
getCutoffFromShuffledBmat <- function(sBmat, scale = 3.0,
                                     noLessThan = 0.04) {
  perct <- Matrix::colSums(sBmat) / nrow(sBmat)
  m <- mean(perct)
  std <- sd(perct)
  message("Mean of bin percentage: ", round(m, 5))
  message("Std of bin percentage: ", round(std, 5))
  r <- m + std * scale
  r <- max(noLessThan, r)
  message("Cutoff: ", round(r, 5))
  return(r)
}
