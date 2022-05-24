#' Sample index based on sequence depth defined by log10UMI
#' without repeat.
#' 
#' @param log10UMI numeric vector
#' @param n integer, number of samples needed
#' @param seed integer, random seed, 2022 (default)
#' @return vector of integer
#' @export
sampleBasedOnDepth2 <- function(log10UMI, n, seed = 2022) {
  depths <- log10UMI
  dens <- stats::density(x = depths, bw = "nrd", adjust = 1)
  samplingProb <- 1 / (stats::approx(x = dens$x, y = dens$y, xout = depths)$y + .Machine$double.eps)
  set.seed(seed)
  idx <- sort(sample(x = seq_along(depths),
                     size = min(n, length(depths)),
                     prob = samplingProb,
                     replace = FALSE))
  return(idx)
}
