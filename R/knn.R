####################
## From Seurat
####################
# Run annoy
#
#' @param data Data to build the index with
#' @param query A set of data to be queried against data
#' @param metric Distance metric; can be one of "euclidean", "cosine", "manhattan",
#' "hamming"
#' @param n.trees More trees gives higher precision when querying
#' @param k Number of neighbors
#' @param search.k During the query it will inspect up to search_k nodes which
#' gives you a run-time tradeoff between better accuracy and speed.
#' @param include.distance Include the corresponding distances
#' @param index optional index object, will be recomputed if not
#' provided
#' @return list(nn.idx, nn.dists)
#' @export
AnnoyNN <- function(data,
                    query = data,
                    metric = "euclidean",
                    n.trees = 50,
                    k,
                    search.k = -1,
                    include.distance = TRUE,
                    index = NULL) {
  idx <- index %||% AnnoyBuildIndex(
    data = data,
    metric = metric,
    n.trees = n.trees
  )
  nn <- AnnoySearch(
    index = idx,
    query = query,
    k = k,
    search.k = search.k,
    include.distance = include.distance
  )
  nn$idx <- idx
  nn$alg.info <- list(metric = metric, ndim = ncol(x = data))
  return(nn)
}

# Build the annoy index
#
#' @param data Data to build the index with
#' @param metric Distance metric; can be one of "euclidean", "cosine", "manhattan",
#' "hamming"
#' @param n.trees More trees gives higher precision when querying
#' @return RccAnny::Annoy[metric] class
#' @export
AnnoyBuildIndex <- function(data, metric = "euclidean", n.trees = 50) {
  f <- ncol(x = data)
  a <- switch(
    EXPR = metric,
    "euclidean" =  new(Class = RcppAnnoy::AnnoyEuclidean, f),
    "cosine" = new(Class = RcppAnnoy::AnnoyAngular, f),
    "manhattan" = new(Class = RcppAnnoy::AnnoyManhattan, f),
    "hamming" = new(Class = RcppAnnoy::AnnoyHamming, f),
    stop("Invalid metric")
  )
  for (ii in seq(nrow(x = data))) {
    a$addItem(ii - 1, data[ii, ])
  }
  a$build(n.trees)
  return(a)
}

#' Search an Annoy approximate nearest neighbor index
#
#' @param Annoy index, built with AnnoyBuildIndex
#' @param query A set of data to be queried against the index
#' @param k Number of neighbors
#' @param search.k During the query it will inspect up to search_k nodes which
#' gives you a run-time tradeoff between better accuracy and speed.
#' @param include.distance Include the corresponding distances in the result
#'
#' @return A list with 'nn.idx' (for each element in 'query', the index of the
#' nearest k elements in the index) and 'nn.dists' (the distances of the nearest
#' k elements)
#
#' @export
AnnoySearch <- function(index, query, k, search.k = -1, include.distance = TRUE) {
  n <- nrow(x = query)
  idx <- matrix(nrow = n, ncol = k)
  dist <- matrix(nrow = n, ncol = k)
  convert <- methods::is(index, "Rcpp_AnnoyAngular")
  if (!inherits(x = future::plan(), what = "multicore")) {
    oplan <- future::plan(strategy = "sequential")
    on.exit(future::plan(oplan), add = TRUE)
  }
  res <- future.apply::future_lapply(X = 1:n, FUN = function(x) {
    res <- index$getNNsByVectorList(query[x, ], k, search.k, include.distance)
    # Convert from Angular to Cosine distance
    if (convert) {
      res$dist <- 0.5 * (res$dist * res$dist)
    }
    list(res$item + 1, res$distance)
  })
  for (i in 1:n) {
    idx[i, ] <- res[[i]][[1]]
    if (include.distance) {
      dist[i, ] <- res[[i]][[2]]
    }
  }
  return(list(nn.idx = idx, nn.dists = dist))
}

#' OR of NULL.
#' If the first element is not NULL, then return the first one,
#' else return the second one.
#' @param lhs R object
#' @param rhs R object
#' @return lhs if lhs is not NULL, else rhs
#' @export
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}
