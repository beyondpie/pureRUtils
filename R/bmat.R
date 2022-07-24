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

#' Get number of diff bins and the pvalue given two groups.
#' @param bmat sparseMatrix, cell by bin
#' @param group vector of integer, cluster index of cells
#' @param ith integer, default 1
#' @param jth integer, default 2
#' @param threshold double, default 0.04
#' @param nperm integer, default 20
#' @return vector
#' first is number of diff bin, second is p-value
#' @importFrom stringr str_glue
#' @export
getPvalueOfNDiff.default <- function(mat,
                                     group,
                                     ith = 1,
                                     jth = 2,
                                     threshold = 0.04,
                                     nperm = 20) {
  if (ith == jth) {
    message(ith, " and ", jth, " are the same group.")
    return(c(0, 1.0))
  }
  if (length(group) != nrow(mat)) {
    stop("Group length and nrow of mat does not match.")
  }
  ## only cells within the two groups are condisered
  ## this is specific for permutation test.
  rowIndex <- group %in% c(ith, jth)
  mat <- mat[rowIndex, ]
  group <- group[rowIndex]
  ni <- sum(group == ith)
  nj <- sum(group == jth)
  si <- which(Matrix::colSums(mat[group == ith, , drop = FALSE]) / ni >= threshold)
  sj <- which(Matrix::colSums(mat[group == jth, , drop = FALSE]) / nj >= threshold)
  if (length(si) < 1) {
    message(ith, " has no enriched bins.")
    return(c(-1, 0.0))
  }
  if (length(sj) < 1) {
    message(jth, " has no enriched bins.")
    return(c(-1, 0.0))
  }
  ## if(length(si) == ncol(mat)) {
  ##   message(ith, " treats all bins as enriched.")
  ##   return(c(-1, 0.0))
  ## }
  ## if(length(sj) == ncol(mat)) {
  ##   message(jth, " treats all bins as enriched.")
  ##   return(c(-1, 0.0))
  ## }
  message("Enriched bin from ", ith, " is ", length(si))
  message("Enriched bin from ", jth, " is ", length(sj))
  if(length(intersect(si, sj)) < 1) {
    message("The two sets have no joint element.")
    return(c(length(si) + length(sj), 0.0))
  }

  nJoint <- length(intersect(si, sj))
  message(str_glue("number of joint between {ith}th and {jth}th: {nJoint}"))
  ni_j <- length(setdiff(si, sj))
  nj_i <- length(setdiff(sj, si))
  message(str_glue("Bins in {ith}th not in {jth}th: {ni_j}"))
  message(str_glue("Bins in {jth}th not in {ith}th: {nj_i}"))
  
  ## total number of diff bins
  r <- length(si) + length(sj) - 2 * length(intersect(si, sj))
  message("Total number of diff bins between ",
    ith, " ith and ", jth, " jth group is: ", r)
  ## permutation test
  ## This can be accelarated by Rcpp.
  ndiffPerm <- vapply(seq(nperm), function(i_) {
    ## randomly shuffle group
    g <- sample(group, replace = F)
    si <- which(Matrix::colSums(mat[g == ith, , drop = FALSE]) / ni >= threshold)
    sj <- which(Matrix::colSums(mat[g == jth, , drop = FALSE]) / nj >= threshold)
    s <- length(si) + length(sj) - 2 * length(intersect(si, sj))
    return(s)
  }, FUN.VALUE = 0)
  message("Quantile of the random permutation: ")
  message(paste(quantile(ndiffPerm), collapse = ";"))
  p <- sum(ndiffPerm >= r) / nperm
  return(c(r, p))
}

#' Get number of diff bins and the pvalue for any two groups
#' @param bmat sparseMatrix, cell by bin
#' @param group vector of integer, cluster index of cells
#' @param uniqueGroupWithName vector, default NULL
#' unique names for group
#' @param threshold double, default 0.04
#' @param nperm integer, default 20
#' @return list of symmetric matrix
#' "diffbin" are number of diff bin, "p" are p-values.
#' each matrix has the names defined by names of uniqueGroupWithName
#' generated by this function.
#' @export
getPvalueOfNDiff.multiGroup <- function(mat, group,
                                        uniqueGroupWithName = NULL,
                                        threshold = 0.04,
                                        nperm = 20) {
  if (is.null(uniqueGroupWithName)) {
    uniqueGroupWithName <- sort(unique(group))
    names(uniqueGroupWithName) <- paste0("c", uniqueGroupWithName)
  }
  n <- length(uniqueGroupWithName)
  s <- matrix(data = 0, nrow = n, ncol = n,
    dimnames = list(names(uniqueGroupWithName),
      names(uniqueGroupWithName)))
  r <- matrix(data = 1, nrow = n, ncol = n,
    dimnames = list(names(uniqueGroupWithName),
      names(uniqueGroupWithName)))
  if (length(uniqueGroupWithName) > 1) {
    for (i in seq(n - 1)) {
      for (j in (i + 1):n) {
        ith <- uniqueGroupWithName[i]
        jth <- uniqueGroupWithName[j]
        message("Calculating p-values for the group ", ith, " and ",
          jth, " pair.")
        t <- getPvalueOfNDiff.default(mat = mat,
          group = group,
          ith = ith,
          jth = jth,
          threshold = threshold,
          nperm = nperm)
        s[i, j] <- t[1]
        r[i, j] <- t[2]
        ## symmetric
        s[j, i] <- s[i, j]
        r[j, i] <- r[i, j]
      }
    }
  }
  return(list(diffbin = s, p = r))
}


#' Based on pval matrix, link the nodes without significant p-values.
#' @param pvalmat symmetric matrix
#' @param pthres double, default is 0.1
#' @param prefix characters, default is "c"
#' @return list of vectors, each vector/scalar is a graph module
#' defined by the pval matrix.
#' @export
getIsolatedGroup <- function(pvalmat,
                             pthres = 0.1,
                             prefix = "c") {
  adj <- 1 - (pvalmat < pthres)

  ## * if node linked with all the others, treat it as isolated one.
  ## degree <- colSums(adj)
  ## fdIds <- which(degree == (nrow(pvalmat) - 1))
  ## if(length(fdIds) > 0) {
  ##   message("Nodes with full degree: ", paste(fdIds, collapse = ","))
  ##   message("Treated as isolated nodes.")
  ##   adj[fdIds, ] <- 0.0
  ##   adj[, fdIds] <- 0.0
  ## }

  diag(adj) <- 1
  graph <- igraph::graph_from_adjacency_matrix(
    adjmatrix = adj,
    mode = "undirected",
    weighted = TRUE,
    diag = TRUE)
  modules <- igraph::components(graph = graph)
  return(igraphModule2List(modules = modules, prefix = prefix))
}

#' Transfer modules from igraph::components to a list.
#' @param modules list, returned from igraph::components
#' @param prefix characters, default is "c"
#' @return list of vector / scalar with names
#' @export
igraphModule2List <- function(modules, prefix = "c") {
  clusterIds <- seq(modules$no)
  r <- lapply(clusterIds, function(i) {
    return(which(modules$membership == i))
  })
  names(r) <- paste0(prefix, clusterIds)
  return(r)
}

#' Transer list of numeric to string
#' @param module list of numeric
#' @param sepkv characters, sep between key and value, default is ":"
#' @param sepv characters, sep between values, default is "."
#' @param sepk characters, sep between keys, default is ","
#' @return characters
#' @export
listOfNumeric2str <- function(module, sepkv = ":",
                              sepv = ".", sepk = ",") {
  if (is.null(names(module))) {
    names(module) <- paste0("c", seq_along(module))
  }
  r <- vapply(seq_along(module), function(i) {
    paste(names(module)[i], paste(module[[i]], collapse = sepv),
          sep = sepkv)
  }, "c1:1.2")
  return(paste(r, collapse = sepk))
}

#' Transfer specific strings to adjacent matrix
#' A typical string like: pval@c1-c2:0,c1-c3:0
#' @param mystring characters, like: pval@c1-c2:0,c1-c3:0
#' @param title characters, default is "pval@"
#' @param sepg characters, sep between groups, default is ","
#' @param sepk characters, sep between keys, default is "-"
#' @param sepkv characters, sep between key and value, default is ":"
#' @param clusterPrefix characters, default is "c"
#' @return Matrix symmetric
#' @export
strToAdjmat <- function(mystring, title = "pval@",
                        sepg = ",", sepk = "-", sepkv = ":",
                        clusterPrefix = "c") {
  mystring <- gsub(title, "", mystring)
  pairs <- strsplit(x = mystring, split = sepg, fixed = TRUE)[[1]]
  twoColDF <- t(vapply(pairs, function(p) {
    t <- strsplit(p, split = sepkv, fixed = TRUE)[[1]]
    return(c(t[1], trimws(t[2])))
  }, c("c1-c2", "0.3")))
  uniqueClusters <- unique(unlist(
    lapply(twoColDF[,1], function(cij) {
      t <- strsplit(cij, split = sepk, fixed = TRUE)[[1]]
      t <- gsub(clusterPrefix, "", t)
      t <- as.integer(t)
      return(t)
    })
  ))
  df <- matrix(data = 0, nrow = max(uniqueClusters),
               ncol = max(uniqueClusters))
  rownames(df) <- paste0(clusterPrefix, seq(nrow(df)))
  colnames(df) <- rownames(df)
  for(i in seq(nrow(twoColDF))) {
    cij <- twoColDF[i,1]
    v <- as.numeric(twoColDF[i,2])
    t <- strsplit(cij, split = sepk, fixed = TRUE)[[1]]
    t <- gsub(clusterPrefix, "", t)
    t <- as.integer(t)
    df[t[1], t[2]] <- v
    df[t[2], t[1]] <- v
  }
  return(df)
}


#' Transfer the upper tri of a matrix to a string
#' @param mat Matrix
#' @param sepg characters, sep between groups, default is ","
#' @param sepk characters, sep between keys, default is "-"
#' @param sepkv characters, sep between key and value, default is ":"
#' @param digit integer, default 2
#' @return characters
#' @export
uptri2str <- function(mat,
                      sepkv = ":", sepk = "-", sepg = ",",
                      digit = 2) {
  if (length(mat) == 1) {
    if (is.null(names(mat))) {
      if (is.null(rownames(mat))) {
        names(mat) <- "c1"
      } else {
        names(mat) <- rownames(mat)
      }
    }
    return(paste(names(mat), 1, sep = sepkv))
  }
  if (is.null(rownames(mat))) {
    rownames(mat) <- paste0("c", seq(nrow(mat)))
  }
  n <- nrow(mat)
  r <- vapply(seq(n - 1), function(i) {
    v <- vapply((i + 1):n, function(j) {
      paste(paste(rownames(mat)[i], rownames(mat)[j], sep = sepk),
        round(mat[i, j], digit), sep = sepkv)
    }, "c1-c2:0.1")
    paste(v, collapse = sepg)
  }, "c1-c2:0.1,c1-c3:0.2")
  paste(r, collapse = sepg)
}
