#' Get Seurat dot plot for snap gmat.
#' @param snap snap object, default is NULL
#' @param snapList list of snap object, default is NULL
#' @param cluster vector of characters, integers, or factors
#' must have names, which are the cell ids. Cluster size can be
#' larger than snap or snapList size, but have to make sure all
#' the cells in snap or snapList have the cluster identification.
#' @param features vector of characters, genes for plotting.
#' @param group.by characters default is cluster
#' @param cellidSep characters default is "."
#' @param axis.text.x.angle numeric, default is 30
#' @param axis.text.x.hjust numeric, default is 1
#' @param ... used for Seurat::DotPlot function
#' @return ggplot object
#' @export
plotDotplot <- function(snap = NULL,
                        snapList = NULL,
                        cluster = NULL,
                        features,
                        groupBy = "cluster",
                        cellidSep = ".",
                        axis.text.x.angle = 30,
                        axis.text.x.hjust = 1,
                        ...) {
  if (!is.null(snapList)) {
    message("Use snapList instead of a single snap object.")
    snap <- SnapATAC::snapListRbind(snapList)
  }
  if (is.null(snap)) {
    stop("snap object is NULL")
  }
  if ((!is.null(cluster)) & (is.null(names(cluster)))) {
    stop("cluster have no names.")
  }
  snapSeurat <- snapGmat2Seurat(
    snap = snap, eigDims = 1:50,
    assay = "GeneScore", pcaPrefix = "SnapATAC_",
    useSnapATACEmbed = FALSE
  )
  snapSeurat <- Seurat::NormalizeData(
    object = snapSeurat,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    margin  = 1)
  if (is.null(cluster) & (!is.null(snapList))) {
    message("cluster is NULL, but snapList is not NULL.")
    message("Then we treat each snap object as one cluster")
    if (!is.null(names(snapList)))  {
      message("Label the cluster id as the names of snapList ",
              "since they are not NULL")
    } else {
      message("Label the cluster id as the order of snapList.")
    }
    clist <- lapply(seq_along(snapList), function(i) {
      s <- snapList[[i]]
      if (is.null(names(snapList))) {
        r <- rep(names(snapList)[i], SnapATAC::nrow(s))
      } else {
        r <- rep(i, SnapATAC::nrow(s))
      }
      names(r) <- paste(s@sample, s@barcode, sep = cellidSep)
      return(r)
    })
    cluster <- unlist(clist)
  }
  if (!is.null(cluster)) {
    cellids <- with(
      snapSeurat@meta.data, paste(sample, barcode, sep = "."))
    if (sum(cellids %in% names(cluster)) != length(cellids)) {
      stop("some cells have no cluster info.")
    }
    snapSeurat[[groupBy]] <- cluster[cellids]
  }
  ## Seurat::Idents(snapSeurat) <- snapSeurat[[group.by]]
  if (any(duplicated(features))) {
    warning("features have replicate elements, and ",
            "only one of the replicate will be used.")
    features <- unique(features)
  }
  p <- Seurat::DotPlot(
    object = snapSeurat,
    assay = "GeneScore",
    features = features,
    group.by = groupBy,
    ...
    ) +
    ggplot2::theme(axis.text.x =
                     ggplot2::element_text(angle = axis.text.x.angle,
                                           hjust = axis.text.x.hjust))
  return(p)
}
