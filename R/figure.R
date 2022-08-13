#' Generate list of ggplot figures of umap/tsne on the features
#' If some feature is not in the meta, will have a NULL.
#'
#' @param embed data.frame cell by (x, y) dim, rownames are cells.
#' @param meta data.frame cell by features, rownames are cells,
#' colnames are features.
#' @param checkRowName bool, default FALSE
#' @param names character vector, features to be draw
#' Default is c("cluster", "tsse", "logumi")
#' @param discretes bool vector, default is c(T, F, F)
#' @param legends bool vector, default is c(F, T, T).
#' @param addLabels bool vectors, default is c(T, F, F).
#' @param n integer, number of down samples, default 10,000
#' @param ... parameters ggPoint, but seems no use there.
#' @return list of ggplot object
#' @export
plot2D <- function(embed,
                   meta,
                   checkRowName = FALSE,
                   names = c("cluster", "tsse", "log10UMI"),
                   discretes = c(TRUE, FALSE, FALSE),
                   legends = c(FALSE, TRUE, TRUE),
                   addLabels = c(TRUE, FALSE, FALSE),
                   n = 10000,
                   ...) {
  if(checkRowName) {
    commonRows <- base::intersect(rownames(embed), rownames(meta))
    if (is.null(commonRows) | (length(commonRows) == 0)) {
      stop("No common rows between embed and meta.")
    }
    embed <- embed[commonRows, ]
    meta <- meta[commonRows, ]
  }
  p <- lapply(seq_along(names), function(i) {
    if (!(names[i] %in% colnames(meta))) {
      warning(names[i], " is not in meta.Skip it.")
      return(NULL)
    }
    message("Drawing 2D image for ", names[i])
    rowIndex <- which(!is.na(meta[, names[i]]))
    if (length(rowIndex) > n) {
      rowIndex <- sample(x = rowIndex, size = n, replace = FALSE)
    }
    r <- ggPoint(
      x = embed[rowIndex, 1],
      y = embed[rowIndex, 2],
      color = meta[rowIndex, names[i]],
      discrete = discretes[i],
      title = names[i],
      labelMeans = addLabels[i],
      ...
    )
    return(r)
  })
  return(p)
}


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
#' @param scaleSizeMin integer, default is 4
#' @param scaleSizeMax integer, default is 6
#' @param scaleColorLow color name, default is "blue"
#' @param scaleColorMid color name, default is "white"
#' @param scaleColorHigh color name, default is "red"
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
                        scaleSizeMin = 4,
                        scaleSizeMax = 6,
                        scaleColorLow = "blue",
                        scaleColorMid = "white",
                        scaleColorHigh = "red",
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
  if ((!is.null(cluster)) & (is.null(snapSeurat[[groupBy]]))) {
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
                                           hjust = axis.text.x.hjust)) +
    ggplot2::scale_size(range = c(scaleSizeMin, scaleSizeMax)) +
    ggplot2::scale_colour_gradient2(low = scaleColorLow,
                                    mid = scaleColorMid,
                                    high = scaleColorHigh)
  return(p)
}

#' plotFeatureSingle from SnapATAC package
#' Allow the object can be a matrix of umap/tsne or a snap
#'
#' @param snap A snap object, default is NULL
#' @param featureValue Feature enrichment value for each cell.
#' Value will be normalized betweeen 0 and 1.
#' @param embedmat a matrix of umap or tsne, default is NULL
#' @param method Visulization method c("tsne", "umap").
#' @param pointSize Point size [1].
#' @param pointShape Point shape type [19].
#' @param downSample Downsample the original cells to
#' down.sample cells to ovoid large dataset [10,000].
#' @param pdfile pdf file name to save the plot [NULL].
#' @param pdfWidth Width of the graphics region in inches [7].
#' @param pdfHeight Height of the graphics region in inches [7].
#' @param quantiles Feature value outside this range will be removed [c(0.01, 0.99)]
#' @param ... Arguments passed to plot method.
#' @importFrom grDevices pdf dev.off
#' @importFrom scales alpha
#' @importFrom graphics plot text title legend
#' @importFrom plot3D scatter2D
#' @importFrom viridis viridis
#' @export
plotFeatureSingle.SnapATAC <- function(snap = NULL,
                                       featureValue,
                                       embedmat = NULL,
                                       method = c("tsne", "umap"),
                                       pointSize = 0.1,
                                       pointShape = 19,
                                       downSample = 10000,
                                       pdfile = NULL,
                                       pdfWidth = 7,
                                       pdfHeight = 7,
                                       quantiles = c(0.01, 0.99),
                                       ...) {
  if( (!is.null(snap)) & (!is.null(embedmat))) {
    stop("Both snap and embedmat are NULL.")
  }
  method <- match.arg(method)
  if (!is.null(snap)) {
    message("Use snap to get embedmat.")
    dataUse <- as.data.frame(methods::slot(snap, method))
  }
  if (!is.null(embedmat)) {
    message("Argument embedmat is not NULL, use it.")
    dataUse <- embedmat
  }
  ncell <- nrow(dataUse)
  colnames(dataUse) <- paste(method, c(1, 2), sep = "-")
  labNames <- colnames(dataUse)
  if (nrow(dataUse) != length(featureValue)) {
    stop("feature.value has different length with dataUse.")
  }

  quantilesLow <- quantile(featureValue, quantiles[1])
  quantilesHigh <- quantile(featureValue, quantiles[2])
  featureValue[featureValue > quantilesHigh] <- quantilesHigh
  featureValue[featureValue < quantilesLow] <- quantilesLow

  if (!is.null(pdfile)) {
    if (file.exists(pdfile)) {
      warning("pdf.file already exists")
      file.remove(pdfile)
    } else {
      if (!file.create(pdfile)) {
        stop("cannot create pdf.file, not a directory")
      }
      file.remove(pdfile)
    }
    pdf(pdfile, width = pdfWidth, height = pdfHeight)
  }

  downSample <- min(downSample, ncell)
  idx <- sort(sample(seq(ncell), downSample))
  dataUse <- dataUse[idx, , drop = FALSE]
  featureValue <- featureValue[idx]
  xlims <- c(-max(abs(dataUse[, 1])) * 1.05, max(abs(dataUse[, 1])) * 1.2)
  ylims <- c(-max(abs(dataUse[, 2])) * 1.05, max(abs(dataUse[, 2])) * 1.05)

  plot3D::scatter2D(
    x = dataUse[, 1],
    y = dataUse[, 2],
    colvar = featureValue,
    cex = pointSize,
    pch = pointShape,
    bty = "l",
    font.lab = 2,
    col.axis = "darkgrey",
    xlim = xlims,
    ylim = ylims,
    xlab = labNames[1],
    ylab = labNames[2],
    col = viridis(256, option = "D"),
    ...
  )
  box(bty = "l", lwd = 2)
  if (!is.null(pdfile)) {
    dev.off()
  }
}

#' plotFeatureSingle
#' Allow the object can be a matrix of umap/tsne or a snap
#'
#' @param snap A snap object, default is NULL
#' @param featureValue Feature enrichment value for each cell.
#' Value will be normalized betweeen 0 and 1.
#' @param embedmat a matrix of umap or tsne, default is NULL
#' @param method Visulization method c("tsne", "umap").
#' @param downSample Downsample the original cells to
#' down.sample cells to ovoid large dataset [10,000].
#' @param pdfile pdf file name to save the plot [NULL].
#' @param pdfWidth Width of the graphics region in inches [7].
#' @param pdfHeight Height of the graphics region in inches [7].
#' @param quantiles Feature value outside this range will be removed [c(0.01, 0.99)]
#' @param pal colors for ploting continuous numbers, default viridis::viridis(n=256)
#' @param colorTitle legend title, default is "gene score"
#' @param pointSize numeric, default is 0.1
#' @param labelSize integer default is 5
#' @param baseSize integer default is 12
#' @param legendSize integer default is 8
#' @param ... Arguments passed to ggPoint method
#' @return ggplot object
#' @importFrom viridis viridis
#' @export
plotFeatureSingle <- function(snap = NULL,
                              featureValue,
                              embedmat = NULL,
                              method = c("tsne", "umap"),
                              downSample = 10000,
                              pdfile = NULL,
                              pdfWidth = 7,
                              pdfHeight = 7,
                              quantiles = c(0.01, 0.99),
                              pal = viridis::viridis(n = 256),
                              colorTitle = "gene score",
                              pointSize = 0.1,
                              labelSize = 5,
                              baseSize = 12,
                              legendSize = 8,
                              discrete = FALSE,
                              labelMeans = FALSE,
                              ...) {
  if( (!is.null(snap)) & (!is.null(embedmat))) {
    stop("Both snap and embedmat are NULL.")
  }
  method <- match.arg(method)
  if (!is.null(snap)) {
    message("Use snap to get embedmat.")
    dataUse <- as.data.frame(methods::slot(snap, method))
  }
  if (!is.null(embedmat)) {
    message("Argument embedmat is not NULL, use it.")
    dataUse <- embedmat
  }
  ncell <- nrow(dataUse)
  colnames(dataUse) <- paste(method, c(1, 2), sep = "-")
  labNames <- colnames(dataUse)
  if (nrow(dataUse) != length(featureValue)) {
    stop("feature.value has different length with dataUse.")
  }
  if (!discrete) {
    quantilesLow <- quantile(featureValue, quantiles[1])
    quantilesHigh <- quantile(featureValue, quantiles[2])
    featureValue[featureValue > quantilesHigh] <- quantilesHigh
    featureValue[featureValue < quantilesLow] <- quantilesLow
  }

  if (!is.null(pdfile)) {
    prepareOutfile(pdfile)
    pdf(pdfile, width = pdfWidth, height = pdfHeight)
  }

  downSample <- min(downSample, ncell)
  idx <- sort(sample(seq(ncell), downSample))
  dataUse <- dataUse[idx, , drop = FALSE]
  featureValue <- featureValue[idx]
  p <- ggPoint(
    x = dataUse[, 1],
    y = dataUse[, 2],
    color = featureValue,
    discrete = discrete,
    labelMeans = labelMeans,
    pal = pal,
    colorTitle = colorTitle,
    size = pointSize,
    labelSize = labelSize,
    baseSize = baseSize,
    legendSize = legendSize,
    ...
  )
  if (!is.null(pdfile)) {
    prepareOutfile(pdfile)
    ggplot2::ggsave(filename = pdfile, plot = p,
                    width = pdfWidth, height = pdfHeight)
  }
  return(p)
}
