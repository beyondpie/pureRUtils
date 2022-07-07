#' Convert snap object to Seurat by gmat
#'
#' @param snap snap object defined in SnapATAC package
#' snap@gmat must be not empty; snap@smat@dmat should not be empty
#' @param eigDims vector, used for choosing PCA components, default 1:50
#' @param assay characters, name used in Seurat object
#' @param pcaPrefix characters, default "SnapATAC_"
#' @param nameDelim characters, default "-"
#' @param useSnapATACEmbed bool, use SnapATAC embedding as pca
#' embeding in Seurat, default TRUE.
#' @return Seurat object
#' @import Seurat
#' @export
snapGmat2Seurat <- function(snap, eigDims = 1:50,
                            assay = "GeneScore",
                            pcaPrefix = "SnapATAC_",
                            nameDelim = "-",
                            useSnapATACEmbed = TRUE) {
  # check snap@gmat
  # check snap@smat@dmat
  metaData <- snap@metaData
  rownames(metaData) <- paste(metaData$sample, metaData$barcode, sep = ".")
  gmatUse <- Matrix::t(snap@gmat)
  colnames(gmatUse) <- paste(metaData$sample, metaData$barcode, sep = ".")
  snapSeurat <- Seurat::CreateSeuratObject(counts = gmatUse, assay =
                                                               assay,
                                           names.delim = nameDelim)
  snapSeurat <- Seurat::AddMetaData(snapSeurat, metadata = metaData)
  if (useSnapATACEmbed) {
    message("Use SnapATAC Embed as pca in Seurat.")
    pcaUse <- snap@smat@dmat[, eigDims]
    rownames(pcaUse) <- paste(metaData$sample, metaData$barcode, sep = ".")
    colnames(pcaUse) <- paste0(pcaPrefix, 1:ncol(pcaUse))
    snapSeurat[["pca"]] <- methods::new(Class = "DimReduc", cell.embeddings = pcaUse,
      feature.loadings = matrix(0, 0, 0),
      feature.loadings.projected = matrix(0, 0, 0),
      assay.used = assay, stdev = rep(1, ncol(pcaUse)),
      key = pcaPrefix,
      jackstraw = new(Class = "JackStrawData"), misc = list())
  }
  return(snapSeurat)
}

#' Consensus plot.
#' @param consensusFile characters
#' @param type characters, cell type name for title
#' @param outPDF characters
#' @param width numeric, default is 7
#' @param height numeric, default is 7
#' @return None.
#' Side effect: generate the figure file.
#' @importFrom graphics axis legend mtext par
#' @importFrom methods new
#' @importFrom utils read.table
#' @importFrom stats as.formula complete.cases
#' @importFrom grDevices dev.off pdf
#' @export
plotConsensus <- function(consensusFile, type,
                          outPDF, width = 7, height = 7) {
  consensus <- read.table(consensusFile)
  colnames(consensus) <- c("file", "Resolution", "Dispersion",
                           "ProportionOfAmbiguousClustering")
  pdf(file = outPDF, width = width, height = height)
  par(mar = c(5, 4, 4, 4) + 0.3)
  plot(consensus$Resolution, consensus$Dispersion,
       type = "l", pch = 16, col = "red", lwd = 4,
       xlab = "Resolution", ylab = "Dispersion", cex.lab = 1.5,
       main =
         paste("Consensus analysis for Leiden-base clustering on",
               type, "cells")
       )
  par(new = TRUE)
  plot(consensus$Resolution, consensus$ProportionOfAmbiguousClustering,
       pch = 17, col = "blue",
       axes = FALSE, type = "l", lwd = 4, xlab = "", ylab = ""
       )
  axis(side = 4,
       at = pretty(range(consensus$ProportionOfAmbiguousClustering)))
  legend("topright", legend = c("Dispersion (left;red)", "PAC (right;blue)"),
         col = c("red", "blue"), cex = 1, lty = c(3, 3))
  mtext("Proportion Of Ambiguous Clustering (PAC)",
        side = 4, line = 3, cex = 1.5)
  dev.off()
}


#' Down sample snap object based on a given cluster meta.
#' @param snap SnapObject from SnapATAC
#' @param cluster vector of characters/integer
#' aligned with the cells in snap in order.
#' @param n integer, number of samples in each unique cluster id.
#' @return SnapObject
#' @export
downSampleOnSnap <- function(snap, cluster, n = 200) {
  snapList <- lapply(unique(cluster), function(i) {
    s <- snap[cluster %in% i, , drop = FALSE]
    if(SnapATAC::nrow(s) <= n) {
      return(s)
    } else {
      return(s[sample(SnapATAC::nrow(s), size = n, replace = F), ])
    }
  })
  return(SnapATAC::snapListRbind(snapList = snapList))
}

#' SnapATAC landmark embedding.
#' @param snap snap object of SnapATAC, default NULL
#' @param snapFile characters name of the snap File to load, default NULL
#' @param blacklistGR GenomicRanges object, default NULL
#' @param blackListFile characters file of the blacklistGR, default NULL
#' @param binSize integer dim of bmat, default 5,000
#' @param nPC integer dim of embeding for embedding, default 50
#' @param ncores integer number of cores to use for parallel, default 1
#' @param outFile characters to which we write snap with embedding, default NULL
#' @param excludeChr characters used to filter bmat, default "random|chrM"
#' @param removeBmat bool default FALSE
#' @param removeJmat bool default FALSE
#' @return SnapObject
#' @export
landmarkEmbedding <- function(snap = NULL,
                              snapFile = NULL,
                              blacklistGR = NULL,
                              blackListFile = NULL,
                              binSize = 5000,
                              nPC = 50,
                              ncores = 1,
                              outFile = NULL,
                              excludeChr = "random|chrM",
                              removeBmat = FALSE,
                              removeJmat = FALSE) {
  if (!is.null(snapFile)) {
    if (grepl("\\.RData", snapFile)) {
      snap <- loadRData(snapFile, check = FALSE)
    } else {
      snap <- readRDS(snapFile)
    }
  }
  if (is.null(snap)) {
    stop("No snap found.")
  }
  if (!is.null(blackListFile)) {
    blackList <- read.table(blackListFile,
      header = FALSE, quote = "")
    blacklistGR <- GenomicRanges::GRanges(                                                                                                                                                      
      seqnames = blackList[, 1],                                                                                                                                                 
      ranges = IRanges::IRanges(blackList[, 2], blackList[, 3])                                                                                                                           
    )
  }
  if(is.null(blacklistGR)) {
    stop("No blacklist found")
  }
  if (!is.null(outFile)) {
    prepareOutfile(outFile)
  }
  ## add bmat
  if (nrow(snap@bmat) > 0) {
    message("Snap has bmat, will skip addBmat.")
  } else {
    snap <- SnapATAC::addBmatToSnap(
      obj = snap, bin.size = binSize, do.par = TRUE,
      num.cores = ncores,
      checkSnap = FALSE
    )
  }
  ## preprocess
  idy <- S4Vectors::queryHits(
    GenomicRanges::findOverlaps(snap@feature, blacklistGR))
  if (length(idy) > 0) {
    message("remove blacklist")
    snap <- snap[, -idy, mat = "bmat"]
  }

  chr.exclude <-  GenomeInfoDb::seqlevels(snap@feature)[
    grep(excludeChr, GenomeInfoDb::seqlevels(snap@feature))
  ]

  idy <- grep(paste(chr.exclude, collapse = "|"), snap@feature)
  if (length(idy) > 0) {
    message("remove ", excludeChr)
    snap <- snap[, -idy, mat = "bmat"]
  }

  ## filter based on bincoverage
  bin.cov <- log10(Matrix::colSums(snap@bmat) + 1)
  binCutoff <- quantile(bin.cov[bin.cov > 0], 0.95)
  idy <- which(bin.cov <= binCutoff & bin.cov > 0)
  snap <- snap[, idy, mat = "bmat"]
  ## binary snap bmat
  if(max(snap@bmat) == 1) {
    message("@bmat is already binarized.")
  } else {
    snap <- SnapATAC::makeBinary(snap, mat = "bmat")
  }

  snap <- SnapATAC::runDiffusionMaps(
    obj = snap,
    input.mat = "bmat",
    num.eigs = nPC,
    method = "RSpectra"
  )
  snap@metaData$landmark <- 1
  if(removeBmat & (nrow(methods::slot(snap, "bmat")) != 0L)) {
    message("Remove bmat.")
    snap <- SnapATAC::rmBmatFromSnap(snap)
  }
  if(removeJmat) {
    message("Remove jmat.")
    snap@jmat <- SnapATAC::newJaccard()
  }
  
  if (!is.null(outFile)) {
    saveRDS(object = snap, file = outFile)
  }
  return(snap)
}

#' SnapATAC query embedding.
#' @param snapLandmark snap object of SnapATAC, default NULL
#' @param snaplandmarkFile characters name of the snap File to load, default NULL
#' @param snapQuery snap object of SnapATAC, default NULL
#' @param snapQueryFile characters name of the snap File to load, default NULL
#' @param binSize integer dim of bmat, default 5,000
#' @param outFile characters to which we write snap with embedding, default NULL
#' @param chunkSize integer default 20,000
#' @param removeBmat bool default FALSE
#' @param removeJmat bool default FALSE
#' @return SnapObject
#' @export
queryEmbedding <- function(snapLandmark = NULL,
                           snapLandmarkFile = NULL,
                           snapQuery = NULL,
                           snapQueryFile = NULL,
                           binSize = 5000,
                           outFile = NULL,
                           chunkSize = 20000,
                           removeBmat = TRUE,
                           removeJmat = TRUE) {
  if (!is.null(snapLandmarkFile)) {
    snapLandmark <- readRDS(snapLandmarkFile)
  }
  if (is.null(snapLandmark)) {
    stop("No snapLandmark found.")
  }
  if (!is.null(snapQueryFile)) {
    snapQuery <- readRDS(snapQueryFile)
  }
  if (is.null(snapQuery)) {
    stop("No snapQuery found.")
  }
  if (!is.null(outFile)) {
    prepareOutfile(outFile)
  }
  n <- SnapATAC::nrow(snapQuery)
  if (n > chunkSize) {
    message(n, " has too many cells than ",
            chunkSize, " chunk.")
    message("now partition them into chunks.")
    chunks <- split(seq(n), ceiling(seq(n) / chunkSize))
  } else {
    chunks <- NULL
  }
  message("Add bmat.")
  if (is.null(chunks)) {
    snapQuery <- queryEmbedding.default(
      snapLandmark, snapQuery, binSize, removeBmat, removeJmat)
  } else {
    snapQueryList <- lapply(chunks, function(i) {
      s <- snapQuery[i, , drop = FALSE]
      s <- queryEmbedding.default(
        snapLandmark, s, binSzie, removeBmat, removeJmat
      )
      return(s)
    })
    snapQuery <- SnapATAC::snapListRbind(snapList = snapQueryList,
                                         checkSnap = FALSE)
  }
  if(!is.null(outFile)) {
    saveRDS(object = snapQuery, file = outFile)
  }
  return(snapQuery)
}

queryEmbedding.default <- function(snapLandmark,
                                   snapQuery,
                                   binSize = 5000,
                                   removeBmat = TRUE,
                                   removeJmat = TRUE) {
  if (nrow(snapQuery@bmat) > 0) {
    message("snapQuery has bmat, will skip addBmat for it.")
  } else {
    snapQuery <- SnapATAC::addBmatToSnap(snapQuery, bin.size = binSize)
  }
  if(max(snapQuery@bmat) == 1) {
    message("@bmat is already binaryized.")
  } else {
    snapQuery <- SnapATAC::makeBinary(snapQuery)
  }
  idy <- unique(S4Vectors::queryHits(
    GenomicRanges::findOverlaps(snapQuery@feature, snapLandmark@feature)))
  snapQuery <- snapQuery[, idy, mat = "bmat"]
  message("Run embedding for query set.")
  snapQuery <- SnapATAC::runDiffusionMapsExtension(
    obj1 = snapLandmark,
    obj2 = snapQuery,
    input.mat = "bmat"
  )
  snapQuery@metaData$landmark <- 0
  if (removeBmat & (nrow(methods::slot(snapQuery, "bmat")) != 0L)) {
    message("Remove query bmat")
    snapQuery <- SnapATAC::rmBmatFromSnap(snapQuery)
  }
  if (removeJmat) {
    message("Remove jmat.")
    snapQuery@jmat <- SnapATAC::newJaccard()
  }
  return(snapQuery)
}

#' Run SnapATAC KNN.
#'
#' It will check snapLandmark and snapQuery firstly,
#' if snapLandmark is not NULL, will use it.
#' if snapQuery then is also NULL, will merge them.
#' if both of them are NULL, then use snapAllFile, snapAll in order.
#' 
#' @param snapAll SnapObject, default NULL
#' @param snapAllFile characters, default NULL
#' @param snapLandmark SnapObject, default NULL
#' @param snapLandmarkFile characters, default NULL
#' @param snapQuery SnapObject, default NULL
#' @param snapQueryFile SnapObject, default NULL
#' @param removeBmat bool, default TRUE
#' @param removeJmat bool, default TRUE
#' @param k integer, default 50
#' @param dims integer or vector, default is 1:30
#' @param method characters, method for KNN (either "RANN" or "Annoy")
#' default "RANN"
#' @param runUMAP bool, default TRUE
#' @param umapNcores integer, default 1
#' @param outmmtxFile characters, output file of mmtx, default NULL
#' @param outSnapFile characters, default NULL
#' @return SnapObject
#' @export
runKNN <- function(snapAll = NULL,
                   snapAllFile = NULL,
                   snapLandmark = NULL,
                   snapLandmarkFile = NULL,
                   snapQuery = NULL,
                   snapQueryFile = NULL,
                   removeBmat = TRUE,
                   removeJmat = TRUE,
                   k = 50,
                   dims = 1:30,
                   method = "RANN",
                   runUMAP = TRUE,
                   umapNcores = 1,
                   outmmtxFile = NULL,
                   outSnapFile = NULL) {
  if (!is.null(snapQueryFile)) {
    snapQuery <- readRDS(snapQueryFile)
  }
  if (!is.null(snapLandmarkFile)) {
    snapLandmark <- readRDS(snapLandmarkFile)
  }

  if (!is.null(snapLandmark)) {
    if (!is.null(snapQuery)) {
      snapAll <- SnapATAC::snapRbind(snapLandmark, snapQuery)
    } else {
      snapAll <- snapLandmark
    }
  }
  if (!is.null(snapAllFile)) {
    snapAll <- readRDS(snapAllFile)
  }

  if (is.null(snapAll)) {
    stop("Final snapAll is NULL.")
  }

  if (!is.null(outmmtxFile)) {
    prepareOutfile(outmmtxFile)
  }
  if (!is.null(outSnapFile)) {
    prepareOutfile(outSnapFile)
  }

  if (removeBmat & (nrow(methods::slot(snapAll, "bmat")) != 0L)) {
    message("Remove snapAll bmat.")
    snapAll <- SnapATAC::rmBmatFromSnap(snapAll)
  }

  if (removeJmat) {
    message("Remove snapAll jmat.")
    snapAll@jmat <- SnapATAC::newJaccard()
  }
  if (length(dims) == 1) {
    message("Dim is a scalar and convert it to vector.")
    dims <- 1:dims
  }

  message("Run KNN with k ", k, " and size of dims ", length(dims))
  message("KNN method: ", method)
  snapAll <- SnapATAC::runKNN(obj = snapAll, eigs.dims = dims, k = k, method = method)
  if (!is.null(outmmtxFile)) {
    message("save graph as mmtx file.")
    Matrix::writeMM(snapAll@graph@mat, file = outmmtxFile)
  }
  if (runUMAP) {
    message("Run Umap.")
    snapAll <- SnapATAC::runViz(
      obj = snapAll, tmp.folder = tempdir(), dims = 2,
      eigs.dims = dims, method = "uwot", seed.use = 2022,
      num.cores = umapNcores
    )
  }
  if (!is.null(outSnapFile)) {
    saveRDS(snapAll, outSnapFile)
  }
  return(snapAll)
}

#' Run Leiden Algorithm.
#' @param snap SnapObject default NULL
#' @param snapFile characters default NULL
#' @param r numeric resolution paramter for leiden, default 0.5
#' @param pt characters partition type for leiden, default "RB"
#' @param seed integer default NULL
#' @param pathToPython characters where the python is, default NULL
#' @param outLeidenFile characters default NULL
#' @param outClusterMetaCSV characters default NULL
#' @param outClusterPDF characters default NULL
#' @param pdfn function used to draw figures for outClusterPDF
#' the inputs are snap object and others
#' default NULL.
#' @param colName characters used for record cluster namae in snap meta,
#' default is "cluster". if snap@metaData has this column,
#' will use [colName].1 instead
#' @param dims integer or vector for UMAP default is 1:30
#' @param umapNcores integer default 1
#' @param ... use for pdfn function
#' @return SnapObject, slot "cluster" is the result
#' @import reticulate
#' @export
runLeiden <- function(snap = NULL,
                      snapFile = NULL,
                      r = 0.5,
                      pt = "RB",
                      seed = 10,
                      pathToPython = NULL,
                      outLeidenFile = NULL,
                      outClusterMetaCSV = NULL,
                      outClusterPDF = NULL,
                      pdfn = NULL,
                      colName = "cluster",
                      dims = 1:30,
                      umapNcores = 1,
                      ...) {
  if (!is.null(snapFile)) {
    snap <- readRDS(snapFile)
  }
  if (is.null(snap)) {
    stop("Snap is NULL.")
  }
  invisible(lapply(c(outLeidenFile, outClusterMetaCSV, outClusterPDF),
    function(i) {if (!is.null(i)) {prepareOutfile(i)}}))
  message("Run Leiden algorithm with: resolution ", r,
          " partition type: ", pt)
  if(!is.null(pathToPython)) {
    reticulate::use_python(pathToPython, required = TRUE)
    message("python path: ", pathToPython)
  }
  setSessionTimeLimit(cpu = Inf, elapsed = Inf)
  message("Loading python smmuty package.")
  pymod <- reticulate::import(module = "smmuty", convert = FALSE)
  message("Clustering.")
  c <- as.factor(reticulate::py_to_r(
    pymod$leiden(knn = reticulate::r_to_py(snap@graph@mat),
                 reso = r,
                 seed = 10L,
                 opt = pt)
  ))
  ## snap only accept factors
  snap@cluster <- c
  ## Label start from 1
  c <- as.integer(c)
  message("Summarize the clustering result:")
  print(table(c))
  m <- snap@metaData
  if (!("CellID" %in% colnames(m))) {
    m$CellID <- paste(snap@sample, snap@barcode, sep = ".")
  }
  if (colName %in% colnames(m)) {
    warning(colName, " is in the column of snap meta data.")
    colName <- paste0(colName, ".1")
    warning("Use ", colName, " instead.")
  }
  m[, colName] <- c

  ## cluster start from 1
  if (!is.null(outLeidenFile)) {
    message("Output Leiden result: ", outLeidenFile)
    write.table(x = c, file = outLeidenFile,
      row.names = FALSE, col.names = FALSE, quote = FALSE)
  } else {
    warning("No output of Leiden cluster file.")
  }
  if (!is.null(outClusterPDF)) {
    if (nrow(methods::slot(snap, "umap")) == 0L) {
      warning("No UMAP found, and now calculate it.")
      if (length(dims) == 1) {
        message("Dim is a scalar and convert it to vector.")
        dims <- 1:dims
      }
      snap <- SnapATAC::runViz(
        obj = snap, tmp.folder = tempdir(), dims = 2,
        eigs.dims = dims, method = "uwot", seed.use = 2022,
        num.cores = umapNcores
      )
    } # end of update umap
    if (is.null(pdfn)) {
      warning("No pdfn is found.")
      warning("No output of custer pdf.")
    } else {
      message("Output cluster pdf: ", outClusterPDF)
      withr::with_pdf(outClusterPDF, code = {
        pList <- pdfn(embed = snap@umap,
             meta = m,
             checkRowName = FALSE,
             names = c(colName, "TSSe", "log10UMI"),
             discretes = c(TRUE, FALSE, FALSE),
             legends = c(FALSE, TRUE, TRUE),
             addLabels = c(TRUE, FALSE, FALSE),
             ...)
        lapply(pList, function(p) {
          if(!is.null(p)) {
            print(p)
          }
        }) }, width = 10, height = 10)
    }
  } else {
    warning("No output of custer pdf.")
  }# end of outClusterPDF
  if (!is.null(outClusterMetaCSV)) {
    message("Output meta data to: ", outClusterMetaCSV)
    if (nrow(methods::slot(snap, "umap")) > 1L) {
      m$UMAP1 <- snap@umap[, 1]
      m$UMAP2 <- snap@umap[, 2]
    }
    write.table(x = m, file = outClusterMetaCSV,
      row.names = FALSE, col.names = TRUE,
      sep = "\t", quote = FALSE)
  } else {
    warning("No output of cluster meta file.")
  }
  return(snap)
}

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
