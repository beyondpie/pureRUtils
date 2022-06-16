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
  legend("topright", legend = c("Dispersion (left)", "PAC (right)"),
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
    blacklistGR <- read.table(blackListFile,
      header = FALSE, quote = "")
  }
  if(is.null(blacklistGR)) {
    stop("No blacklist found")
  }
  if (!is.null(outFile)) {
    prepareOutfile(outFile)
  }
  ## add bmat
  snap <- SnapATAC::addBmatToSnap(
    obj = snap, bin.size = binSize, do.par = TRUE,
    num.cores = ncores,
    checkSnap = FALSE
  )
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
  snap <- SnapATAC::makeBinary(snap, mat = "bmat")

  snap <- SnapATAC::runDiffusionMaps(
    obj = snap,
    input.mat = "bmat",
    num.eigs = nPC,
    method = "RSpectra"
  )
  snap@metaData$landmark <- 1
  if(removeBmat) {
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
  message("Add bmat.")
  snapQuery <- SnapATAC::addBmatToSnap(snapQuery, bin.size = binSize)
  snapQuery <- SnapATAC::makeBinary(snapQuery)
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
  if(removeBmat) {
    message("Remove query bmat")
    snapQuery <- SnapATAC::rmBmatFromSnap(snapQuery)
  }
  if(removeJmat) {
    message("Remove jmat.")
    snapQuery@jmat <- SnapATAC::newJaccard()
  }
  if(!is.null(outFile)) {
    saveRDS(object = snapQuery, file = outFile)
  }
  return(snapQuery)
}
