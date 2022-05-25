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

