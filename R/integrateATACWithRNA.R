#' Integration with single-cell RNA sequencing data
#'
#' @param snapSeurat Seurat object transformed from snap object
#' Assume meta.data has no column named ClusterName
#' @param rnaSeurat Seurat obejct, reference scRNA sequencing ddta
#' Assume meta.data has a column for rnaType
#' @param eigDims dims used for snapSeurat pca,
#' if NULL (default), use all the dims.
#' @param snapAssay characters, snapSeurat assay name, default "GeneScore"
#' @param preprocessSnap bool, if need to normalize and scale snap, default TRUE
#' @param preprocessRNA bool, if need to normalize, scale rnaSeurat, default TRUE
#' @param rnaTypeColnm characters, colname for rnaSeurat's type, default "ClusterName"
#' @param reso numeric, clustering resolution parameter, default 0.5
#' @return Seurat object, merge snapSeurat and rnaSeurat
#' Add predictId  and predictMaxScore from rnaSeurat for snapSeurat on metaData,
#' Add imputed gene expression from rnaSeurat for snapSeurat with Assay: "RNA"
#' Add tech as either ATAC or RNA in the meta.data
#' Clustering on the merged Seurat object as Co-Embedding Seurat object
#' @import ggplot2
#' @import Seurat
#' @export
integrateWithScRNASeq <- function(snapSeurat,
                                  rnaSeurat,
                                  eigDims = NULL,
                                  snapAssay = "GeneScore",
                                  preprocessSnap = TRUE,
                                  preprocessRNA = TRUE,
                                  rnaTypeColnm = "ClusterName",
                                  reso = 0.5) {
  if (nrow(snapSeurat) < 1) {
    stop("No cells in snapSeurat.")
  }
  if (nrow(rnaSeurat) < 1) {
    stop("No cells in rnaSeurat.")
  }
  Seurat::DefaultAssay(snapSeurat) <- snapAssay
  if (is.null(eigDims)) {
    eigDims <- 1:nrow(snapSeurat[["pca"]])
  }
  if (preprocessSnap) {
    snapSeurat <- NormalizeData(snapSeurat)
    snapSeurat <- ScaleData(snapSeurat, features = rownames(snapSeurat))
  }
  if (preprocessRNA) {
    rnaSeurat <- NormalizeData(rnaSeurat)
    rnaSeurat <- FindVariableFeatures(rnaSeurat)
    rnaSeurat <- ScaleData(rnaSeurat, features = rownames(rnaSeurat))
    rnaSeurat <- RunPCA(rnaSeurat, features = VariableFeatures(rnaSeurat))
  }
  anchors <- FindTransferAnchors(
    reference = rnaSeurat,
    query = snapSeurat,
    features = VariableFeatures(rnaSeurat),
    reference.assay = "RNA",
    query.assay = snapAssay,
    reduction = "cca"
  )
  ## predict type
  transferLabel <- TransferData(
    anchorset = anchors,
    refdata = rnaSeurat@meta.data[, rnaTypeColnm],
    weight.reduction = snapSeurat[["pca"]],
    dims = eigDims
  )
  rownames(transferLabel) <- colnames(snapSeurat)
  t <- data.frame(
    predictId = transferLabel$predicted.id,
    predictMaxScore = apply(transferLabel[, -1], 1, max),
    row.names = colnames(snapSeurat)
  )
  snapSeurat <- AddMetaData(snapSeurat, metadata = t)
  ## impute gene expression
  refdata <- GetAssayData(
    object = rnaSeurat,
    assay = "RNA",
    slot = "data"
  )
  geneImpute <- TransferData(
    anchorset = anchors,
    refdata = refdata,
    weight.reduction = snapSeurat[["pca"]],
    dims = eigDims
  )
  snapSeurat[["RNA"]] <- CreateAssayObject(counts = geneImpute@data)
  rm(geneImpute)

  ## coEmbed
  coEmbed <- merge(x = snapSeurat, y = rnaSeurat)
  DefaultAssay(coEmbed) <- "RNA"
  coEmbed$tech <- ifelse(!is.na(coEmbed@meta.data[, rnaTypeColnm]), "RNA", "ATAC")
  coEmbed <- ScaleData(coEmbed,
    features = VariableFeatures(rnaSeurat),
    do.scale = FALSE)
  coEmbed <- RunPCA(coEmbed,
    features = VariableFeatures(rnaSeurat),
    verbose = FALSE)
  coEmbed <- RunUMAP(coEmbed, dims = eigDims)
  coEmbed <- FindNeighbors(coEmbed, dims = eigDims)
  coEmbed <- FindClusters(coEmbed, resolution = reso)
  return(coEmbed)
}

#' Calculate Overlap Score matrix for ATAC and RNA jointly embedding metaTable
#' @description
#' Based on the co-embeding clustering result, we summation the minmum scores of
#' the percentages of cells on either ATAC or RNA group in each cluster.
#' @param meta data.frame,
#' three columns as ident, atacCluster, rnaCluster defined as in the following
#' parameters
#' @param ident characters, name of column ident, i.e. the cluster Ids for
#' the co-embedded seurat, "coembed.idents" as default.
#' @param atacCol characters, name of column atac,
#' "MajorType"as default
#' @param rnaCol characters, name of column rna,
#' "ClusterName" as default
#' @return data.frame of numeric, atac by rna
#' @export
getOverlapMatrix <- function(meta,
                             ident = "coembed.idents",
                             atacCol = "MajorType",
                             rnaCol = "ClusterName") {
  ident2rna <- data.frame(idents = meta[[ident]], rna_label = meta[[rnaCol]])
  ident2rna <- ident2rna[stats::complete.cases(ident2rna), ]
  ident2atac <- data.frame(idents = meta[[ident]], atac_label = meta[[atacCol]])
  ident2atac <- ident2atac[stats::complete.cases(ident2atac), ]
  rnaTable <- table(ident2rna)
  atacTable <- table(ident2atac)
  rnaPct <- apply(rnaTable, 2, function(x) {
    x / sum(x)
  })
  atacPct <- apply(atacTable, 2, function(x) {
    x / sum(x)
  })
  rnaClusterName <- colnames(rnaPct)
  atacClusterName <- colnames(atacPct)
  calOvlpScoreElement <- function(t1l, t2l) {
    t1PctDF <- data.frame(rnaPct[, t1l])
    colnames(t1PctDF) <- "t1"
    t1PctDF$ident <- rownames(t1PctDF)
    t2PctDF <- data.frame(atacPct[, t2l])
    colnames(t2PctDF) <- "t2"
    t2PctDF$ident <- rownames(t2PctDF)
    comp <- plyr::join(t1PctDF, t2PctDF, by = "ident", type = "full")
    comp[is.na(comp)] <- 0
    comp$ident <- NULL
    comp <- t(comp)
    return(sum(apply(comp, 2, min)))
  }
  rna2atacType <- outer(rnaClusterName, atacClusterName, FUN = paste, sep = "|")
  ovlpScore <- apply(rna2atacType, MARGIN = c(1, 2), FUN = function(i) {
    t <- unlist(strsplit(i, split = "|", fixed = TRUE))
    calOvlpScoreElement(t[1], t[2])
  })
  rownames(ovlpScore) <- rnaClusterName
  colnames(ovlpScore) <- atacClusterName
  return(t(ovlpScore))
}
