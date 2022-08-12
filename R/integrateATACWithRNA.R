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
#' @param reduction characters, method used in FindTransferAnchors,
#' @param nVariableFeature integer, default is 3000
#' default "cca".
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
                                  reso = 0.5,
                                  reduction = "cca",
                                  nVariableFeature = 3000) {
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
    rnaSeurat <- FindVariableFeatures(rnaSeurat, nfeatures = nVariableFeature)
    rnaSeurat <- ScaleData(rnaSeurat, features = rownames(rnaSeurat))
    rnaSeurat <- RunPCA(rnaSeurat, features = VariableFeatures(rnaSeurat))
  }
  anchors <- FindTransferAnchors(
    reference = rnaSeurat,
    query = snapSeurat,
    features = VariableFeatures(rnaSeurat),
    reference.assay = "RNA",
    query.assay = snapAssay,
    reduction = reduction
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
  rna2atacType <- outer(rnaClusterName, atacClusterName,
                        FUN = paste, sep = "|")
  ovlpScore <- apply(rna2atacType, MARGIN = c(1, 2), FUN = function(i) {
    t <- unlist(strsplit(i, split = "|", fixed = TRUE))
    calOvlpScoreElement(t[1], t[2])
  })
  rownames(ovlpScore) <- rnaClusterName
  colnames(ovlpScore) <- atacClusterName
  return(t(ovlpScore))
}

#' Calculate overlap scores.
#' Old function, not expored and deprecated
calOvlpScoreOld <- function(t1, t2) {
  t1.table <- table(t1)
  t2.table <- table(t2)
  t1.pct <- apply(t1.table, 2, function(x) {
    x / sum(x)
  })
  t2.pct <- apply(t2.table, 2, function(x) {
    x / sum(x)
  })
  t1.labels <- colnames(t1.pct)
  t2.labels <- colnames(t2.pct)
  ovlpScore.df <- data.frame(
    anno1 = as.character(),
    anno2 = as.character(), ovlpScore = as.numeric()
  )
  for (t1.label in t1.labels) {
    for (t2.label in t2.labels) {
      t1.pct.df <- data.frame(t1.pct[, t1.label])
      colnames(t1.pct.df) <- "t1"
      t1.pct.df$ident <- rownames(t1.pct.df)
      t2.pct.df <- data.frame(t2.pct[, t2.label])
      colnames(t2.pct.df) <- "t2"
      t2.pct.df$ident <- rownames(t2.pct.df)
      comp.df <- plyr::join(t1.pct.df, t2.pct.df, by = "ident",
                      type = "full")
      comp.df[is.na(comp.df)] <- 0
      comp.df$ident <- NULL
      comp.df <- t(comp.df)
      ovlpScore <- sum(apply(comp.df, 2, min))
      out <- data.frame(anno1 = t1.label, anno2 = t2.label,
                        ovlpScore = ovlpScore)
      ovlpScore.df <- rbind(ovlpScore.df, out)
    }
  }
  return(ovlpScore.df)
}



#' Calculate Overlap Score matrix for ATAC and RNA jointly embedding
#' metaTable. [Deprecated]
#' @description
#' Based on the co-embeding clustering result, we summation the minmum scores of
#' the percentages of cells on either ATAC or RNA group in each cluster.
#' This function is provided by Yang Li with minor modification,
#' and now is deprecated.
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
getOverlapMatrixOld <- function(meta,
                                ident = "coembed.idents",
                                atacCol = "MajorType",
                                rnaCol = "ClusterName") {
  ident2rna <- data.frame(idents = meta[[ident]],
                          rna_label = meta[[rnaCol]])
  ident2rna <- ident2rna[complete.cases(ident2rna), ]
  ident2atac <- data.frame(idents = meta[[ident]], atac_label = meta[[atacCol]])
  ident2atac <- ident2atac[complete.cases(ident2atac), ]
  ovlp <- calOvlpScoreOld(ident2rna, ident2atac)
  colnames(ovlp) <- c(rnaCol, atacCol, "ovlpScore")
  ovlp <- reshape2::dcast(ovlp,
                          as.formula(paste(atacCol, "~", rnaCol, sep = " ")),
                          value.var = "ovlpScore",
                          fun.aggregate =  identity, fill = 0.0)
  ovlp <- as.data.frame(ovlp)
  rnm <- ovlp[[atacCol]]
  rownames(ovlp) <- rnm
  ## remove a name vector
  ovlp[[atacCol]] <- NULL
  ovlp <- as.matrix(sapply(ovlp, as.numeric))
  rownames(ovlp) <- rnm
  return(ovlp)
}
